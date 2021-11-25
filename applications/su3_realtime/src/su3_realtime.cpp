#include "hila.h"
#include "gauge/staples.h"
#include "gauge/degauss.h"

using SU3 = SUmatrix<3, double>;


template <typename group>
void update_E(const GaugeField<group> &U, VectorField<Algebra<group>> &E,
              double delta) {

    Field<group> staple;

    foralldir(d) {
        staplesum(U, staple, d);

        onsites(ALL) {
            E[d][X] -= delta * (U[d][X] * staple[X].dagger()).project_to_algebra();
        }
    }
}

template <typename group>
void update_U(GaugeField<group> &U, const VectorField<Algebra<group>> &E,
              double delta) {

    foralldir(d) {
        onsites(ALL) U[d][X] = exp(E[d][X] * delta) * U[d][X];
    }
}

template <typename group>
void regroup_gauge(GaugeField<group> &U) {
    foralldir(d) {
        onsites(ALL) U[d][X].reunitarize();
    }
}

template <typename group>
double measure_plaq(const GaugeField<group> &U) {
    Reduction<double> plaq;
    plaq = 0;
    plaq.allreduce(false).delayed(true);

    foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {
        onsites(ALL) {
            plaq += real(trace(U[dir1][X] * U[dir2][X + dir1] *
                               U[dir1][X + dir2].dagger() * U[dir2][X].dagger()));
        }
    }
    return plaq.value();
}


static double degauss_quality = 1e-12;


template <typename group>
void thermalize(GaugeField<group> &U, VectorField<Algebra<group>> &E, int iterations,
                double dt) {

    regroup_gauge(U);

    static hila::timer therm_timer("Thermalization");

    therm_timer.start();
    for (int loop = 0; loop < iterations; loop++) {
        foralldir(d) onsites(ALL) E[d][X].gaussian_random(0.5);

        degauss(U, E, degauss_quality);

            auto pl = measure_plaq(U);
            double e2 = 0;
            foralldir(d) {
                e2 += E[d].squarenorm();
            }

            output0 << " Before: Plaq 1: " << pl << " E^2 " << e2 << " action " << e2/4 - pl 
                    << '\n';

        // 1/2 timestep for U first
        update_U(U, E, dt / 2);
        // do 2 time units of evolution with leapfrog
        for (int steps = 0; steps < 1.0 / dt - 1; steps++) {
            update_E(U, E, dt);
            update_U(U, E, dt);
        }
        // and bring U and E to the same time value
        update_E(U, E, dt);
        update_U(U, E, dt / 2);

         pl = measure_plaq(U);
        e2 = 0;
            foralldir(d) {
                e2 += E[d].squarenorm();
            }

            output0 << " After: Plaq 1: " << pl << " E^2 " << e2 << " action " << e2/4 - pl 
                    << '\n';

        regroup_gauge(U);

    }
    therm_timer.stop();
}


int main(int argc, char **argv) {

    // hila::initialize should be called as early as possible
    hila::initialize(argc, argv);

    // hila provides an input class hila::input, which is
    // a convenient way to read in parameters from input files.
    // parameters are presented as key - value pairs, as an example
    //  " lattice size  64, 64, 64 "
    // is read below.
    //
    // Values are broadcast to all MPI nodes.
    //
    // .get() -method can read many different input types,
    // see file "input.h" for documentation

    hila::input par("parameters");

    CoordinateVector lsize;
    lsize = par.get("lattice size"); // reads NDIM numbers

    int loops = par.get("updates");
    long seed = par.get("random seed");

    par.close(); // file is closed also when par goes out of scope

    // setting up the lattice is convenient to do after reading
    // the parameter
    lattice->setup(lsize);

    // We need random number here
    hila::seed_random(seed);

    // Alloc gauge field and momenta (E)
    GaugeField<SU3> U;
    VectorField<Algebra<SU3>> E;

    foralldir(d) {
        onsites(ALL) U[d][X].gaussian_random(0.3).reunitarize();
    }

    thermalize(U, E, 20, 0.02);

    double n = 0;
    foralldir(d) n += E[d].squarenorm();

    output0 << " E square norm " << n / (NDIM * lattice->volume()) << '\n';


    hila::finishrun();
}
