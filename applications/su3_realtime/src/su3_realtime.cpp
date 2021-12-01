#include "hila.h"
#include "gauge/staples.h"
#include "gauge/degauss.h"

using SU3 = SUmatrix<3, double>;


template <typename group>
void update_E(const GaugeField<group> &U, VectorField<Algebra<group>> &E,
              double delta) {

    Field<group> staple;

    foralldir (d) {
        staplesum(U, staple, d);

        onsites (ALL) {
            E[d][X] -= delta * (U[d][X] * staple[X].dagger()).project_to_algebra();
        }
    }
}

template <typename group>
void update_U(GaugeField<group> &U, const VectorField<Algebra<group>> &E,
              double delta) {

    foralldir (d) {
        onsites (ALL)
            U[d][X] = exp(E[d][X] * delta) * U[d][X];
    }
}

template <typename group>
void regroup_gauge(GaugeField<group> &U) {
    foralldir (d) {
        onsites (ALL)
            U[d][X].reunitarize();
    }
}

template <typename group>
double measure_plaq(const GaugeField<group> &U) {

    Reduction<double> plaq;
    plaq.allreduce(false).delayed(true);


    foralldir (dir1)
        foralldir (dir2)
            if (dir1 < dir2) {
                onsites (ALL) {
                    plaq += group::size() - real(trace(U[dir1][X] * U[dir2][X + dir1] *
                                                       U[dir1][X + dir2].dagger() *
                                                       U[dir2][X].dagger()));
                }
            }
    return plaq.value();
}


static double degauss_quality = 1e-12;

template <typename group>
void measure_stuff(GaugeField<group> &U, VectorField<Algebra<group>> &E, int trajectory,
                   int n, double dt) {

    auto plaq = measure_plaq(U);

    double e2;
    foralldir (d)
        e2 += E[d].squarenorm();
    e2 /= 2;

    Field<Algebra<group>> g;
    get_gauss_violation(U, E, g);
    auto viol = g.squarenorm();

    //output0 << "Measure_start " << n << "\n";
    output0 << "MEAS " << trajectory << ' ' << n << ' ' << plaq << ' ' << e2 << ' '
            << viol << '\n';
    //output0 << "Measure_end " << n << "\n";
}


template <typename group>
void thermalize(GaugeField<group> &U, VectorField<Algebra<group>> &E, double g2Ta,
                int iterations, double dt) {

    regroup_gauge(U);

    static hila::timer therm_timer("Thermalization");

    therm_timer.start();
    for (int loop = 0; loop < iterations; loop++) {
        foralldir (d)
            onsites (ALL)
                E[d][X].gaussian_random(sqrt(g2Ta));

        degauss(U, E, degauss_quality);

        // 1/2 timestep for U first
        update_U(U, E, dt / 2);
        // do 1 time unit of evolution with leapfrog
        for (int steps = 0; steps < 1.0 / dt - 1; steps++) {
            update_E(U, E, dt);
            update_U(U, E, dt);
        }
        // and bring U and E to the same time value
        update_E(U, E, dt);
        update_U(U, E, dt / 2);

        double pl = measure_plaq(U);
        double e2 = 0;
        foralldir (d) { e2 += E[d].squarenorm(); }

        output0 << "THERM: Plaq: " << pl << " E^2 " << e2 << " action "
                << e2 / 2 + 2 * pl << '\n';

        regroup_gauge(U);
    }
    therm_timer.stop();
}

////////////////////////////////////////////////////////////////

template <typename group>
void do_trajectory(GaugeField<group> &U, VectorField<Algebra<group>> &E, int trajectory,
                   int trajlen, int measure_interval, double dt) {


    for (int n = 0; n < trajlen; n += measure_interval) {
        update_U(U, E, dt / 2);
        // do 2 time units of evolution with leapfrog
        for (int steps = 0; steps < measure_interval - 1; steps++) {
            update_E(U, E, dt);
            update_U(U, E, dt);
        }
        // and bring U and E to the same time value
        update_E(U, E, dt);
        update_U(U, E, dt / 2);

        measure_stuff(U, E, trajectory, n, dt);
    }
}

/////////////////////////////////////////////////////////////////

   #pragma hila ast_dump
    template <typename T>
    class testc {
        using base_type = hila::number_type<T>;
        using argument_type = T;
        double v;
    };



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

    double g2Ta = par.get("g^2 Ta");
    double dt = par.get("dt");
    int trajlen = par.get("trajectory length");
    int n_traj = par.get("number of trajectories");
    int measure_interval = par.get("measurement interval");
    int n_thermal = par.get("thermalisation");
    int n_thermal_start = par.get("thermalisation start");
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

    // some initial noise for gauge field
    foralldir (d) {
        onsites (ALL)
            U[d][X].gaussian_random(0.3).reunitarize();
    }

    testc<double> td;

    Field<Complex<double>> cd(1);
    Complex<double> sum(0);
    onsites(ALL) sum += cd[X];
     

    thermalize(U, E, g2Ta, n_thermal_start, dt);

    for (int trajectory = 0; trajectory < n_traj; trajectory++) {
        thermalize(U, E, g2Ta, n_thermal, dt);
        do_trajectory(U, E, trajectory, trajlen, measure_interval, dt);
    }


    hila::finishrun();
}
