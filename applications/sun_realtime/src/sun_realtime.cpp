#include "hila.h"
#include "gauge/staples.h"
#include "gauge/degauss.h"

#ifndef NSU
    #error "!!! Specify which SU(N) in Makefile: eg. -DNSU=3"
#endif

using SUN = SUmatrix<NSU, double>;


///////////////////////////////////////////////////////////////////////////////////////////////
/* Hamiltonian time evolution for gauge and electric fields. 'delta' is the time difference. */
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

///////////////////////////////////////////////////////////////////////////////////////////////

template <typename group>
void regroup_gauge(GaugeField<group> &U) {
    foralldir (d) {
        onsites (ALL)
            U[d][X].reunitarize();
    }
}

// Total plaquette action (without beta prefactor)
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


/* Calculate the lattice counterpart of B_i(x) = -1/2 eps_{ijk} F_jk(x) everywhere. */
/* We extract the lattice field strength tensor from minimal clover terms:
* Q_{jk} = P_{jk} + P_{k,-j} + P_{-j,-k} + P_{-k,j}, then
* i g F_{jk} = (Q_{jk} - Q_{kj}) / (8a^2). 
* Because of antisymmetry we only need Q_{jk} for j<k. 
* Here I do NOT project the result to algebra in order to avoid unnecessary operations
* => B here is of the GaugeField type. 
* Note that I do not subtract the trace from F_{jk}, so it's not necessarily traceless. 
* But the trace part drops out from the topological charge anyway since Tr E_i = 0. */
template <typename group>
void get_magnetic_field(const GaugeField<group> &U, GaugeField<group> &B) {

    foralldir(d1) onsites(ALL) B[d1][X] = 0;

    // "i,j,k"
    foralldir(d1) foralldir(d2) foralldir(d3) if (d1 != d2 && d1 != d3 && d2 < d3) {
       
        /* Get clovers in the plane orthogonal to 'd1': */
        Field<group> Q;
        onsites(ALL) {
            group P1, P2, P3, P4;
            P1 = U[d2][X] * U[d3][X + d2] * U[d2][X + d3].dagger() * U[d3][X].dagger();
            P2 = U[d3][X] * U[d2][X + d3 - d2].dagger() * U[d3][X - d2].dagger() * U[d2][X - d2];
            P3 = U[d2][X - d2].dagger() * U[d3][X - d2 - d3].dagger() * U[d2][X - d3 - d2] * U[d3][X - d3]; 
            P4 = U[d3][X - d3].dagger() * U[d2][X - d3] * U[d3][X + d2 - d3] * U[d2][X].dagger();

            Q[X] = P1 + P2 + P3 + P4;
        }

        /* B_i lives on the links, so improve the definition by calculating F_{jk} at the midpoint:
        * F_{jk}(x + 0.5i) = 1/2 (F_{jk}(x) + U_i(x) F_{jk}(x+i) U_i(x).dagger() ),
        * where the second term is a parallel transport to the next lattice site. 
        * Then i g a^2 B_i(x + 0.5i) = -1/8 Sum_{j<k} e_{ijk} ( Q_{jk}(x) + U_i(x) Q_{jk}(x+i) U_i(x).dagger() ).
        * The code uses anti-Hermitian generators for the algebra, so we actually compute g a^2 B_i. */
        
        // int sign = (d1-d2)*(d2-d3)*(d3-d1)/2; // epsilon tensor in 3 dimensions
        int sign;
        // Here d2 < d3 always
        if (d2 > d1) {
            sign = 1;
        } else if (d3 > d1) {
            sign = -1;
        } else {
            sign = 1;
        }

        onsites(ALL) {
            B[d1][X] += -(1.0/8.0) * sign * (Q[X] + U[d1][X] * Q[X+d1] * U[d1][X].dagger());
            // Without averaging over the two clovers:
            // B[d1][X] += (-2 * (1.0/8.0) * sign * (Q[X]));
        }
    } // end foralldir d1 d2 d3
 
}


/* Measure classical topological charge density chi everywhere */ 
template <typename group>
void calc_topoCharge(const GaugeField<group> &U, const VectorField<Algebra<group>> &E, Field<double> &result) {

    double c_chi = 1.0 / (64.0 * M_PI*M_PI);
    result[ALL] = 0;

    // Magnetic field, but not projected to algebra for better performance:
    GaugeField<group> B;
    get_magnetic_field(U, B);

    /* Now a^4 chi = -16 c_chi a^4 g^2 Tr (E_i B_i),
    * and the lattice fields are actually g a^2 E_i and g a^2 B_i. */
    foralldir(d1) onsites(ALL) {
        result[X] -= 16.0 * c_chi * real( trace(E[d1][X].expand() * B[d1][X]) );
    }

}


static double degauss_quality = 1e-12;

// Do the measurements. Here 'n' just labels the time step, so the actual time is n*dt. 
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

    Field<double> chi;
    calc_topoCharge(U, E, chi);

    Reduction<double> chi_total;
    chi_total.allreduce(true).delayed(true);
    onsites(ALL) 
         chi_total += chi[X];

    double chi_avg = chi_total.value() / lattice->volume();

    //output0 << "Measure_start " << n << "\n";
    // Print the actual time (in lattice units) instead of just 'n'. Also more precision, needed for long trajectories
    char buf[1024]; sprintf(buf, "MEAS %d %.10g %.8g %.8g %.8g %.8g", trajectory, n*dt, plaq, e2, viol, chi_avg);
    output0 << std::string(buf) << "\n";

    // output0 << "MEAS " << trajectory << ' ' << n*dt << ' ' << plaq << ' ' << e2 << ' ' << viol << ' ' << chi_avg << '\n';
    //output0 << "Measure_end " << n << "\n";
}


/* Thermalization according to ref. hep-ph/9603384 (see Appendix A).
* Gauss constraint is enforced by degauss() routine 
* that operates as described around eq. (38) in the reference. */
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



int main(int argc, char **argv) {

    // hila::initialize should be called as early as possible
    hila::initialize(argc, argv);

    output0 << "---- Hamiltonian time evolution for SU(N) theory, using N=" << NSU <<  " ----\n";

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
    GaugeField<SUN> U;
    VectorField<Algebra<SUN>> E;

    // some initial noise for gauge field
    foralldir (d) {
        onsites (ALL) {
            U[d][X].gaussian_random(0.3).reunitarize();
        }
    }
     

    thermalize(U, E, g2Ta, n_thermal_start, dt);

    for (int trajectory = 0; trajectory < n_traj; trajectory++) {
        thermalize(U, E, g2Ta, n_thermal, dt);
        do_trajectory(U, E, trajectory, trajlen, measure_interval, dt);
    }


    hila::finishrun();
}
