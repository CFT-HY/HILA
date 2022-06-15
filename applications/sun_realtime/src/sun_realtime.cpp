#include "hila.h"
#include "gauge/staples.h"
#include "gauge/degauss.h"

#ifndef NSU
    #error "!!! Specify which SU(N) in Makefile: eg. -DNSU=3"
#endif

using SUN = SUmatrix<NSU, double>;

// Output stream for results
std::ofstream measureFile;

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


/* Calculate the lattice counterpart of B_i(x) = -1/2 eps_{ijk} F_jk(x) = -sum_{j<k} eps_{ijk} F_jk(x) everywhere. */
/* We extract the lattice field strength tensor from minimal clover terms:
* Q_{jk} = P_{jk} + P_{k,-j} + P_{-j,-k} + P_{-k,j}, then
* i g F_{jk} = (Q_{jk} - Q_{kj}) / (8a^2). 
* Because Q_{kj} = Q^+_{jk}, we only need Q_{jk} for j<k. 
* Here I do NOT project the result to algebra in order to avoid unnecessary operations
* => B here is of the GaugeField type. 
* Note that I do not subtract the trace from F_{jk}, so it's not necessarily traceless. 
* But the trace part drops out from the topological charge anyway since Tr E_i = 0. */
template <typename group>
void get_magnetic_field(const GaugeField<group> &U, GaugeField<group> &B) {

    foralldir(d1) onsites(ALL) B[d1][X] = 0;

    // "i,j,k" with j<k only
    foralldir(d1) {

        foralldir(d2) foralldir(d3) if (d1 != d2 && d1 != d3 && d2 < d3) {
       
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
            
            // Get local magnetic field (actually g a^2 B_i, anti-Hermitian algebra)
            onsites(ALL) {
                B[d1][X] += -(1.0/ 8.0) * sign * (Q[X] - Q[X].dagger());
            }
            
        } // end foralldir d2 d3
        
    } // end d1; mag. field done
 
}


/* Measure classical topological charge density chi everywhere */ 
template <typename group>
void calc_topoCharge(const GaugeField<group> &U, const VectorField<Algebra<group>> &E, Field<double> &result) {

    double c_chi = 1.0 / (64.0 * M_PI*M_PI);
    onsites(ALL) result[X] = 0;

    // Magnetic field, but not projected to algebra for better performance:
    GaugeField<group> B;
    get_magnetic_field(U, B);

    /* Now a^4 chi = 16 c_chi a^4 g^2 Tr (E_i B_i) with antiHermitian E_i, B_i, 
    * and the lattice fields are actually g a^2 E_i and g a^2 B_i. */
    foralldir(d1) onsites(ALL) {
        result[X] += 16.0 * c_chi * real( trace(E[d1][X].expand() * B[d1][X]) );
        // this trace is automatically real
    } 

}


static double degauss_quality = 1e-12;

// Do the measurements. Here 't' labels the Hamiltonian time 
template <typename group>
void measure_stuff(GaugeField<group> &U, VectorField<Algebra<group>> &E, int trajectory, double t) {

    auto plaq = measure_plaq(U); // N - Tr Re P_ij

    double e2 = 0.0;
    foralldir (d)
        e2 += E[d].squarenorm();
    e2 /= 2; // now e2 = Tr E_i^2

    // total energy ("action") times g^2 a T
    double energy = e2 + 2.0 * plaq;

    Field<Algebra<group>> g;
    get_gauss_violation(U, E, g);
    auto viol = g.squarenorm();

    /* Measure 'improved' or 'symmetrized' charge density: The way E_i appears in the EOM suggests
    * that we should identify E_i(x) as living at link midpoint, while the mag. field B_i(x) is local to x.
    * To measure E.B at x, we take covariant average of E: E^{imp}_i(x) = 0.5 * (E_i(x) + U^+_i(x-i)E_i(x-i)U_i(x-i)) */
    Field<double> chi;
    VectorField<Algebra<group>> E_imp;
    foralldir(d1) onsites(ALL) { 
        E_imp[d1][X] = 0.5 * (E[d1][X] + 
            (U[d1][X-d1].dagger() * E[d1][X-d1].expand() * U[d1][X-d1]).project_to_algebra() );
    }
    calc_topoCharge(U, E_imp, chi);

    double chi_avg = 0.0;
    onsites(ALL) 
         chi_avg += chi[X];

    chi_avg /= lattice->volume();


    char buf[1024]; 
    sprintf(buf, "%d %.10g %.8g %.8g %.8g %.8g %.8g", trajectory, t, plaq, e2, viol, energy, chi_avg);
    if (hila::myrank() == 0) measureFile << std::string(buf) << "\n";
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

        /*
        double pl = measure_plaq(U);
        double e2 = 0;
        foralldir (d) { e2 += E[d].squarenorm(); }

        output0 << "THERM: Plaq: " << pl << " E^2 " << e2 << " action "
                << e2 / 2 + 2 * pl << '\n';
        */

        regroup_gauge(U);
    }
    therm_timer.stop();
}

////////////////////////////////////////////////////////////////

template <typename group>
void do_trajectory(GaugeField<group> &U, VectorField<Algebra<group>> &E, int trajectory,
                   int trajlen, int measure_interval, double dt) {

    // Measure first at time t=0:
    double t = 0.0;
    measure_stuff(U, E, trajectory, t);

    // Then evolve until we reach t = trajlen*dt
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
        
        t += dt * measure_interval;
        measure_stuff(U, E, trajectory, t);
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
    std::string meas_fname = par.get("measurement file");

    par.close(); // file is closed also when par goes out of scope

    // setting up the lattice is convenient to do after reading
    // the parameter
    lattice->setup(lsize);

    // We need random number here
    hila::seed_random(seed);


    if (hila::myrank() == 0) {
        measureFile.open(meas_fname, std::ios_base::app);
        if (!measureFile) {
            output0 << "!!! Error opening measurement file\n";
            hila::finishrun();
        }
    }

    // Print measurement labels
    if (hila::myrank() == 0) {
        std::ofstream labelFile;
        std::string labelFileName = "labels_" + meas_fname;
        labelFile.open(labelFileName);
        if (!labelFile) {
            output0 << "!!! Error opening file " << labelFileName << "\n";
            hila::finishrun();
        }
        labelFile << "1 trajectory\n" << "2 time (lattice units)\n" << "3 plaq avg: sum_i<j (N - Tr Re P_ij)\n" 
            << "4 Tr E_i^2\n" << "5 Gauss violation (G^a G^a over the whole system)\n" << "6 total energy\n"
            << "7 average chi\n";

        labelFile.close();
    }

    // Alloc gauge field and momenta (E)
    GaugeField<SUN> U;
    VectorField<Algebra<SUN>> E;
    foralldir(d) onsites(ALL) {
        U[d][X] = 0;
        E[d][X] = 0;
    }

    // some initial noise for gauge field
    foralldir (d) {
        onsites (ALL) {
            U[d][X].gaussian_random(0.3).reunitarize();
        }
    }
     

    thermalize(U, E, g2Ta, n_thermal_start, dt);

    for (int trajectory = 0; trajectory < n_traj; trajectory++) {
        thermalize(U, E, g2Ta, n_thermal, dt);
        if (trajectory % 500 == 0) {
            output0 << "Trajectory " << trajectory << "\n";
        }
        do_trajectory(U, E, trajectory, trajlen, measure_interval, dt);
    }


    // done
    if (hila::myrank() == 0) {
        measureFile.close();
    }

    hila::finishrun();
}
