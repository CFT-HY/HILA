#include "hila.h"
#include "gauge/staples.h"
#include "gauge/sun_heatbath.h"
#include "gauge/sun_overrelax.h"
#include "gauge/gradient_flow.h"
#include "tools/string_format.h"
#include "tools/floating_point_epsilon.h"

#ifndef NCOLOR
#define NCOLOR 3
#endif

#ifndef BCOPEN
#define BCOPEN -1
#endif

using ftype = double;
using mygroup = SU<NCOLOR, ftype>;

// define a struct to hold the input parameters: this
// makes it simpler to pass the values around
struct parameters {
    ftype beta;         // inverse gauge coupling
    int n_traj;         // number of trajectories to generate
    int n_therm;        // number of thermalization trajectories (counts only accepted traj.)
    int n_update;       // number of heat-bath sweeps per "trajectory"
    int n_overrelax;    // number of overrelaxation sweeps between heat-batch sweeps
    int gflow_freq;     // number of trajectories between gflow measurements
    ftype gflow_max_l;  // flow scale at which gradient flow stops
    ftype gflow_l_step; // flow scale interval between flow measurements
    ftype gflow_a_accu; // desired absolute accuracy of gradient flow integration steps
    ftype gflow_r_accu; // desired relative accuracy of gradient flow integration steps
    int n_save;         // number of trajectories between config. check point
    std::string config_file;
    ftype time_offset;
};


///////////////////////////////////////////////////////////////////////////////////
// heat-bath functions

/**
 * @brief Wrapper update function
 * @details Gauge Field update sweep with randomly chosen parities and directions
 *
 * @tparam group
 * @param U GaugeField to update
 * @param p Parameter struct
 * @param relax If true evolves GaugeField with over relaxation if false then with heat bath
 * @param plaqw plaquette weights
 */
template <typename group>
void update(GaugeField<group> &U, const parameters &p, bool relax, const plaqw_t<ftype> &plaqw) {
    for (int i = 0; i < 2 * NDIM; ++i) {
        int tdp = hila::broadcast((int)(hila::random() * 2 * NDIM));
        int tdir = tdp / 2;
        int tpar = 1 + (tdp % 2);
        //hila::out0 << "   " << Parity(tpar) << " -- " << Direction(tdir);
        update_parity_dir(U, p, Parity(tpar), Direction(tdir), relax, plaqw);
    }
    //hila::out0 << "\n";
}

/**
 * @brief Wrapper function to updated GaugeField per direction
 * @details Computes first staplesum, then uses computed result to evolve GaugeField either with
 * over relaxation or heat bath
 *
 * @tparam group
 * @param U GaugeField to evolve
 * @param p parameter struct
 * @param par Parity
 * @param d Direction to evolve
 * @param relax If true evolves GaugeField with over relaxation if false then with heat bath
 * @param plaqw plaquette weights
 */
template <typename group>
void update_parity_dir(GaugeField<group> &U, const parameters &p, Parity par, Direction d,
                       bool relax, const plaqw_t<ftype> &plaqw) {

    static hila::timer hb_timer("Heatbath");
    static hila::timer or_timer("Overrelax");
    static hila::timer staples_timer("Staplesum");

    Field<group> staples;

    staples_timer.start();

    staplesum(U, staples, d, plaqw, par);

    staples_timer.stop();

    if (relax) {

        or_timer.start();

        onsites(par) {
            if(plaqw[d][d][X] != 0) {
#ifdef SUN_OVERRELAX_dFJ
                suN_overrelax_dFJ(U[d][X], staples[X], p.beta);
#else
                suN_overrelax(U[d][X], staples[X]);
#endif
            }
        }
        or_timer.stop();

    } else {

        hb_timer.start();
        onsites(par) {
            if (plaqw[d][d][X] != 0) {
                suN_heatbath(U[d][X], staples[X], p.beta);
            }
        }
        hb_timer.stop();
    }
}

/**
 * @brief Evolve gauge field
 * @details Evolution happens by means of heat bath and over relaxation. For each heatbath update
 * (p.n_update) we update on average also p.n_overrelax times with over relaxation.
 *
 * @tparam group
 * @param U
 * @param p
 */
template <typename group>
void do_trajectory(GaugeField<group> &U, const plaqw_t<ftype> &plaqw, const parameters &p) {

    for (int n = 0; n < p.n_update; n++) {
        for (int i = 0; i <= p.n_overrelax; i++) {
            bool relax = hila::broadcast((int)(hila::random() * (1 + p.n_overrelax)) != 0);
            update(U, p, relax, plaqw);
            //hila::out0 << relax << "\n";
        }
    }
    U.reunitarize_gauge();
}

// heat-bath functions
///////////////////////////////////////////////////////////////////////////////////
// measurement functions

template <typename group>
std::vector<double> measure_s_wplaq_ps(const GaugeField<group> &U, Direction obd,
                                       const plaqw_t<ftype> &plaqw) {
    // measure the total Wilson plaquette action for the gauge field U
    int obdlsize = lattice.size(obd);
    ReductionVector<double> plaql(obdlsize);
    plaql = 0.0;
    plaql.allreduce(false).delayed(true);
    double normf = 2.0 / (NDIM * (NDIM - 1) * lattice.volume() / obdlsize);

    foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {
        U[dir2].start_gather(dir1, ALL);
        U[dir1].start_gather(dir2, ALL);
        onsites(ALL) {
            double tplq = (1.0 - real(trace(U[dir1][X] * U[dir2][X + dir1] *
                                            (U[dir2][X] * U[dir1][X + dir2]).dagger())) /
                                     group::size()) *
                          plaqw[dir1][dir2][X] * normf;
            int obdind = X.coordinate(obd);
            if (dir1 != obd && dir2 != obd) {
                plaql[obdind] += tplq;
            } else {
                if (obdind > 0 && obdind < obdlsize - 1) {
                    if (obdind == 1) {
                        plaql[obdind] += tplq;
                        plaql[obdind + 1] += 0.5 * tplq;
                    } else if (obdind == obdlsize - 2) {
                        plaql[obdind] += 0.5 * tplq;
                        plaql[obdind + 1] += tplq;
                    } else {
                        plaql[obdind] += 0.5 * tplq;
                        plaql[obdind + 1] += 0.5 * tplq;
                    }
                }
            }
        }
    }
    plaql.reduce();
    return plaql.vector();
}

/**
 * @brief Measure Polyakov lines to direction dir
 * @details Naive implementation, includes extra communication
 * @tparam T GaugeField Group
 * @param U GaugeField to measure
 * @param dir Direction
 * @return Complex<double>
 */
template <typename T>
std::vector<Complex<double>> measure_polyakov_ps(const GaugeField<T> &U, Direction obd,
                                                 Direction dir = Direction(NDIM - 1)) {

    Field<T> polyakov = U[dir];

    // mult links so that polyakov[X.dir == 0] contains the polyakov loop
    for (int plane = lattice.size(dir) - 2; plane >= 0; plane--) {

        // safe_access(polyakov) pragma allows the expression below, otherwise
        // hilapp would reject it because X and X+dir can refer to the same
        // site on different "iterations" of the loop.  However, here this
        // is restricted on single dir-plane so it works but we must tell it to hilapp.

#pragma hila safe_access(polyakov)
        onsites(ALL) {
            if (X.coordinate(dir) == plane) {
                polyakov[X] = U[dir][X] * polyakov[X + dir];
            }
        }
    }

    ReductionVector<Complex<double>> ploopl(lattice.size(obd));
    ploopl = 0.0;
    ploopl.allreduce(false).delayed(true);

    int obdlsize = lattice.size(obd);
    double normf = 1.0 / (lattice.volume() / (lattice.size(dir) * obdlsize));

    onsites(ALL) if (X.coordinate(dir) == 0) {
        ploopl[X.coordinate(obd)] += trace(polyakov[X]) * normf;
    }
    ploopl.reduce();

    return ploopl.vector();
}

template <typename group>
void measure_stuff(const GaugeField<group> &U, const plaqw_t<ftype> &plaqw) {
    // perform measurements on current gauge and momentum pair (U, E) and
    // print results in formatted form to standard output
    static bool first = true;
#if BCOPEN >= 0 && BCOPEN < NDIM
    if (first) {
        // print legend for measurement output
        hila::out0 << "LMPLAQPS:     plaq[1]       plaq[2]       plaq[3]       plaq[4]          ...\n";
        hila::out0 << "LMPOLPS:      P[1].re       P[1].im       P[2].re       P[2].im          ...\n";
        first = false;
    }
    Direction obd = Direction(BCOPEN);
    auto plaql = measure_s_wplaq_ps(U, obd, plaqw);
    auto polyl = measure_polyakov_ps(U, obd, e_t);
    hila::out0 << "MPLAQPS";
    for (int i = 1; i < lattice.size(obd); ++i) {
        hila::out0 << string_format(" % 0.6e", plaql[i]);
    }
    hila::out0 << '\n';
    hila::out0 << "MPOLPS ";
    for (int i = 1; i < lattice.size(obd); ++i) {
        hila::out0 << string_format(" % 0.6e % 0.6e", polyl[i].real(), polyl[i].imag());
    }
    hila::out0 << '\n';
#else
    if (first) {
        // print legend for measurement output
        hila::out0 << "LMEAS:        plaq        P.real        P.imag\n";
        first = false;
    }
    auto plaq = measure_s_wplaq(U) / (lattice.volume() * NDIM * (NDIM - 1) / 2);
    auto poly = measure_polyakov(U, e_t);
    hila::out0 << string_format("MEAS % 0.6e % 0.6e % 0.6e", plaq, poly.real(), poly.imag())
            << '\n';
#endif
}

// end measurement functions
///////////////////////////////////////////////////////////////////////////////////
// load/save config functions

template <typename group>
void checkpoint(const GaugeField<group> &U, int trajectory, const parameters &p) {
    double t = hila::gettime();
    // name of config with extra suffix
    std::string config_file =
        p.config_file + "_" + std::to_string(abs((trajectory + 1) / p.n_save) % 2);
    // save config
    U.config_write(config_file);
    // write run_status file
    if (hila::myrank() == 0) {
        std::ofstream outf;
        outf.open("run_status", std::ios::out | std::ios::trunc);
        outf << "trajectory  " << trajectory + 1 << '\n';
        outf << "seed        " << static_cast<uint64_t>(hila::random() * (1UL << 61)) << '\n';
        outf << "time        " << hila::gettime() << '\n';
        // write config name to status file:
        outf << "config name  " << config_file << '\n';
        outf.close();
    }
    std::stringstream msg;
    msg << "Checkpointing, time " << hila::gettime() - t;
    hila::timestamp(msg.str().c_str());
}

template <typename group>
bool restore_checkpoint(GaugeField<group> &U, int &trajectory, parameters &p) {
    uint64_t seed;
    bool ok = true;
    p.time_offset = 0;
    hila::input status;
    if (status.open("run_status", false, false)) {
        hila::out0 << "RESTORING FROM CHECKPOINT:\n";
        trajectory = status.get("trajectory");
        seed = status.get("seed");
        p.time_offset = status.get("time");
        // get config name with suffix from status file:
        std::string config_file = status.get("config name");
        status.close();
        hila::seed_random(seed);
        U.config_read(config_file);
        ok = true;
    } else {
        std::ifstream in;
        in.open(p.config_file, std::ios::in | std::ios::binary);
        if (in.is_open()) {
            in.close();
            hila::out0 << "READING initial config\n";
            U.config_read(p.config_file);
            ok = true;
        } else {
            ok = false;
        }
    }
    return ok;
}

// end load/save config functions
///////////////////////////////////////////////////////////////////////////////////


int main(int argc, char **argv) {

    // hila::initialize should be called as early as possible
    hila::initialize(argc, argv);

    // hila provides an input class hila::input, which is
    // a convenient way to read in parameters from input files.
    // parameters are presented as key - value pairs, as an example
    //  " lattice size  64, 64, 64, 64"
    // is read below.
    //
    // Values are broadcast to all MPI nodes.
    //
    // .get() -method can read many different input types,
    // see file "input.h" for documentation

    hila::out0 << "SU(" << mygroup::size() << ") heat-bath with overrelaxation\n";

    hila::out0 << "Using floating point epsilon: " << fp<ftype>::epsilon << "\n";

    parameters p;

    hila::input par("parameters");

    CoordinateVector lsize;
    // reads NDIM numbers
    lsize = par.get("lattice size");
    // inverse gauge coupling
    p.beta = par.get("beta");
    // number of trajectories
    p.n_traj = par.get("number of trajectories");
    // number of heat-bath (HB) sweeps per trajectory
    p.n_update = par.get("updates in trajectory");
    // number of overrelaxation sweeps petween HB sweeps
    p.n_overrelax = par.get("overrelax steps");
    // number of thermalization trajectories
    p.n_therm = par.get("thermalization trajs");
    // wilson flow frequency (number of traj. between w. flow measurement)
    p.gflow_freq = par.get("gflow freq");
    // wilson flow max. flow-distance
    p.gflow_max_l = par.get("gflow max lambda");
    // wilson flow flow-distance step size
    p.gflow_l_step = par.get("gflow lambda step");
    // wilson flow absolute accuracy (per integration step)
    p.gflow_a_accu = par.get("gflow abs. accuracy");
    // wilson flow relative accuracy (per integration step)
    p.gflow_r_accu = par.get("gflow rel. accuracy");
    // random seed = 0 -> get seed from time
    long seed = par.get("random seed");
    // save config and checkpoint
    p.n_save = par.get("trajs/saved");
    // measure surface properties and print "profile"
    p.config_file = par.get("config name");

    par.close(); // file is closed also when par goes out of scope

    // set up the lattice
    lattice.setup(lsize);

    // We need random number here
    hila::seed_random(seed);

    // Alloc gauge field and momenta (E)
    GaugeField<mygroup> U;
    VectorField<Algebra<mygroup>> E;

    plaqw_t<ftype> plaqw;
    foralldir(d1) {
        foralldir(d2) {
            onsites(ALL) {
                plaqw[d1][d2][X] = 1.0;
            }
        }
    }

#if BCOPEN>=0 && BCOPEN<NDIM
    Direction obd = Direction(BCOPEN);
    hila::out0 << "Using open boundary conditions in " << obd << "-direction\n";
    foralldir(d1) {
        foralldir(d2) {
            onsites(ALL) {
                if (X.coordinate(obd) == 0) {
                    plaqw[d1][d2][X] = 0;
                }
                if (X.coordinate(obd) == 1) {
                    if (d1 != obd && d2 != obd) {
                        //plaqw[d1][d2][X] = 0.5;
                    }
                }
                if (X.coordinate(obd) == lattice.size(obd) - 1) {
                    if (d1 != obd && d2 != obd) {
                        //plaqw[d1][d2][X] = 0.5;
                    } else {
                        plaqw[d1][d2][X] = 0;
                    }
                }
            }
        }
    }
#endif //END BCOPEN

    // use negative trajectory for thermal
    int start_traj = -p.n_therm;

    if (!restore_checkpoint(U, start_traj, p)) {
        U = 1;
    }


    hila::timer update_timer("Updates");
    hila::timer measure_timer("Measurements");
    hila::timer gf_timer("Gradient Flow");

    ftype t_step0 = 0;
    for (int trajectory = start_traj; trajectory < p.n_traj; ++trajectory) {

        ftype ttime = hila::gettime();

        update_timer.start();

        do_trajectory(U, plaqw, p);

        // put sync here in order to get approx gpu timing
        hila::synchronize_threads();
        update_timer.stop();

        measure_timer.start();

        hila::out0 << "Measure_start " << trajectory << '\n';

        measure_stuff(U, plaqw);

        hila::out0 << "Measure_end " << trajectory << '\n';

        measure_timer.stop();

        if (trajectory >= 0) {

            if (p.gflow_freq > 0 && trajectory % p.gflow_freq == 0) {


                int gtrajectory = trajectory / p.gflow_freq;
                if (p.gflow_l_step > 0) {

                    gf_timer.start();

                    int nflow_steps = (int)(p.gflow_max_l / p.gflow_l_step);
                    ftype gftime = hila::gettime();
                    hila::out0 << "Gflow_start " << gtrajectory << '\n';

                    GaugeField<mygroup> V = U;

                    ftype t_step = t_step0;
                    measure_gradient_flow_stuff(V, (ftype)0.0, t_step);
                    t_step = do_gradient_flow_adapt(V, (ftype)0.0, p.gflow_l_step, p.gflow_a_accu,
                                                    p.gflow_r_accu, t_step);
                    measure_gradient_flow_stuff(V, p.gflow_l_step, t_step);
                    t_step0 = t_step;

                    for (int i = 1; i < nflow_steps; ++i) {

                        t_step =
                            do_gradient_flow_adapt(V, i * p.gflow_l_step, (i + 1) * p.gflow_l_step,
                                                   p.gflow_a_accu, p.gflow_r_accu, t_step);

                        measure_gradient_flow_stuff(V, (i + 1) * p.gflow_l_step, t_step);
                    }

                    hila::out0 << "Gflow_end " << gtrajectory << "    time " << std::setprecision(3)
                               << hila::gettime() - gftime << '\n';

                    gf_timer.stop();
                }
            }
        }

        if (p.n_save > 0 && (trajectory + 1) % p.n_save == 0) {
            checkpoint(U, trajectory, p);
        }
    }

    hila::finishrun();
}
