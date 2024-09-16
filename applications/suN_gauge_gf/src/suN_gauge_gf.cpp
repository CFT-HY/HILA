#include "hila.h"
#include "gauge/polyakov.h"
#include "gauge/wilson_plaquette_action.h"
#include "gauge/sun_heatbath.h"
#include "gauge/sun_overrelax.h"
#include "gauge/staples.h"
#include "gauge/gradient_flow.h"
#include "tools/string_format.h"
#include "tools/floating_point_epsilon.h"

#ifndef NCOLOR
#define NCOLOR 3
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
 */
template <typename group>
void update(GaugeField<group> &U, const parameters &p, bool relax) {
    for (int i = 0; i < 2 * NDIM; ++i) {
        int tdp = hila::broadcast((int)(hila::random() * 2 * NDIM));
        int tdir = tdp / 2;
        int tpar = 1 + (tdp % 2);
        //hila::out0 << "   " << Parity(tpar) << " -- " << Direction(tdir);
        update_parity_dir(U, p, Parity(tpar), Direction(tdir), relax);
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
 */
template <typename group>
void update_parity_dir(GaugeField<group> &U, const parameters &p, Parity par, Direction d,
                       bool relax) {

    static hila::timer hb_timer("Heatbath");
    static hila::timer or_timer("Overrelax");
    static hila::timer staples_timer("Staplesum");

    Field<group> staples;

    staples_timer.start();

    staplesum(U, staples, d, par);
    staples_timer.stop();

    if (relax) {

        or_timer.start();

        onsites(par) {
#ifdef SUN_OVERRELAX_dFJ
            suN_overrelax_dFJ(U[d][X], staples[X], p.beta);
#else
            suN_overrelax(U[d][X], staples[X]);
#endif
        }
        or_timer.stop();

    } else {

        hb_timer.start();
        onsites(par) {
            suN_heatbath(U[d][X], staples[X], p.beta);
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
void do_trajectory(GaugeField<group> &U, const parameters &p) {

    for (int n = 0; n < p.n_update; n++) {
        for (int i = 0; i <= p.n_overrelax; i++) {
            bool relax = hila::broadcast((int)(hila::random() * (1 + p.n_overrelax)) != 0);
            update(U, p, relax);
            //hila::out0 << relax << "\n";
        }
    }
    U.reunitarize_gauge();
}

// heat-bath functions
///////////////////////////////////////////////////////////////////////////////////
// measurement functions

template <typename group>
void measure_stuff(const GaugeField<group> &U) {
    // perform measurements on current gauge and momentum pair (U, E) and
    // print results in formatted form to standard output
    static bool first = true;
    if (first) {
        // print legend for measurement output
        hila::out0 << "LMEAS:        plaq        P.real        P.imag\n";
        first = false;
    }
    auto plaq = measure_s_wplaq(U) / (lattice.volume() * NDIM * (NDIM - 1) / 2);
    auto poly = measure_polyakov(U, e_t);
    hila::out0 << string_format("MEAS % 0.6e % 0.6e % 0.6e", plaq, poly.real(), poly.imag())
               << '\n';
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

        do_trajectory(U, p);

        // put sync here in order to get approx gpu timing
        hila::synchronize_threads();
        update_timer.stop();

        measure_timer.start();

        hila::out0 << "Measure_start " << trajectory << '\n';

        measure_stuff(U);

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
