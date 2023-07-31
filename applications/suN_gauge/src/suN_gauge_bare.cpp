#include "suN_gauge_bare.h"
#include "parameters.h"
#include "checkpoint.h"

/**
 * @brief Helper function to get valid z-coordinate index
 *
 * @param z
 * @return int
 */
int z_ind(int z) {
    return (z + lattice.size(e_z)) % lattice.size(e_z);
}

template <typename group>
void measure_stuff(const GaugeField<group> &U, const parameters &p) {

    static bool first = true;

    if (first) {
        hila::out0 << "Legend:";
        hila::out0 << " plaq  P.real  P.imag\n";

        first = false;
    }

    auto poly = measure_polyakov(U);

    auto plaq = U.measure_plaq() / (lattice.volume() * NDIM * (NDIM - 1) / 2);

    hila::out0 << "MEAS " << std::setprecision(8);

    // write the -(polyakov potential) first, this is used as a weight factor in aa

    hila::out0 << plaq << ' ' << poly << '\n';
}

template <typename group>
void update(GaugeField<group> &U, const parameters &p, bool relax) {

    foralldir(d) {
        for (Parity par : {EVEN, ODD}) {

            update_parity_dir(U, p, par, d, relax);
        }
    }
}

template <typename group>
void update_parity_dir(GaugeField<group> &U, const parameters &p, Parity par, Direction d,
                       bool relax) {

    static hila::timer hb_timer("Heatbath");
    static hila::timer or_timer("Overrelax");
    static hila::timer staples_timer("Staplesum");

    Field<group> staples;

    staples_timer.start();
    if (p.deltab != 0)
        staplesum_db(U, staples, d, par, p.deltab);
    else
        staplesum(U, staples, d, par);
    staples_timer.stop();

    if (relax) {

        or_timer.start();
        onsites(par) {
            suN_overrelax(U[d][X], staples[X]);
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


template <typename group>
void do_trajectory(GaugeField<group> &U, const parameters &p) {

    for (int n = 0; n < p.n_update; n++) {
        for (int i = 0; i < p.n_overrelax; i++) {
            update(U, p, true);
        }
        update(U, p, false);
    }
    U.reunitarize_gauge();
}


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

    parameters p;

    hila::out0 << "SU(" << mygroup::size() << ") heat bath + overrelax update\n";

    hila::input par("parameters");

    CoordinateVector lsize;
    lsize = par.get("lattice size"); // reads NDIM numbers

    p.beta = par.get("beta");
    // deltab sets system to different beta on different sides, by beta*(1 +- deltab)
    // use for initial config generation only
    p.deltab = par.get("delta beta fraction");
    // trajectory length in steps
    p.n_overrelax = par.get("overrelax steps");
    p.n_update = par.get("updates in trajectory");
    p.n_trajectories = par.get("trajectories");
    p.n_thermal = par.get("thermalization");

    // random seed = 0 -> get seed from time
    long seed = par.get("random seed");
    // save config and checkpoint
    p.n_save = par.get("traj/saved");
    // measure surface properties and print "profile"
    p.config_file = par.get("config name");

    // if polyakov range is off, do nothing with
    int p_item = par.get_item("polyakov potential", {"off", "min", "range"});
    if (p_item == 0) {
        p.polyakov_pot = poly_limit::OFF;
    } else if (p_item == 1) {
        p.polyakov_pot = poly_limit::PARABOLIC;
        p.poly_min = par.get();
        p.poly_m2 = par.get("mass2");
    } else {
        p.polyakov_pot = poly_limit::RANGE;
        Vector<2, double> r;
        r = par.get();
        p.poly_min = r[0];
        p.poly_max = r[1];
    }

    if (par.get_item("updates/profile meas", {"off", "%i"}) == 1) {
        p.n_profile = par.get();
    } else {
        p.n_profile = 0;
    }

    if (p.n_profile) {
        p.n_smear = par.get("smearing steps");
        p.smear_coeff = par.get("smear coefficient");
        p.z_smear = par.get("z smearing steps");
        p.n_surface = par.get("traj/surface");
        p.n_dump_polyakov = par.get("traj/polyakov dump");

        if (p.n_smear.size() != p.z_smear.size()) {
            hila::out0 << "Error in input file: number of values in 'smearing steps' != 'z "
                          "smearing steps'\n";
            hila::terminate(0);
        }

    } else {
        p.n_dump_polyakov = 0;
    }

    par.close(); // file is closed also when par goes out of scope

    // setting up the lattice is convenient to do after reading
    // the parameter
    lattice.setup(lsize);

    // Alloc gauge field
    GaugeField<mygroup> U;

    U = 1;

    // use negative trajectory for thermal
    int start_traj = -p.n_thermal;

    hila::timer update_timer("Updates");
    hila::timer measure_timer("Measurements");

    // restore_checkpoint(U, start_traj, p);

    // We need random number here
    if (!hila::is_rng_seeded())
        hila::seed_random(seed);

    for (int trajectory = start_traj; trajectory < p.n_trajectories; trajectory++) {

        double ttime = hila::gettime();

        update_timer.start();

        double acc = 0;

        do_trajectory(U, p);

        // put sync here in order to get approx gpu timing
        hila::synchronize_threads();
        update_timer.stop();

        // trajectory is negative during thermalization
        if (trajectory >= 0) {
            measure_timer.start();

            hila::out0 << "Measure_start " << trajectory << '\n';

            measure_stuff(U, p);

            hila::out0 << "Measure_end " << trajectory << std::endl;

            measure_timer.stop();
        }

        // if (p.n_save > 0 && (trajectory + 1) % p.n_save == 0) {
        //     checkpoint(U, trajectory, p);
        // }
    }

    hila::finishrun();
}
