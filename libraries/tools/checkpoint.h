/**
 * @file checkpoint.h
 * @brief Checkpointing functions for saving and restoring lattice configurations
 *
 */
#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include "hila.h"

/// Functions checkpoint / restore_checkpoint allow one to save lattice config periodically
/// Checkpoint keeps file "run_status" which holds the current trajectory.
/// By modifying "run status" the number of trajectories can be changed

/// parameters-type must


template <typename group>
void checkpoint(const GaugeField<group> &U, const std::string &config_file, int &n_trajectories,
                int trajectory, bool save_old = true) {

    double t = hila::gettime();

    if (save_old && hila::myrank() == 0 && filesys_ns::exists(config_file)) {
        filesys_ns::rename(config_file, config_file + ".prev");
        // rename config to config.prev
    }

    // save config
    U.config_write(config_file);

    if (hila::myrank() == 0) {

        // check if n_trajectories has changed
        hila::input status;
        status.quiet();
        if (status.open("run_status", false, false)) {
            int ntraj = status.get("trajectories");
            if (ntraj != n_trajectories) {
                hila::out0 << "* NUMBER OF TRAJECTORIES " << n_trajectories << " -> " << ntraj
                           << '\n';
            }
            n_trajectories = ntraj;
            status.close();
        }

        // write the status file
        std::ofstream outf;
        outf.open("run_status", std::ios::out | std::ios::trunc);
        outf << "trajectories " << n_trajectories
             << "   # CHANGE TO ADJUST NUMBER OF TRAJECTORIES IN THIS RUN\n";
        outf << "trajectory   " << trajectory + 1 << '\n';
        outf << "seed         " << static_cast<uint64_t>(hila::random() * (1UL << 61)) << '\n';
        outf << "time         " << hila::gettime() << '\n';
        outf.close();

        std::stringstream msg;
        msg << "Checkpointing, time " << hila::gettime() - t;
        hila::timestamp(msg.str());
    }

    // sync trajectory numbers
    hila::broadcast(n_trajectories);
}


template <typename group>
bool restore_checkpoint(GaugeField<group> &U, const std::string &config_file, int &n_trajectories,
                        int &trajectory) {
    uint64_t seed;
    bool ok = true;
    hila::input status;
    if (status.open("run_status", false, false)) {
        hila::out0 << "RESTORING FROM CHECKPOINT:\n";
        int traj = status.get("trajectories");
        if (traj != n_trajectories) {
            hila::out0 << "Number of trajectories set in 'run_status' " << n_trajectories << " -> "
                       << traj << '\n';
            n_trajectories = traj;
        }

        trajectory = status.get("trajectory");
        seed = status.get("seed");
        // p.time_offset = status.get("time");

        status.close();
        hila::seed_random(seed);

        U.config_read(config_file);
        ok = true;

    } else {

        bool exists = hila::myrank() == 0 && filesys_ns::exists(config_file);
        hila::broadcast(exists);
        if (exists) {
            hila::out0 << "READING initial config\n";
            U.config_read(config_file);
            ok = true;
        } else {
            ok = false;
        }
    }
    return ok;
}

#endif