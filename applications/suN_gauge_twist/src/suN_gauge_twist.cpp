#include "gauge/staples.h"
#include "gauge/polyakov.h"
#include "gauge/stout_smear.h"
#include "gauge/sun_heatbath.h"
#include "gauge/sun_overrelax.h"
#include "hila.h"
#include "multicanonical.h"

//#include "gauge/polyakov.h"

#include <fftw3.h>
#include <numeric>
#include <stdio.h>
#include <stdlib.h>

// local includes
#include "checkpoint.h"
#include "parameters.h"
#include "polyakov_surface.h"
#include "twist_specific_methods.h"
#include "utility.h"
// /**
//  * @brief Helper function to get valid z-coordinate index
//  *
//  * @param z
//  * @return int
//  */
// int z_ind(int z) { return (z + lattice.size(e_z)) % lattice.size(e_z); }

/**
 * @brief Measures Polyakov lines and Wilson action
 *
 * @tparam group
 * @param U GaugeField to measure
 * @param p Parameter struct
 */
template <typename group>
void measure_plaq(const GaugeField<group> &U, const parameters &p) {


  auto plaq_z = measure_plaq_with_z(
      U, p.twist_coeff);

  print_formatted_numbers(plaq_z, "plaquette z", false, true);
  hila::out0 << "plaquette: " << plaq_z.back() << '\n';
}

/**
 * @brief Measures Polyakov lines and Wilson action with multicanonical weights
 *
 * @tparam group
 * @param U GaugeField to measure
 * @param p Parameter struct
 */
template <typename group>
void measure_plaq_multicanonical(const GaugeField<group> &U,
                                 const parameters &p) {

  auto plaq = measure_plaq_with_z(
      U, p.twist_coeff); /// (lattice.volume() * NDIM * (NDIM - 1) / 2);

  print_formatted_numbers(plaq, "plaquette", false, true);
  hila::out0 << "muca_plaquette: " << hila::muca::weight(plaq.back()) << '\n';
}

template <typename group>
void measure_poly(const GaugeField<group> &U, const parameters &p) {
  auto poly_z = measure_polyakov_with_z(U);
  auto poly_abs = measure_polyakov_with_z_abs(U);

  print_formatted_numbers(poly_z, "polyakov z", false, true);
  print_formatted_numbers(poly_abs, "polyakov abs z", false, true);
  hila::out0 << "polyakov: " << poly_z.back() << '\n';
  hila::out0 << "polyakov abs: " << poly_abs.back() << '\n';
}

template <typename group>
void measure_poly_multicanonical(const GaugeField<group> &U,
                                 const parameters &p) {
  auto poly = measure_polyakov_with_z(U);
  auto poly_abs = measure_polyakov_with_z_abs(U);

  print_formatted_numbers(poly, "polyakov z", false, true);
  print_formatted_numbers(poly_abs, "polyakov abs z", false, true);
  hila::out0 << "polyakov: " << poly.back() << '\n';
  hila::out0 << "polyakov abs: " << poly_abs.back() << '\n';

  hila::out0 << "muca_polyakov: " << hila::muca::weight(poly_abs.back())
             << '\n';
}

/**
 * @brief Wrapper update function
 * @details Updates Gauge Field one direction at a time first EVEN then ODD
 * parity
 *
 * @tparam group
 * @param U GaugeField to update
 * @param p Parameter struct
 * @param relax If true evolves GaugeField with over relaxation if false then
 * with heat bath
 */
template <typename group>
void update(GaugeField<group> &U, const parameters &p, bool relax) {
  foralldir(d) {
    for (Parity par : {EVEN, ODD}) {

      update_parity_dir(U, p, par, d, relax);
    }
  }
}

template <typename group>
void update_multicanonical(GaugeField<group> &U, const parameters &p,
                           bool relax) {
  foralldir(d) {
    for (Parity par : {EVEN, ODD}) {
      if (p.muca_poly) {
        if (d == e_t) {
          auto U_old = U;
          update_parity_dir(U, p, par, d, relax);
          auto OP = measure_polyakov_with_z_abs(U);
          auto OP_old = measure_polyakov_with_z_abs(U_old);
          // hila::out0 << "OP old: " << OP_old.back() << " OP new: " <<
          // OP.back() << std::endl;
          if (!hila::muca::accept_reject(OP_old.back(), OP.back())) {
            U = U_old;
          }
        } else {
          update_parity_dir(U, p, par, d, relax);
        }
      }
      if (p.muca_action) {
        auto U_old = U;
        update_parity_dir(U, p, par, d, relax);
        auto OP = measure_plaq_with_z(U, p.twist_coeff);
        auto OP_old = measure_plaq_with_z(U_old, p.twist_coeff);
        // hila::out0 << "OP old: " << OP_old.back() << " OP new: " << OP.back()
        // << std::endl;
        if (!hila::muca::accept_reject(OP_old.back(), OP.back())) {
          U = U_old;
        }
      }
    }
  }
}

/**
 * @brief Wrapper function to updated GaugeField per direction
 * @details Computes first staplesum, then uses computed result to evolve
 * GaugeField either with over relaxation or heat bath
 *
 * @tparam group
 * @param U GaugeField to evolve
 * @param p parameter struct
 * @param par Parity
 * @param d Direction to evolve
 * @param relax If true evolves GaugeField with over relaxation if false then
 * with heat bath
 */
template <typename group>
void update_parity_dir(GaugeField<group> &U, const parameters &p, Parity par,
                       Direction d, bool relax) {

  static hila::timer hb_timer("Heatbath");
  static hila::timer or_timer("Overrelax");
  static hila::timer staples_timer("Staplesum");

  Field<group> staples;

  staples_timer.start();

  // staplesum_twist(U, staples, d, p.twist_coeff,par);
  staplesum_twist(U, staples, d, p.twist_coeff, par);

  staples_timer.stop();

  if (relax) {

    or_timer.start();
    onsites(par) { suN_overrelax(U[d][X], staples[X]); }
    or_timer.stop();

  } else {

    hb_timer.start();
    onsites(par) { suN_heatbath(U[d][X], staples[X], p.beta); }
    hb_timer.stop();
  }
}

/**
 * @brief Evolve gauge field
 * @details Evolution happens by means of heat bath and over relaxation. For
 * each heatbath update (p.n_update) we update p.n_overrelax times with over
 * relaxation.
 *
 * @tparam group
 * @param U
 * @param p
 */
template <typename group>
void do_trajectory(GaugeField<group> &U, const parameters &p) {
  auto U_old = U;
  for (int n = 0; n < p.n_update; n++) {
    for (int i = 0; i < p.n_overrelax; i++) {
      update(U, p, true);
    }

    update(U, p, false);
  }
  U.reunitarize_gauge();
}

/**
 * @brief Evolve gauge field with multicanonical update
 * @details Evolution happens by means of heat bath and over relaxation. For
 * each heatbath update (p.n_update) we update p.n_overrelax times with over
 * relaxation.
 *
 * @tparam group
 * @param U
 * @param p
 */
template <typename group>
void do_trajectory_multicanonical(GaugeField<group> &U, const parameters &p) {
  for (int n = 0; n < p.n_update; n++) {
    for (int i = 0; i < p.n_overrelax; i++) {
      update_multicanonical(U, p, true);
    }

    update_multicanonical(U, p, false);
  }
  U.reunitarize_gauge();
}

/**
 * @brief Create weight function with multicanonical method
 *
 * @param U Gauge field
 * @param p Parameters struct
 */
void iterate_weights_multicanonical(GaugeField<mygroup> U,
                                    const parameters &p) {

  hila::out0 << "here also\n";
  bool iterate_status = true;
  while (iterate_status) {
    for (int i = 0; i < p.muca_updates; i++) {
      do_trajectory_multicanonical(U, p);
    }
    if (p.muca_action) {
      auto OP = measure_plaq_with_z(U, p.twist_coeff);
      auto P = measure_polyakov_with_z(U);
      hila::out0 << "Order parameter: " << OP.back() << std::endl;
      hila::out0 << "polyakov muca: " << P.back() << std::endl;
      iterate_status = hila::muca::iterate_weights(OP.back());
    } else {
      auto OP = measure_polyakov_with_z_abs(U);
      auto Plaq = measure_plaq_with_z(U, p.twist_coeff);
      auto Poly = measure_polyakov_with_z(U);
      hila::out0 << "Order parameter: " << abs(OP.back()) << std::endl;
      hila::out0 << "plaquette muca: "<< Plaq.back() << std::endl;
      hila::out0 << "polyakov muca: "<< Poly.back() << std::endl;

      iterate_status = hila::muca::iterate_weights(OP.back());
    }
  }
  hila::muca::write_weight_function(hila::muca::generate_outfile_name());
}

int main(int argc, char **argv) {
  parameters p;
  CoordinateVector lsize;
  hila::cmdline.add_flag("-beta", "Temparature");
  hila::cmdline.add_flag("-twist", "Twist coefficient");
  hila::cmdline.add_flag("-config", "Config filename");
  hila::cmdline.add_flag("-ntraj", "Number of trajectories");
  hila::cmdline.add_flag("-lsize", "Lattice size");
  hila::initialize(argc, argv);

  hila::out0 << "SU(" << mygroup::size() << ") heat bath + overrelax update\n";

  hila::input par("parameters");

  lsize = par.get("lattice size"); // reads NDIM numbers

  p.beta = par.get("beta");

  // deltab sets system to different beta on different sides, by beta*(1 +-
  // deltab) use for initial config generation only
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
  p.twist_coeff = par.get("twist_coeff");

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
      hila::out0 << "Error in input file: number of values in 'smearing "
                    "steps' != 'z "
                    "smearing steps'\n";
      hila::terminate(0);
    }

  } else {
    p.n_dump_polyakov = 0;
  }
  p.muca_action = par.get("muca_action");
  p.muca_poly = par.get("muca_poly");
  p.muca_updates = par.get("muca_updates");
  p.measure_surface = par.get("measure_surface");
  std::string initial_state = par.get("initial_state");

  par.close(); // file is closed also when par goes out of scope

  if (hila::cmdline.flag_present("-beta")) {
    p.beta = hila::cmdline.get_double("-beta");
    hila::out0 << "beta            " << p.beta << std::endl;
  }
  if (hila::cmdline.flag_present("-twist")) {
    p.twist_coeff = hila::cmdline.get_int("-twist");
    hila::out0 << "twist_coeff     " << p.twist_coeff << std::endl;
  }
  if (hila::cmdline.flag_present("-config")) {
    p.config_file = hila::cmdline.get_string("-config");
  }
  if (hila::cmdline.flag_present("-ntraj")) {
    p.n_trajectories = hila::cmdline.get_int("-ntraj");
    hila::out0 << "Trajectories    " << p.twist_coeff << std::endl;
  }
  if (hila::cmdline.flag_present("-lsize")) {
    for (int i = 0; i < NDIM; i++) {
      lsize[i] = hila::cmdline.get_int("-lsize", i);
    }
    hila::out0 << "Lattice size    ";
    for (int i = 0; i < NDIM; i++) {
      hila::out0 << lsize[i] << " ";
    }
    hila::out0 << std::endl;
  }
  // setting up the lattice is convenient to do after reading
  // the parameter
  lattice.setup(lsize);
  //  Alloc gauge field

  // U = initial_state;
  //  use negative trajectory for thermal
  int start_traj = -p.n_thermal;

  hila::timer update_timer("Updates");
  hila::timer measure_timer("Measurements");
  hila::timer muca_timer("Muca");
  GaugeField<mygroup> U;

  // We need random number here
  if (!hila::is_rng_seeded())
    hila::seed_random(seed);

  if (initial_state == "cold") {
    foralldir(d) onsites(ALL) { U[d][X].gaussian_random(); }
  } else if (initial_state == "hot") {
    U = 1;
  } else if (initial_state == "both") {
    foralldir(d) onsites(EVEN) { U[d][X].gaussian_random(); }
    foralldir(d) onsites(ODD) { U[d][X] = 1; }
  }

  // restore_checkpoint(U, start_traj, p);

  //  muca_timer.start();

  if (p.muca_action || p.muca_poly) {
    hila::muca::initialise("muca_parameters");
    std::ofstream MFile;
    MFile.open("muca_measurements", std::ios_base::app);
    iterate_weights_multicanonical(U, p);
  }
  // muca_timer.stop();
  hila::out0 << "MEASURE start\n";
  void (*do_trajectory_ptr)(GaugeField<mygroup> &, const parameters &);
  void (*measure_plaquette_ptr)(const GaugeField<mygroup> &,
                                const parameters &);
  void (*measure_polyakov_ptr)(const GaugeField<mygroup> &, const parameters &);
  if (p.muca_action) {
    do_trajectory_ptr = do_trajectory_multicanonical<mygroup>;
    measure_plaquette_ptr = measure_plaq_multicanonical<mygroup>;
    measure_polyakov_ptr = measure_poly<mygroup>;
  } else if (p.muca_poly) {
    do_trajectory_ptr = do_trajectory_multicanonical<mygroup>;
    measure_plaquette_ptr = measure_plaq<mygroup>;
    measure_polyakov_ptr = measure_poly_multicanonical<mygroup>;
  } else {
    do_trajectory_ptr = do_trajectory<mygroup>;
    measure_plaquette_ptr = measure_plaq<mygroup>;
    measure_polyakov_ptr = measure_poly<mygroup>;
  }

  for (int trajectory = start_traj; trajectory < p.n_trajectories;
       trajectory++) {
    double ttime = hila::gettime();

    update_timer.start();

    double acc = 0;
    do_trajectory_ptr(U, p);

    // put sync here in order to get approx gpu timing
    hila::synchronize_threads();
    update_timer.stop();

    // trajectory is negative during thermalization
    if (trajectory >= 0) {
      measure_timer.start();

      measure_plaquette_ptr(U, p);

      measure_polyakov_ptr(U, p);

      if (p.measure_surface)
        measure_polyakov_surface(U, p, trajectory);

      // hila::out0 << "Measure_end " << trajectory << std::endl;

      measure_timer.stop();
    }

    if (p.n_save > 0 && (trajectory + 1) % p.n_save == 0) {
      checkpoint(U, trajectory, p);
    }
  }
  
  hila::out0 << "MEASURE end\n";
  hila::out0 << expi(4.0 / 3.0 * M_PI) << std::endl;
  hila::finishrun();
}
