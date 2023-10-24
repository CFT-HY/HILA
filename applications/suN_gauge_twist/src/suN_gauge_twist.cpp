#include "gauge/staples.h"
#include "gauge/stout_smear.h"
#include "gauge/sun_heatbath.h"
#include "gauge/sun_overrelax.h"
#include "hila.h"
//#include "gauge/polyakov.h"

#include <fftw3.h>
#include <numeric>

// local includes
#include "checkpoint.h"
#include "parameters.h"
#include "twist_specific_methods.h"

/**
 * @brief Helper function to get valid z-coordinate index
 *
 * @param z
 * @return int
 */
int z_ind(int z) { return (z + lattice.size(e_z)) % lattice.size(e_z); }

/**
 * @brief Measures Polyakov lines and Wilson action
 *
 * @tparam group
 * @param U GaugeField to measure
 * @param p Parameter struct
 */
template <typename group>
void measure_plaq(const GaugeField<group> &U, const parameters &p) {

  static bool first = true;

  auto plaq = measure_plaq_with_z(
      U, p.twist_coeff); /// (lattice.volume() * NDIM * (NDIM - 1) / 2);

  if (first) {
    hila::out0 << "plaquette: ";
    for (int i = 0; i < plaq.size() - 1; i++)
      hila::out0 << i << " ";
    hila::out0 << "sum \n";
    first = false;
  }

  // write the -(polyakov potential) first, this is used as a weight factor in
  // aa auto sum = std::reduce(plaq.begin(), plaq.end() - 1);
  hila::out0 << "plaquette: ";
  for (int i = 0; i < plaq.size(); i++) {
    hila::out0 << plaq[i] << " ";
  }
  hila::out0 << "\n";
}

template <typename group>
void measure_poly(const GaugeField<group> &U, const parameters &p) {
  static bool first = true;
  auto poly = measure_polyakov_twist(U);
  if (first) {
    hila::out0 << "polyakov: ";
    for (int i = 0; i < poly.size() - 1; i++)
      hila::out0 << i << ", ";
    hila::out0 << "sum \n";
    first = false;
  }
  hila::out0 << "polyakov: ";
  for (int i = 0; i < poly.size(); i++) {
    hila::out0 << poly[i] << ", ";
  }
  hila::out0 << "\n";
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

  restore_checkpoint(U, start_traj, p);

  // We need random number here
  if (!hila::is_rng_seeded())
    hila::seed_random(seed);
  hila::out0 << "MEASURE start\n";

  for (int trajectory = start_traj; trajectory < p.n_trajectories;
       trajectory++) {

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

      measure_plaq(U, p);

      measure_poly(U, p);

      // hila::out0 << "Measure_end " << trajectory << std::endl;

      measure_timer.stop();
    }

    if (p.n_save > 0 && (trajectory + 1) % p.n_save == 0) {
      checkpoint(U, trajectory, p);
    }
  }
  hila::out0 << "MEASURE end\n";

  hila::finishrun();
}
