/*********************************************
 * Simulates Wilson fermions with a gauge *
 * interaction                               *
 *********************************************/

#define DEBUG_CG

#include "Wilson.h"




int main(int argc, char **argv){

  input parameters = input();
  parameters.import("parameters");
  double beta = parameters.get("beta");
  double kappa = parameters.get("kappa");
  int seed = parameters.get("seed");
	double hmc_steps = parameters.get("hmc_steps");
	double traj_length = parameters.get("traj_length");

  lattice->setup( nd[0], nd[1], nd[2], nd[3], argc, argv );
  seed_random(seed);

  // Define gauge field and momentum field
  field<SUN> gauge[NDIM];
  field<SUN> momentum[NDIM];

  // Define gauge and momentum action terms
  gauge_momentum_action ma(gauge, momentum);
  gauge_action ga(gauge, momentum, beta);

  ga.set_unity();

  // Define a Dirac operator
  Dirac_Wilson_evenodd<VEC, SUN> D(kappa, gauge);
  fermion_action fa(D, gauge, momentum);


  // Build two integrator levels. Gauge is on the lowest level and
  // the fermions are on higher level
  integrator integrator_level_1(ga, ma);
  integrator integrator_level_2(fa, integrator_level_1);
  
  // Initialize the gauge field
  ga.set_unity();

  // Run HMC using the integrator
  for(int step = 0; step < 5; step ++){
    update_hmc(integrator_level_2, hmc_steps, traj_length);
    double plaq = plaquette(ga.gauge);
    output0 << "Plaq: " << plaq << "\n";
  }

  finishrun();

  return 0;
}
