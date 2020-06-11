/*********************************************
 * Simulates Wilson fermions with a gauge *
 * interaction                               *
 *********************************************/

#include "Wilson.h"




int main(int argc, char **argv){

  input parameters = input();
  parameters.import("parameters");
  double beta = parameters.get("beta");
  double kappa = parameters.get("kappa");
  int seed = parameters.get("seed");
	int n_trajectories = parameters.get("n_trajectories");
	double hmc_steps = parameters.get("hmc_steps");
	double traj_length = parameters.get("traj_length");
	std::string configfile = parameters.get("configuration_file");

  lattice->setup( nd[0], nd[1], nd[2], nd[3], argc, argv );
  seed_random(seed);

  // Define gauge field and momentum field
  field<SUN> gauge[NDIM];
  field<SUN> momentum[NDIM];
  field<adjoint<N,double>> adj_gauge[NDIM];

  // Define gauge and momentum action terms
  gauge_momentum_action ma(gauge, momentum);
  gauge_action ga(gauge, momentum, beta);

  ga.set_unity();

  // Define a Dirac operator
  Dirac_Wilson_evenodd<N*N-1, double, adjoint<N,double>> D(kappa, adj_gauge);
  high_representation_fermion_action fa(D, momentum, gauge, adj_gauge);


  // Build two integrator levels. Gauge is on the lowest level and
  // the fermions are on higher level
  integrator integrator_level_1(ga, ma);
  integrator integrator_level_2(fa, integrator_level_1);
  
  // Initialize the gauge field
  ga.set_unity();

  int config_found = (bool) std::ifstream(configfile);
  broadcast(config_found);
  if( config_found )
  {
    output0 << "Found configuration file, reading\n";
    read_fields(configfile, gauge[0], gauge[1], gauge[2], gauge[3]);
  } else {
    output0 << "No config file " << configfile << ", starting new run\n";
  }

  // Run HMC using the integrator
  for(int step = 0; step < n_trajectories; step ++){
    update_hmc(integrator_level_2, hmc_steps, traj_length);
    double plaq = plaquette(ga.gauge);
    output0 << "Plaq: " << plaq << "\n";

    write_fields(configfile, gauge[0], gauge[1], gauge[2], gauge[3]);

  }

  finishrun();

  return 0;
}
