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
  gauge_field<SU<N,double>> gauge;
  adjoint_gauge_field<N,double> adj_gauge(gauge);

  // Initialize the gauge field
  gauge.set_unity();

  // Define gauge and momentum action terms
  gauge_action ga(gauge, beta);

  // Define a Dirac operator
  Dirac_Wilson_evenodd D(kappa, adj_gauge);
  Hasenbusch_action_1 fa1(D, adj_gauge, 1.0);
  Hasenbusch_action_2 fa2(D, adj_gauge, 1.0);

  // Build two integrator levels. Gauge is on the lowest level and
  // the fermions are on higher level
  integrator integrator_level_2(fa1, ga);
  integrator integrator_level_3(fa2, integrator_level_2);
  
  int config_found = (bool) std::ifstream(configfile);
  broadcast(config_found);
  if( config_found )
  {
    output0 << "Found configuration file, reading\n";
    gauge.read_file(configfile);
  } else {
    output0 << "No config file " << configfile << ", starting new run\n";
  }

  // Run HMC using the integrator
  for(int step = 0; step < n_trajectories; step ++){
    update_hmc(integrator_level_3, hmc_steps, traj_length);
    double plaq = gauge.plaquette();
    output0 << "Plaq: " << plaq << "\n";

    gauge.write_file(configfile);

  }

  finishrun();

  return 0;
}
