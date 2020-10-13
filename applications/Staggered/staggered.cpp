/*********************************************
 * Simulates staggered fermions with a gauge *
 * interaction                               *
 *********************************************/

//#define DEBUG_CG
#include "staggered.h"




int main(int argc, char **argv){

  input parameters = input();
  parameters.import("parameters");
  double beta = parameters.get("beta");
  double mass = parameters.get("mass");
  int seed = parameters.get("seed");
	int n_trajectories = parameters.get("n_trajectories");
	double hmc_steps = parameters.get("hmc_steps");
	double traj_length = parameters.get("traj_length");
	std::string configfile = parameters.get("configuration_file");

  lattice->setup( nd[0], nd[1], nd[2], nd[3], argc, argv );
  seed_random(seed);

  // Define gauge field and momentum field
  gauge_field<SU<N,double>> gauge;
  symmetric_gauge_field<N,double> sym_gauge(gauge);
  antisymmetric_gauge_field<N,double> antisym_gauge(gauge);
  adjoint_gauge_field<N,double> adj_gauge(gauge);

  // Initialize the gauge field
  gauge.set_unity();

  // Define gauge and momentum action terms
  gauge_action ga(gauge, beta);
  gauge_momentum_action ma(gauge, beta);

  // Define a Dirac operator
  dirac_staggered_evenodd D(mass, gauge);
  fermion_action fa(D, gauge);

  // A second fermion, for checking that addition works
  dirac_staggered_evenodd D2(mass, antisym_gauge);
  fermion_action fa2(D2, antisym_gauge);

  dirac_staggered_evenodd D3(mass, sym_gauge);
  fermion_action fa3(D3, sym_gauge);

  dirac_staggered_evenodd D4(mass, adj_gauge);
  fermion_action fa4(D4, adj_gauge); 

  // Build two integrator levels. Gauge is on the lowest level and
  // the fermions are on higher level
  O2_integrator integrator_level_1(ga, ma);
  action_base fsum = fa+fa2+fa3+fa4;
  O2_integrator integrator_level_2(fsum, integrator_level_1);


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
    update_hmc(integrator_level_2, hmc_steps, traj_length);
    double plaq = gauge.plaquette();
    output0 << "Plaq: " << plaq << "\n";

    gauge.write_file(configfile);

  }



  finishrun();

  return 0;
}
