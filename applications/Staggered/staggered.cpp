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
  field<SU<N,double>> gauge[NDIM];
  field<antisymmetric<N, double>> antisym_gauge[NDIM];
  field<symmetric<N, double>> sym_gauge[NDIM];
  field<adjoint<N, double>> adj_gauge[NDIM];
  field<SUN> momentum[NDIM];

  foralldir(dir){
    onsites(ALL){
      antisym_gauge[dir][X] = 1;
      sym_gauge[dir][X] = 1;
      adj_gauge[dir][X] = 1;
    }
  }

  // Define gauge and momentum action terms
  gauge_momentum_action ma(gauge, momentum);
  gauge_action ga(gauge, momentum, beta);
  
  // Define a Dirac operator
  dirac_staggered_evenodd<SU_vector<N,double>, SU<N,double>> D(mass, gauge);
  fermion_action fa(D, momentum);

  // A second fermion, for checking that addition works
  dirac_staggered_evenodd<SU_vector<N*(N-1)/2,double>, antisymmetric<N, double>> D2(mass, antisym_gauge);
  high_representation_fermion_action fa2(D2, momentum, gauge, antisym_gauge);

  dirac_staggered_evenodd<SU_vector<N*(N+1)/2,double>, symmetric<N, double>> D3(mass, sym_gauge);
  high_representation_fermion_action fa3(D3, momentum, gauge, sym_gauge);

  dirac_staggered_evenodd<SU_vector<N*N-1,double>, adjoint<N, double>> D4(mass, adj_gauge);
  high_representation_fermion_action fa4(D4, momentum, gauge, adj_gauge);

  // Build two integrator levels. Gauge is on the lowest level and
  // the fermions are on higher level
  integrator integrator_level_1(ga, ma);
  integrator integrator_level_2(fa+fa2+fa3+fa4, integrator_level_1);
  
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
