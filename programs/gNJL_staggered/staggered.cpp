/*********************************************
 * Simulates staggered fermions with a gauge *
 * interaction                               *
 *********************************************/


#include "staggered.h"





void measure(field<SU<N>> (&gauge)[NDIM], field<double> &sigma, field<double> &pi){
  double plaq = plaquette(gauge);
  output0 << "Plaq: " << plaq << "\n";

  double sigmasq = 0, sigma_ave = 0;
  double pisq = 0, pi_ave = 0;
  onsites(ALL){
    sigma_ave += sigma[X];
    sigmasq += sigma[X]*sigma[X];
    pi_ave += pi[X];
    pisq += pi[X]*pi[X];
  }
  sigma_ave /= lattice->volume();
  sigmasq /= lattice->volume();
  pi_ave /= lattice->volume();
  pisq /= lattice->volume();
  output0 << "Sigma: " << sigma_ave << "\n";
  output0 << "Sigma sq: " << sigmasq << "\n";
  output0 << "Pi: " << pi_ave << "\n";
  output0 << "Pi sq: " << pisq << "\n";
}




int main(int argc, char **argv){

  // Read parameters
  input parameters = input();
  parameters.import("parameters");
  double beta = parameters.get("beta");
  double gamma = parameters.get("gamma");
  double mass = parameters.get("mass");
  int seed = parameters.get("seed");
	double hmc_steps = parameters.get("hmc_steps");
	double traj_length = parameters.get("traj_length");

  lattice->setup( nd[0], nd[1], nd[2], nd[3], argc, argv );
  seed_random(seed);


  // Define gauge field and momentum field
  field<SUN> gauge[NDIM];
  field<SUN> momentum[NDIM];

  // Define gauge and momentum action terms
  gauge_momentum_action gma(gauge, momentum);
  gauge_action ga(gauge, momentum, beta);

  // Initialize the gauge field
  ga.set_unity();
  
  // Define auxiliary fields and momentum
  field<double> sigma, pi, sigma_mom, pi_mom;

  // And action
  auxiliary_momentum_action ama(sigma, pi, sigma_mom, pi_mom);
  auxiliary_action aa(sigma, pi, sigma_mom, pi_mom, gamma);

  // Initialize
  sigma[ALL] = 0.05;
  pi[ALL] = 0;

  // Define a Dirac operator (2 flavors)
  dirac_staggered_gNJL<VEC, SUN> D(mass, gauge, sigma, pi);
  gNJL_fermion_action fa(D, momentum, sigma_mom, pi_mom);

  // Build two integrator levels. Gauge is on the lowest level and
  // the fermions are on higher level
  integrator integrator_level_1(ga+aa, gma+ama);
  integrator integrator_level_2(fa, integrator_level_1);


  // Run HMC using the integrator
  for(int step = 0; step < 100; step ++){
    // Run update
    update_hmc(integrator_level_2, hmc_steps, traj_length);
    measure(gauge, sigma, pi);
  }



  finishrun();

  return 0;
}
