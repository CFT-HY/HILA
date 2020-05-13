/*********************************************
 * Simulates staggered fermions with a gauge *
 * interaction                               *
 *********************************************/


#include "staggered.h"





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
  sigma[ALL] = 0;
  pi[ALL] = 0;

  // Define a Dirac operator (2 flavors)
  dirac_staggered_gNJL<N> D(mass, gauge, sigma, pi);
  gNJL_fermion_action fa(D, momentum, sigma_mom, pi_mom);


  // Build two integrator levels. Gauge is on the lowest level and
  // the fermions are on higher level
  integrator integrator_level_1(ga+aa, gma+ama);
  integrator integrator_level_2(fa, integrator_level_1);

  sigma_mom[ALL] = 0;
  pi_mom[ALL] = 0;

  double eps = 1e-6;
  double s = sigma.get_value_at(50);
  static field<SU_vector<N>> psi, chi, tmp, tmp2;
  onsites(ALL){
    psi[X].gaussian(); chi[X].gaussian();
  }

  double s1 = 0;
  D.apply(psi,tmp);
  onsites(ALL){
    s1 += chi[X].rdot(tmp[X]);
  }
  
  if(mynode()==0){
    sigma.set_value_at(s+eps, 50);
  }
  sigma.mark_changed(ALL);

  double s2 = 0;
  D.apply(psi,tmp);
  onsites(ALL){
    s2 += chi[X].rdot(tmp[X]);
  }

  if(mynode()==0){
    sigma.set_value_at(s, 50);
  }
  sigma.mark_changed(ALL);

  D.force(chi, psi, momentum, sigma_mom, pi_mom);
  double f = sigma_mom.get_value_at(50);
  double diff = f - (s2-s1)/eps;

  if(mynode()==0) {
    hila::output << "Action 1 " << s1 << "\n";
    hila::output << "Action 2 " << s2 << "\n";
    hila::output << "Calculated deriv " << f << "\n";
    hila::output << "Actual deriv " << (s2-s1)/eps << "\n";
    hila::output << "deriv diff " << diff << "\n";
    assert( diff*diff < eps*eps*1000 && "Fermion deriv" );
  }



  foralldir(dir){
    momentum[dir][ALL] = 0;
    gauge[dir][ALL] = 1;
  }
  sigma_mom[ALL] = 0;
  pi_mom[ALL] = 0;
  sigma[ALL] = 0;
  pi[ALL] = 0;

  s1 = fa.action();
  output0 << "Action 1 " << s1 << "\n";

  if(mynode()==0){
    sigma.set_value_at(s+eps,50);
  }
  sigma.mark_changed(ALL);
  s2 = fa.action();
  output0 << "Action 2 " << s2 << "\n";

  if(mynode()==0)
    sigma.set_value_at(s, 50);
  sigma.mark_changed(ALL);

  fa.force_step(1.0);
  f = sigma_mom.get_value_at(50);
  diff = diff = f - (s2-s1)/eps;

  if(mynode()==0) {
    hila::output << "Action 1 " << s1 << "\n";
    hila::output << "Action 2 " << s2 << "\n";
    hila::output << "Calculated deriv " << f << "\n";
    hila::output << "Actual deriv " << (s2-s1)/eps << "\n";
    hila::output << "deriv diff " << diff << "\n";
    assert( diff*diff < eps*eps*1000 && "Fermion deriv" );
  }



  // Run HMC using the integrator
  for(int step = 0; step < 100; step ++){
    // Run update
    update_hmc(integrator_level_2, hmc_steps, traj_length);
    
    // Measurements

    double plaq = plaquette(ga.gauge);
    output0 << "Plaq: " << plaq << "\n";

    double sigmasq = 0, sigma_ave = 0;
    onsites(ALL){
      sigma_ave += sigma[X];
      sigmasq += sigma[X]*sigma[X];
    }
    output0 << "Sigma: " << sigma_ave << "\n";
    output0 << "Sigma sq: " << sigmasq << "\n";
  }



  finishrun();

  return 0;
}
