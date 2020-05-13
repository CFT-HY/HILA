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
  sigma[ALL] = gaussian_ran();
  pi[ALL] = gaussian_ran();

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
    hila::output << "Calculated deriv " << f << "\n";
    hila::output << "Actual deriv " << (s2-s1)/eps << "\n";
    hila::output << "deriv diff " << diff << "\n";
    assert( diff*diff < 10*eps && "Fermion deriv" );
  }

  
  
  double p = pi.get_value_at(50);
  s1 = 0;
  D.apply(psi,tmp);
  onsites(ALL){
    s1 += chi[X].rdot(tmp[X]);
  }
  
  if(mynode()==0){
    pi.set_value_at(p+eps, 50);
  }
  pi.mark_changed(ALL);

  s2 = 0;
  D.apply(psi,tmp);
  onsites(ALL){
    s2 += chi[X].rdot(tmp[X]);
  }

  if(mynode()==0){
    pi.set_value_at(p, 50);
  }
  pi.mark_changed(ALL);

  f = pi_mom.get_value_at(50);
  diff = f - (s2-s1)/eps;

  if(mynode()==0) {
    hila::output << "Calculated deriv " << f << "\n";
    hila::output << "Actual deriv " << (s2-s1)/eps << "\n";
    hila::output << "deriv diff " << diff << "\n";
    assert( diff*diff < 10*eps && "Fermion deriv" );
  }



  p = pi.get_value_at(50);
  s1 = 0;
  D.dagger(psi,tmp);
  onsites(ALL){
    s1 += chi[X].rdot(tmp[X]);
  }
  
  if(mynode()==0){
    pi.set_value_at(p+eps, 50);
  }
  pi.mark_changed(ALL);
 
  s2 = 0;
  D.dagger(psi,tmp);
  onsites(ALL){
    s2 += chi[X].rdot(tmp[X]);
  }

  if(mynode()==0){
    pi.set_value_at(p, 50);
  }
  pi.mark_changed(ALL);

  sigma_mom[ALL] = 0;
  pi_mom[ALL] = 0;
  D.force(chi, psi, momentum, sigma_mom, pi_mom);
  f = -pi_mom.get_value_at(50);
  diff = f - (s2-s1)/eps;

  if(mynode()==0) {
    hila::output << "Calculated deriv " << f << "\n";
    hila::output << "Actual deriv " << (s2-s1)/eps << "\n";
    hila::output << "deriv diff " << diff << "\n";
    assert( diff*diff < 10*eps && "Fermion deriv" );
  }


  sigma_mom[ALL] = 0; pi_mom[ALL] = 0;
  sigma[ALL] = gaussian_ran(); pi[ALL] = gaussian_ran();
  foralldir(dir){
    momentum[dir][ALL] = 0;
    gauge[dir][ALL] = 1;
  }

  s = pi.get_value_at(50);

  if(mynode()==0){
    sigma.set_value_at(s+eps,50);
  }
  sigma.mark_changed(ALL);
  s2 = aa.action();

  if(mynode()==0)
    sigma.set_value_at(s, 50);
  sigma.mark_changed(ALL);
  s1 = aa.action();

  aa.force_step(1.0);
  f = sigma_mom.get_value_at(50);
  diff = f - (s2-s1)/eps;

  if(mynode()==0) {
    hila::output << "Calculated force " << f << "\n";
    hila::output << "Actual force " << (s2-s1)/eps << "\n";
    hila::output << "Sigma force diff " << diff << "\n";
    assert( diff*diff < 10*eps && "Fermion force" );
  }


  fa.draw_gaussian_fields();
  foralldir(dir){
    momentum[dir][ALL] = 0;
    gauge[dir][ALL] = 1;
  }
  sigma_mom[ALL] = 0; pi_mom[ALL] = 0;
  sigma[ALL] = 0; pi[ALL] = 0;
  s = pi.get_value_at(50);

  if(mynode()==0){
    sigma.set_value_at(s+eps,50);
  }
  sigma.mark_changed(ALL);
  s2 = fa.action();

  if(mynode()==0)
    sigma.set_value_at(s, 50);
  sigma.mark_changed(ALL);
  s1 = fa.action();

  fa.force_step(1.0);
  f = sigma_mom.get_value_at(50);
  diff = f - (s2-s1)/eps;

  if(mynode()==0) {
    hila::output << "Calculated force " << f << "\n";
    hila::output << "Actual force " << (s2-s1)/eps << "\n";
    hila::output << "Sigma force diff " << diff << "\n";
    assert( diff*diff < 10*eps && "Fermion force" );
  }



  sigma_mom[ALL] = 0; pi_mom[ALL] = 0;
  sigma[ALL] = 0; pi[ALL] = 0;
  p = pi.get_value_at(50);
  
  if(mynode()==0){
    pi.set_value_at(p+eps,50);
  }
  pi.mark_changed(ALL);
  s2 = fa.action();

  if(mynode()==0)
    pi.set_value_at(p, 50);
  pi.mark_changed(ALL);
  s1 = fa.action();

  fa.force_step(1.0);
  f = pi_mom.get_value_at(50);
  diff = f - (s2-s1)/eps;

  if(mynode()==0) {
    hila::output << "Calculated force " << f << "\n";
    hila::output << "Actual force " << (s2-s1)/eps << "\n";
    hila::output << "Pi force diff " << diff << "\n";
    assert( diff*diff < 10*eps && "Fermion force" );
  }
  output0 << "\n";




  sigma[ALL] = 0; pi[ALL] = 0;
  sigma_mom[ALL] = 0; pi_mom[ALL] = 0;
  foralldir(dir){
    momentum[dir][ALL] = 0;
    gauge[dir][ALL] = 1;
  }


  // Run HMC using the integrator
  for(int step = 0; step < 100; step ++){
    // Run update
    update_hmc(integrator_level_2, hmc_steps, traj_length);
    measure(gauge, sigma, pi);
  }



  finishrun();

  return 0;
}
