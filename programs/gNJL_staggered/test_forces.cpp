/*********************************************
 * Simulates staggered fermions with a gauge *
 * interaction                               *
 *********************************************/


#include "staggered.h"





int main(int argc, char **argv){

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




  // Test the derivative of the Dirac operator
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
    //hila::output << "Calculated deriv " << f << "\n";
    //hila::output << "Actual deriv " << (s2-s1)/eps << "\n";
    //hila::output << "deriv diff " << diff << "\n";
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
    //hila::output << "Calculated deriv " << f << "\n";
    //hila::output << "Actual deriv " << (s2-s1)/eps << "\n";
    //hila::output << "deriv diff " << diff << "\n";
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
    //hila::output << "Calculated deriv " << f << "\n";
    //hila::output << "Actual deriv " << (s2-s1)/eps << "\n";
    //hila::output << "deriv diff " << diff << "\n";
    assert( diff*diff < 10*eps && "Fermion deriv" );
  }


  // Test the derivatives of the auxfield action
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
    //hila::output << "Calculated force " << f << "\n";
    //hila::output << "Actual force " << (s2-s1)/eps << "\n";
    //hila::output << "Sigma force diff " << diff << "\n";
    assert( diff*diff < 10*eps && "Fermion force" );
  }


  // Test the derivatives of the fermion action
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
    //hila::output << "Calculated force " << f << "\n";
    //hila::output << "Actual force " << (s2-s1)/eps << "\n";
    //hila::output << "Sigma fermion force diff " << diff << "\n";
    assert( diff*diff < 10*eps && "Fermion force" );
  }



  fa.draw_gaussian_fields();
  sigma_mom[ALL] = 0; pi_mom[ALL] = 0;
  sigma[ALL] = 0; pi[ALL] = 0;
  foralldir(dir){
    momentum[dir][ALL] = 0;
    gauge[dir][ALL] = 1;
  }
  p = pi.get_value_at(0);
  
  if(mynode()==0){
    pi.set_value_at(p+eps,0);
  }
  pi.mark_changed(ALL);
  s2 = fa.action();

  if(mynode()==0)
    pi.set_value_at(p, 0);
  pi.mark_changed(ALL);
  s1 = fa.action();

  fa.force_step(1.0);
  if(mynode()==0) {
    f = pi_mom.get_value_at(0);
    diff = f - (s2-s1)/eps;

    //hila::output << "Action  " << s1 << "\n";
    //hila::output << "Action  " << s2 << "\n";
    //hila::output << "Calculated force " << f << "\n";
    //hila::output << "Actual force " << (s2-s1)/eps << "\n";
    //hila::output << "Pi force diff " << diff << "\n";
    assert( diff*diff < 10*eps && "Fermion force" );
  }
  output0 << "\n";


  finishrun();

  return 0;
}
