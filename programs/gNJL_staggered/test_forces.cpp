/*********************************************
 * Simulates staggered fermions with a gauge *
 * interaction                               *
 *********************************************/


#include "staggered.h"





int main(int argc, char **argv){

  lattice->setup( nd[0], nd[1], nd[2], nd[3], argc, argv );
  seed_random(12345);

  // Define gauge field and momentum field
  field<SUN> gauge[NDIM];
  field<SUN> momentum[NDIM];

  // Define gauge and momentum action terms
  gauge_momentum_action gma(gauge, momentum);
  gauge_action ga(gauge, momentum, 8);

  // Initialize the gauge field
  ga.set_unity();

  // Define auxiliary fields and momentum
  field<double> sigma, pi, sigma_mom, pi_mom;

  // And action
  auxiliary_momentum_action ama(sigma, pi, sigma_mom, pi_mom);
  auxiliary_action aa(sigma, pi, sigma_mom, pi_mom, 0.375);

  // Initialize
  sigma[ALL] = 0.1;
  pi[ALL] = 0;

  // Define a Dirac operator (2 flavors)
  dirac_staggered_gNJL<VEC,SUN> D(1.5, gauge, sigma, pi);
  gNJL_fermion_action fa(D, momentum, sigma_mom, pi_mom);

  for(int ng = 0; ng < SUN::generator_count(); ng++){
    // Test the derivative of the Dirac operator
    foralldir(dir){
      momentum[dir][ALL] = 0;
    }
    sigma_mom[ALL] = 0;
    pi_mom[ALL] = 0;

    double eps = 1e-5;
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

    D.force(chi, psi, momentum, sigma_mom, pi_mom, 1);
    double f = sigma_mom.get_value_at(50);
    double diff = f - (s2-s1)/eps;

    if(mynode()==0) {
      hila::output << "Calculated deriv " << f << "\n";
      hila::output << "Actual deriv " << (s2-s1)/eps << "\n";
      hila::output << "Sigma deriv diff " << diff << "\n";
      assert( diff*diff < 10*eps && "Sigma deriv" );
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
      hila::output << "Pi deriv diff " << diff << "\n";
      assert( diff*diff < 10*eps && "Pi deriv" );
    }

    SUN g1 = gauge[0].get_value_at(50);
    SUN h = SUN(1) + eps * SUN::generator(ng);
    SUN g12 = h*g1;

    s1 = 0;
    D.apply(psi,tmp);
    onsites(ALL){
      s1 += chi[X].rdot(tmp[X]);
    }

    if(mynode()==0){
      gauge[0].set_value_at(g12,50);
    }
    gauge[0].mark_changed(ALL);

    s2 = 0;
    D.apply(psi,tmp);
    onsites(ALL){
      s2 += chi[X].rdot(tmp[X]);
    }

    if(mynode()==0)
      gauge[0].set_value_at(g1, 50);
    gauge[0].mark_changed(ALL);

    SUN m = momentum[0].get_value_at(50);
    double f1 = (s2-s1)/eps;
    double f2 = (m*SUN::generator(ng)).trace().re;
    diff = (f2-f1)/(f1+f2);

    if(mynode()==0) {
      hila::output << "Action 1 " << s1 << "\n";
      hila::output << "Action 2 " << s2 << "\n";
      hila::output << "Calculated deriv " << (m*SUN::generator(ng)).trace().re << "\n";
      hila::output << "Actual deriv " << (s2-s1)/eps << "\n";
      hila::output << "Gauge deriv " << ng << " diff " << diff << "\n";
      assert( diff*diff < eps*10 && "Gauge deriv" );
    }



    sigma_mom[ALL] = 0; pi_mom[ALL] = 0;
    foralldir(dir){
      momentum[dir][ALL] = 0;
    }
    
    s1 = 0;
    D.apply(psi,tmp);
    onsites(ALL){
      s1 += chi[X].rdot(tmp[X]);
    }

    if(mynode()==0){
      sigma.set_value_at(s+eps, 50);
    }
    sigma.mark_changed(ALL);

    s2 = 0;
    D.apply(psi,tmp);
    onsites(ALL){
      s2 += chi[X].rdot(tmp[X]);
    }

    if(mynode()==0){
      sigma.set_value_at(s, 50);
    }
    sigma.mark_changed(ALL);

    D.force(chi, psi, momentum, sigma_mom, pi_mom, -1);
    f = sigma_mom.get_value_at(50);
    diff = f - (s2-s1)/eps;

    if(mynode()==0) {
      hila::output << "Calculated deriv " << f << "\n";
      hila::output << "Actual deriv " << (s2-s1)/eps << "\n";
      hila::output << "Sigma dagger diff " << diff << "\n";
      assert( diff*diff < 10*eps && "Sigma dagger deriv" );
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

    f = pi_mom.get_value_at(50);
    diff = f - (s2-s1)/eps;

    if(mynode()==0) {
      hila::output << "Calculated deriv " << f << "\n";
      hila::output << "Actual deriv " << (s2-s1)/eps << "\n";
      hila::output << "Pi dagger diff " << diff << "\n";
      assert( diff*diff < 10*eps && "Pi dagger deriv" );
    }


    s1 = 0;
    D.dagger(psi,tmp);
    onsites(ALL){
      s1 += chi[X].rdot(tmp[X]);
    }

    if(mynode()==0){
      gauge[0].set_value_at(g12,50);
    }
    gauge[0].mark_changed(ALL);

    s2 = 0;
    D.dagger(psi,tmp);
    onsites(ALL){
      s2 += chi[X].rdot(tmp[X]);
    }

    if(mynode()==0)
      gauge[0].set_value_at(g1, 50);
    gauge[0].mark_changed(ALL);

    m = momentum[0].get_value_at(50);
    f1 = (s2-s1)/eps;
    f2 = (m*SUN::generator(ng)).trace().re;
    diff = (f2-f1)/(f1+f2);

    if(mynode()==0) {
      hila::output << "Action 1 " << s1 << "\n";
      hila::output << "Action 2 " << s2 << "\n";
      hila::output << "Calculated deriv " << (m*SUN::generator(ng)).trace().re << "\n";
      hila::output << "Actual deriv " << (s2-s1)/eps << "\n";
      hila::output << "Gauge dagger deriv " << ng << " diff " << diff << "\n";
      assert( diff*diff < eps*10 && "Gauge dagger deriv" );
    }



    // Test the derivatives of the auxfield action
    sigma_mom[ALL] = 0; pi_mom[ALL] = 0;
    foralldir(dir){
      momentum[dir][ALL] = 0;
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
      hila::output << "Auxfield sigma force diff " << diff << "\n";
      assert( diff*diff < 10*eps && "Auxfield sigma force" );
    }

    // Test the derivatives of the auxfield action
    if(mynode()==0){
      pi.set_value_at(p+eps,50);
    }
    pi.mark_changed(ALL);
    s2 = aa.action();

    if(mynode()==0)
      pi.set_value_at(p, 50);
    pi.mark_changed(ALL);
    s1 = aa.action();

    f = pi_mom.get_value_at(50);
    diff = f - (s2-s1)/eps;

    if(mynode()==0) {
      hila::output << "Calculated force " << f << "\n";
      hila::output << "Actual force " << (s2-s1)/eps << "\n";
      hila::output << "Auxfield pi force diff " << diff << "\n";
      assert( diff*diff < 10*eps && "Auxfield pi force" );
    }



    // Test the derivatives of the fermion action
    fa.draw_gaussian_fields();
    foralldir(dir){
      momentum[dir][ALL] = 0;
      gauge[dir][ALL] = 1;
    }
    sigma_mom[ALL] = 0; pi_mom[ALL] = 0;
    sigma[ALL] = 0.05; pi[ALL] = 0;
    
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
      hila::output << "Sigma fermion force diff " << diff << "\n";
      assert( diff*diff < 10*eps && "Sigma fermion force" );
    }


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

    f = pi_mom.get_value_at(50);
    diff = f - (s2-s1)/eps;

    if(mynode()==0) {

      hila::output << "Action  " << s1 << "\n";
      hila::output << "Action  " << s2 << "\n";
      hila::output << "Calculated force " << f << "\n";
      hila::output << "Actual force " << (s2-s1)/eps << "\n";
      hila::output << "Pi fermion force diff " << diff << "\n";
      assert( diff*diff < 10*eps && "Pi fermion force" );
    }


    s1 = fa.action();

    if(mynode()==0){
      gauge[0].set_value_at(g12,50);
    }
    gauge[0].mark_changed(ALL);

    s2 = fa.action();

    if(mynode()==0){
      gauge[0].set_value_at(g1,50);
    }
    gauge[0].mark_changed(ALL);

    m = momentum[0].get_value_at(50);
    f1 = (s2-s1)/eps;
    f2 = (m*SUN::generator(ng)).trace().re;
    diff = (f2-f1)/(f1+f2);

    if(mynode()==0) {
      hila::output << "Action 1 " << s1 << "\n";
      hila::output << "Action 2 " << s2 << "\n";
      hila::output << "Calculated deriv " << (m*SUN::generator(ng)).trace().re << "\n";
      hila::output << "Actual deriv " << (s2-s1)/eps << "\n";
      hila::output << "Gauge dagger deriv " << ng << " diff " << diff << "\n";
      assert( diff*diff < 10*eps && "Gauge fermion force" );
    }

    output0 << "\n";
  }

  finishrun();

  return 0;
}
