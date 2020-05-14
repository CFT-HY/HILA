#include "test.h"
#include "../datatypes/sun_vector.h"
#include "../datatypes/sun.h"
#include "../plumbing/field.h"
#include "../plumbing/fermion/staggered.h"
#include "../plumbing/fermion/wilson.h"
#include "../plumbing/algorithms/hmc.h"
#include "../plumbing/gauge_field.h"
#include "../plumbing/fermion/fermion_field.h"

#define N 3

int main(int argc, char **argv){

  /* Use a smaller lattice size since
   * the inversion takes a while */
  #if NDIM==1
  lattice->setup( 8, argc, argv );
  #elif NDIM==2
  lattice->setup( 8, 8, argc, argv );
  #elif NDIM==3
  lattice->setup( 8, 8, 8, argc, argv );
  #elif NDIM==4
  lattice->setup( 8, 8, 8, 8, argc, argv );
  #endif
  seed_random(2);

  // Test the force calculation by varying one gauge link
  // (Needs to be moved to tests)
  field<SU<N>> gauge[NDIM];
  field<SU<N>> momentum[NDIM];
  double eps = 1e-5;

  gauge_action<N> ga(gauge, momentum, 1.0);
  foralldir(dir){
    onsites(ALL){
      gauge[dir][X].random();
    }
  }

  
  for(int ng = 0; ng < ga.n_generators(); ng++){
    foralldir(dir){
      onsites(ALL){
        momentum[dir][X] = 0;
      }
    }
    SU<N> g1 = gauge[0].get_value_at(50);
    SU<N> h = 1;
    h += eps * ga.generator(ng);
    SU<N> g12 = h*g1;

    double s1 = plaquette_sum(gauge);

    if(mynode()==0)
      gauge[0].set_value_at(g12,50);
    gauge[0].mark_changed(ALL);

    double s2 = plaquette_sum(gauge);

    if(mynode()==0)
      gauge[0].set_value_at(g1, 50);
    gauge[0].mark_changed(ALL);

    gauge_force(gauge, momentum, 1.0/N);
    SU<N> f = momentum[0].get_value_at(50);
    double diff = (f*ga.generator(ng)).trace().re - (s2-s1)/eps;

    if(mynode()==0) {
      //hila::output << "Force 1 " << (f*ga.generator(ng)).trace().re << "\n";
      //hila::output << "Force 2 " << (s2-s1)/eps << "\n";
      //hila::output << "Force " << ng << " diff " << diff << "\n";
      h = ga.generator(ng);
      assert( diff*diff < eps*100 && "Gauge force" );
    }
  }


  for(int ng = 0; ng < ga.n_generators(); ng++){
    dirac_staggered<SU_vector<N>, SU<N>> D(1.0, gauge);
    fermion_action fa(D, gauge, momentum);
    fa.draw_gaussian_fields();
    foralldir(dir){
      onsites(ALL){
        gauge[dir][X].random();
        momentum[dir][X] = 0;
      }
    }

    SU<N> g1 = gauge[0].get_value_at(50);
    SU<N> h = SU<N>(1) + eps * ga.generator(ng);
    SU<N> g12 = h*g1;

    static field<SU_vector<N>> psi, chi, tmp, tmp2;
    onsites(ALL){
      psi[X].gaussian();
      chi[X].gaussian();
    }
    
    double s1 = 0;
    D.apply(psi,tmp);
    onsites(ALL){
      s1 += chi[X].rdot(tmp[X]);
    }

    if(mynode()==0){
      gauge[0].set_value_at(g12,50);
    }
    gauge[0].mark_changed(ALL);

    double s2 = 0;
    D.apply(psi,tmp);
    onsites(ALL){
      s2 += chi[X].rdot(tmp[X]);
    }

    if(mynode()==0)
      gauge[0].set_value_at(g1, 50);
    gauge[0].mark_changed(ALL);

    D.force(chi, psi, momentum);
    SU<N> f = momentum[0].get_value_at(50);
    double diff = (f*ga.generator(ng)).trace().re - (s2-s1)/eps;

    if(mynode()==0) {
      //hila::output << "Calculated deriv " << (f*ga.generator(ng)).trace().re << "\n";
      //hila::output << "Actual deriv " << (s2-s1)/eps << "\n";
      //hila::output << "Staggered deriv " << ng << " diff " << diff << "\n";
      assert( diff*diff < eps*100 && "Staggered fermion deriv" );
    }


    foralldir(dir){
      momentum[dir][ALL] = 0;
    }

    if(mynode()==0){
      gauge[0].set_value_at(g12,50);
    }
    gauge[0].mark_changed(ALL);
    s2 = fa.action();

    if(mynode()==0)
      gauge[0].set_value_at(g1, 50);
    gauge[0].mark_changed(ALL);
    s1 = fa.action();

    fa.force_step(1.0);
    f = momentum[0].get_value_at(50);
    diff = (f*ga.generator(ng)).trace().re - (s2-s1)/eps;

    if(mynode()==0) {
      //hila::output << "Calculated force " << (f*ga.generator(ng)).trace().re << "\n";
      //hila::output << "Actual force " << (s2-s1)/eps << "\n";
      //hila::output << "Staggered force " << ng << " diff " << diff << "\n";
      assert( diff*diff < eps*100 && "Staggered force" );
    }
  }


  for(int ng = 0; ng < ga.n_generators(); ng++){
    dirac_wilson<N> D(0.05, gauge);
    fermion_action fa(D, gauge, momentum);
    fa.draw_gaussian_fields();
    foralldir(dir){
      onsites(ALL){
        gauge[dir][X].random();
        momentum[dir][X] = 0;
      }
    }

    SU<N> g1 = gauge[0].get_value_at(50);
    SU<N> h = SU<N>(1) + eps * ga.generator(ng);
    SU<N> g12 = h*g1;

    static field<Wilson_vector<N>> psi, chi, tmp, tmp2;
    onsites(ALL){
      psi[X].gaussian();
      chi[X].gaussian();
    }
    double s1 = 0;
    D.apply(psi,tmp);
    onsites(ALL){
      s1 += chi[X].rdot(tmp[X]);
    }

    if(mynode()==0){
      gauge[0].set_value_at(g12,50);
    }
    gauge[0].mark_changed(ALL);
    double s2 = 0;
    D.apply(psi,tmp);
    onsites(ALL){
      s2 += chi[X].rdot(tmp[X]);
    }

    if(mynode()==0)
      gauge[0].set_value_at(g1, 50);
    gauge[0].mark_changed(ALL);

    D.force(chi, psi, momentum, 1);
    SU<N> f = momentum[0].get_value_at(50);
    double diff = (f*ga.generator(ng)).trace().re - (s2-s1)/eps;

    if(mynode()==0) {
      //hila::output << "Calculated deriv " << (f*ga.generator(ng)).trace().re << "\n";
      //hila::output << "Actual deriv " << (s2-s1)/eps << "\n";
      //hila::output << "Wilson deriv " << ng << " diff " << diff << "\n";
      assert( diff*diff < eps*100 && "Wilson fermion deriv" );
    }

    foralldir(dir){
      momentum[dir][ALL] = 0;
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

    D.force(chi, psi, momentum, -1);
    f = momentum[0].get_value_at(50);
    diff = (f*ga.generator(ng)).trace().re - (s2-s1)/eps;

    if(mynode()==0) {
      //hila::output << "Calculated deriv " << (f*ga.generator(ng)).trace().re << "\n";
      //hila::output << "Actual deriv " << (s2-s1)/eps << "\n";
      //hila::output << "Wilson deriv " << ng << " diff " << diff << "\n";
      assert( diff*diff < eps*100 && "Wilson fermion deriv" );
    }


    foralldir(dir){
      momentum[dir][ALL] = 0;
    }

    if(mynode()==0)
      gauge[0].set_value_at(g12,50);
    gauge[0].mark_changed(ALL);
    s2 = fa.action();

    if(mynode()==0)
      gauge[0].set_value_at(g1, 50);
    gauge[0].mark_changed(ALL);
    s1 = fa.action();

    fa.force_step(1.0);
    f = momentum[0].get_value_at(50);
    diff = (f*ga.generator(ng)).trace().re - (s2-s1)/eps;

    if(mynode()==0) {
      //hila::output << "Calculated force " << (f*ga.generator(ng)).trace().re << "\n";
      //hila::output << "Actual force " << (s2-s1)/eps << "\n";
      //hila::output << "Wilson force " << ng << " diff " << diff << "\n";
      assert( diff*diff < eps*100 && "Wilson fermion force" );
    }
  }


  for(int ng = 0; ng < ga.n_generators(); ng++){
    // Check also the momentum action and derivative
    ga.draw_gaussian_fields();

    double s1 = momentum_action(momentum);
    SU<N> h = momentum[0].get_value_at(0);
    h += eps * ga.generator(ng);
    if(mynode()==0)
      momentum[0].set_value_at(h, 0);
    double s2 = momentum_action(momentum);

    double diff = (h*ga.generator(ng)).trace().re + (s2-s1)/eps;
    if(mynode()==0) {
      //hila::output << "Mom 1 " << (h*ga.generator(ng)).trace().re << "\n";
      //hila::output << "Mom 2 " << (s2-s1)/eps << "\n";
      //hila::output << "Mom " << ng << " diff " << diff << "\n";
      h = ga.generator(ng);
      assert( diff*diff < eps*100 && "Momentum derivative" );
    }
  }
}

