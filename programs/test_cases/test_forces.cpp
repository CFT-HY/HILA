#include "test.h"
#include "../datatypes/sun_vector.h"
#include "../datatypes/sun.h"
#include "../plumbing/field.h"
#include "../plumbing/dirac.h"
#include "../plumbing/algorithms/hmc.h"
#include "../plumbing/gauge_field.h"
#include "../plumbing/fermion_field.h"

#define N 3

int main(int argc, char **argv){

  #if NDIM==1
  lattice->setup( 16, argc, argv );
  #elif NDIM==2
  lattice->setup( 8, 8, argc, argv );
  #elif NDIM==3
  lattice->setup( 8, 8, 8, argc, argv );
  #elif NDIM==4
  lattice->setup( 8, 8, 8, 6, argc, argv );
  #endif
  seed_random(2);

  // Test the force calculation by varying one gauge link
  // (Needs to be moved to tests)
  field<SU<N>> gauge[NDIM];
  field<SU<N>> momentum[NDIM];

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

    double eps = 1e-4;
    SU<N> g1 = gauge[0].get_value_at(50);
    SU<N> h = 1;
    h += eps * ga.generator(ng);
    SU<N> g12 = h*g1;

    double s1 = plaquette_sum(gauge);

    if(mynode()==0){
      gauge[0].set_value_at(g12,50);
      //gauge[0].set_value_at(g12,lattice->neighb[1][50]);
    }
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
      assert( diff*diff < eps*eps*1000 && "Gauge force" );
    }


    dirac_staggered<SU_vector<N>, SU<N>> D(1.0, gauge);
    fermion_action fa(D, gauge, momentum);
    for(int ng = 0; ng < ga.n_generators(); ng++){
      fa.draw_gaussian_fields();
      foralldir(dir){
        onsites(ALL){
          gauge[dir][X].random();
        }
      }
      foralldir(dir){
        onsites(ALL){
          momentum[dir][X] = 0;
        }
      }

      g1 = gauge[0].get_value_at(50);
      h = 1;
      h += eps * ga.generator(ng);
      g12 = h*g1;


      static field<SU_vector<N>> psi, chi, tmp, tmp2;
      onsites(ALL){
        psi[X].gaussian();
        chi[X].gaussian();
      }
      s1 = 0;
      D.apply(psi,tmp);
      onsites(ALL){
        s1 += (chi[X]*tmp[X]).re;
      }

      if(mynode()==0){
        gauge[0].set_value_at(g12,50);
      }
      gauge[0].mark_changed(ALL);
      s2 = 0;
      D.apply(psi,tmp);
      onsites(ALL){
        s2 += (chi[X]*tmp[X]).re;
      }

      if(mynode()==0)
        gauge[0].set_value_at(g1, 50);
      gauge[0].mark_changed(ALL);

      D.force(chi, psi, momentum);
      f = momentum[0].get_value_at(50);
      diff = (f*ga.generator(ng)).trace().re - (s2-s1)/eps;

      if(mynode()==0) {
        //hila::output << "Action 1 " << s1 << "\n";
        //hila::output << "Action 2 " << s2 << "\n";
        //hila::output << "Calculated deriv " << (f*ga.generator(ng)).trace().re << "\n";
        //hila::output << "Actual deriv " << (s2-s1)/eps << "\n";
        //hila::output << "deriv " << ng << " diff " << diff << "\n";
        assert( diff*diff < eps*eps*1000 && "Fermion deriv" );
      }


      foralldir(dir){
        onsites(ALL){
          momentum[dir][X] = 0;
        }
      }

      s1 = fa.action();

      if(mynode()==0){
        gauge[0].set_value_at(g12,50);
      }
      gauge[0].mark_changed(ALL);
      s2 = fa.action();

      if(mynode()==0)
        gauge[0].set_value_at(g1, 50);
      gauge[0].mark_changed(ALL);

      fa.force_step(1.0);
      f = momentum[0].get_value_at(50);
      diff = (f*ga.generator(ng)).trace().re - (s2-s1)/eps;

      if(mynode()==0) {
        //hila::output << "Action 1 " << s1 << "\n";
        //hila::output << "Action 2 " << s2 << "\n";
        //hila::output << "Calculated force " << (f*ga.generator(ng)).trace().re << "\n";
        //hila::output << "Actual force " << (s2-s1)/eps << "\n";
        //hila::output << "Force " << ng << " diff " << diff << "\n";
        assert( diff*diff < eps*eps*1000 && "Fermion force" );
      }
    }



    // Check also the momentum action and derivative
    ga.draw_gaussian_fields();

    s1 = momentum_action(momentum);
    h = momentum[0].get_value_at(0);
    h += eps * ga.generator(ng);
    if(mynode()==0)
      momentum[0].set_value_at(h, 0);
    s2 = momentum_action(momentum);

    diff = (h*ga.generator(ng)).trace().re + (s2-s1)/eps;
    if(mynode()==0) {
      //hila::output << "Mom 1 " << (h*ga.generator(ng)).trace().re << "\n";
      //hila::output << "Mom 2 " << (s2-s1)/eps << "\n";
      //hila::output << "Mom " << ng << " diff " << diff << "\n";
      h = ga.generator(ng);
      assert( diff*diff < eps*eps*1000 && "Momentum derivative" );
    }
  }
}

