/*********************************************
 * Simulates staggered fermions with a gauge *
 * interaction                               *
 *********************************************/


#include "staggered.h"
#include "../plumbing/inputs.h"
#include "../plumbing/algorithms/hmc.h"
#include "../plumbing/gauge_field.h"
#include "../plumbing/fermion_field.h"






void test_forces(){
  // Test the force calculation by varying one gauge link
  // (Needs to be moved to tests)
  field<SUN> gauge[NDIM];
  field<SUN> momentum[NDIM];

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
    SUN g1 = gauge[0].get_value_at(50);
    SUN h = 1;
    h += eps * ga.generator(ng);
    SUN g12 = h*g1;

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
    SUN f = momentum[0].get_value_at(50);
    double diff = (f*ga.generator(ng)).trace().re - (s2-s1)/eps;

    if(mynode()==0) {
      //hila::output << "Force 1 " << (f*ga.generator(ng)).trace().re << "\n";
      //hila::output << "Force 2 " << (s2-s1)/eps << "\n";
      //hila::output << "Force " << ng << " diff " << diff << "\n";
      h = ga.generator(ng);
      assert( diff*diff < eps*eps*1000 && "Gauge force" );
    }


    dirac_staggered<VEC, SUN> D(1.0, gauge);
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


      static field<VEC> psi, chi, tmp, tmp2;
      onsites(ALL){
        psi[X].random();
        chi[X].random();
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





int main(int argc, char **argv){

  input parameters = input();
  parameters.import("parameters");
  double beta = parameters.get("beta");
  double mass = parameters.get("mass");
  int seed = parameters.get("seed");
	double hmc_steps = parameters.get("hmc_steps");
	double traj_length = parameters.get("traj_length");

  lattice->setup( nd[0], nd[1], nd[2], nd[3], argc, argv );
  seed_random(seed);

  test_forces();

  // Define gauge field and momentum field
  field<SUN> gauge[NDIM];
  field<SUN> momentum[NDIM];

  // Define gauge and momentum action terms
  gauge_momentum_action ma(gauge, momentum);
  gauge_action ga(gauge, momentum, beta);
  
  // Define a Dirac operator
  dirac_staggered<VEC, SUN> D(mass, gauge);
  fermion_action fa(D, gauge, momentum);

  // A second fermion, for checking that addition works
  fermion_action fa2(D, gauge, momentum);

  // Build two integrator levels. Gauge is on the lowest level and
  // the fermions are on higher level
  integrator integrator_level_1(ga, ma);
  integrator integrator_level_2(fa+fa2, integrator_level_1);
  
  // Initialize the gauge field
  ga.set_unity();

  // Run HMC using the integrator
  for(int step = 0; step < 5; step ++){
    update_hmc(integrator_level_2, hmc_steps, traj_length);
    double plaq = plaquette(ga.gauge);
    output0 << "Plaq: " << plaq << "\n";
  }



  finishrun();

  return 0;
}
