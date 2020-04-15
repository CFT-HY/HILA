/*********************************************
 * Simulates staggered fermions with a gauge *
 * interaction                               *
 *********************************************/


#include "staggered.h"
#include "../plumbing/inputs.h"
#include "../plumbing/algorithms/hmc.h"
#include "../plumbing/gauge_field.h"
#include "../plumbing/fermion_field.h"




using SG = gauge_action<SUN,NMAT>;
using Dtype = dirac_staggered<field<VEC>,field<SUN>>;
using SF = fermion_action<SG, VEC, Dtype>;





class full_action{
  public:
    SF fa;
    field<SUN> *gauge;

    full_action(SF f) : fa(f) {gauge = f.gauge;}

    //The gauge action
    double action(){
      return fa.action();
    }

    /// Gaussian random momentum for each element
    void generate_momentum(){
      fa.generate_momentum();
    }

    // Update the momentum with the gauge field
    void integrate(int steps, double dt){
      for(int step=0; step < steps; step++){
        fa.integrator_step(dt/steps);
      }
    }
};










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


  field<SUN> gauge[NDIM];

  // Starting with a unit configuration
  foralldir(dir){
    onsites(ALL){
      gauge[dir][X].random();
    }
  }


  field<NMAT> momentum[NDIM];
  foralldir(dir){
    onsites(ALL){
      momentum[dir][X] = 0;
    }
  }


  // Test the force calculation by varying one gauge link
  // (Needs to be moved to tests)
  double eps = 1e-6;
  SUN g1 = gauge[0].get_value_at(50);
  SUN h = 1;
  h.c[0][1].re += eps;
  h.c[1][0].re -= eps;
  SUN g12 = h*g1;

  gauge_action<SUN,NMAT> ga(gauge, momentum, 1.0);

  if(mynode()==0)
    gauge[0].set_value_at(g1, 50);
  gauge[0].mark_changed(ALL);
  double s1 = plaquette_sum(gauge);

  if(mynode()==0)
    gauge[0].set_value_at(g12,50);
  gauge[0].mark_changed(ALL);
  double s2 = plaquette_sum(gauge);

  if(mynode()==0)
    gauge[0].set_value_at(g1, 50);
  gauge[0].mark_changed(ALL);

  gauge_force(gauge, momentum, 1.0/N);
  NMAT f = momentum[0].get_value_at(50);
  double diff = 2*f.c[0][1].re + (s2-s1)/eps;
  if(mynode()==0) 
    assert( diff*diff < eps*eps*100 );


  // Check also the momentum action and derivative
  ga.generate_momentum();

  s1 = momentum_action(momentum);
  h = momentum[0].get_value_at(0);
  h.c[0][0].im += eps;
  if(mynode()==0)
    momentum[0].set_value_at(h, 0);
  s2 = momentum_action(momentum);

  diff = h.c[0][0].im - (s2-s1)/eps;
  if(mynode()==0) 
    assert( diff*diff < eps*eps*100 );




  // Now the actual simulation
  ga = SG(gauge, momentum, beta);
  Dtype D(mass, gauge);
  SF fa(ga, D);
  full_action action = full_action(fa);
  for(int step = 0; step < 100; step ++){
    update_hmc(action, hmc_steps, traj_length);
    double plaq = plaquette(gauge);
    output0 << "Plaq: " << plaq << "\n";
  }



  finishrun();

  return 0;
}
