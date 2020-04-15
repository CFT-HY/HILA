/*********************************************
 * Simulates staggered fermions with a gauge *
 * interaction                               *
 *********************************************/


#include "staggered.h"
#include "../plumbing/inputs.h"
#include "../plumbing/algorithms/hmc.h"
#include "../plumbing/gauge_field.h"













class fermion_term{
  public:
    gauge_term<SUN,NMAT> gt;
    field<SUN> *gauge;

    fermion_term(gauge_term<SUN,NMAT> g) : gt(g) {
      gauge = g.gauge;
    }

    //The gauge action
    double action(){
      return gt.action();
    }

    /// Gaussian random momentum for each element
    void generate_momentum(){
      gt.generate_momentum();
    }

    // Update the momentum with the gauge field
    void force_step(double eps){
      // Fermion force here
     }

    // Update the gauge field with momentum
    void momentum_step(double eps){
      gt.integrator_step(eps);
    }

    // A single gauge update
    void integrator_step(double eps){
      O2_step(*this, eps);
    }
};


class full_action{
  public:
    fermion_term ft;
    field<SUN> *gauge;

    full_action(fermion_term f) : ft(f) {gauge = f.gauge;}

    //The gauge action
    double action(){
      return ft.action();
    }

    /// Gaussian random momentum for each element
    void generate_momentum(){
      ft.generate_momentum();
    }

    // Update the momentum with the gauge field
    void integrate(int steps, double dt){
      for(int step=0; step < steps; step++){
        ft.integrator_step(dt/steps);
      }
    }
};










int main(int argc, char **argv){

  input parameters = input();
  parameters.import("parameters");
  double beta = parameters.get("beta");
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

  gauge_term<SUN,NMAT> gt(gauge, momentum, 1.0);

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
  gt.generate_momentum();

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
  gt = gauge_term<SUN,NMAT>(gauge, momentum, beta);
  fermion_term ft = fermion_term(gt);
  full_action action = full_action(ft);
  for(int step = 0; step < 100; step ++){
    update_hmc(action, hmc_steps, traj_length);
    double plaq = plaquette(gauge);
    output0 << "Plaq: " << plaq << "\n";
  }



  finishrun();

  return 0;
}
