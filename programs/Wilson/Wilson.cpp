/*********************************************
 * Simulates Wilson fermions with a gauge *
 * interaction                               *
 *********************************************/


#include "Wilson.h"


void test_gamma_matrices(){
  Wilson_vector<N> w1, w2;
  half_Wilson_vector<N> h1, h2;
  w1.gaussian();

  w2 = w1-gamma5*(gamma5*w1);
  assert(w2.norm_sq() < 0.0001 && "g5*g5 = 1");

  w2 = w1-gamma0*(gamma0*w1);
  assert(w2.norm_sq() < 0.0001 && "g0*g0 = 1");

  w2 = w1-gamma1*(gamma1*w1);
  assert(w2.norm_sq() < 0.0001 && "g1*g1 = 1");

  w2 = w1-gamma2*(gamma2*w1);
  assert(w2.norm_sq() < 0.0001 && "g2*g2 = 1");

  w2 = w1-gamma3*(gamma3*w1);
  assert(w2.norm_sq() < 0.0001 && "g3*g3 = 1");

  foralldir(d){
    w2 = w1 + gamma_matrix[d]*w1;
    h2 = half_Wilson_vector(w1,d);
    double diff = w2.norm_sq() - h2.norm_sq();
    assert(diff*diff < 0.0001 && "half_Wilson_vector projection norm to direction XUP");

    w2 = w1 + gamma_matrix[d]*w1;
    h2 = half_Wilson_vector(w1,d);
    diff = w2.norm_sq() - h2.norm_sq();
    assert(diff*diff < 0.0001 && "half_Wilson_vector projection norm to direction XUP");

    //output0 << w2.str();
    //output0 << h2.str();
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

  // Define gauge field and momentum field
  field<SUN> gauge[NDIM];
  field<SUN> momentum[NDIM];

  // Define gauge and momentum action terms
  gauge_momentum_action ma(gauge, momentum);
  gauge_action ga(gauge, momentum, beta);
  
  test_gamma_matrices();

  finishrun();

  return 0;
}
