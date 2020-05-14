/*********************************************
 * Simulates Wilson fermions with a gauge *
 * interaction                               *
 *********************************************/


#include "Wilson.h"




int main(int argc, char **argv){

  input parameters = input();
  parameters.import("parameters");
  double beta = parameters.get("beta");
  double kappa = parameters.get("kappa");
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

  ga.set_unity();

  // Define a Dirac operator
  dirac_wilson<N> D(kappa, gauge);
  fermion_action fa(D, gauge, momentum);

  // Check conjugate of the Dirac operator
  field<Wilson_vector<N>> a, b, Db, Ddaggera, DdaggerDb;
  onsites(ALL){
    a[X].gaussian();
    b[X].gaussian();
  }
  double diffre = 0, diffim = 0;
  D.apply(b, Db);
  D.dagger(a, Ddaggera);
  onsites(ALL){
    diffre += a[X].dot(Db[X]).re - Ddaggera[X].dot(b[X]).re;
    diffim += a[X].dot(Db[X]).im - Ddaggera[X].dot(b[X]).im;
  }

  assert(diffre*diffre < 1e-16 && "test Wilson dirac conjugate");
  assert(diffim*diffim < 1e-16 && "test Wilson dirac conjugate");
  

  finishrun();

  return 0;
}
