#include "gauge/staples.h"
#include "gauge/polyakov.h"
#include "gauge/stout_smear.h"
#include "gauge/sun_heatbath.h"
#include "gauge/sun_overrelax.h"
#include "hila.h"
#include "multicanonical.h"

//#include "gauge/polyakov.h"

#include <numeric>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {

  hila::initialize(argc, argv);

  hila::input par("parameters");

  CoordinateVector lsize = par.get("lattice size"); // reads NDIM numbers
  int iterations = par.get("iterations");
  
  par.close(); // file is closed also when par goes out of scope

  
  lattice.setup(lsize);

  hila::timer reduction_timer("Reduction timer");
  hila::seed_random(-1);
  Field<double> f;
  onsites(ALL) f[X] = X.z()+1;

  reduction_timer.start();
  ReductionVector<double> reduction_vec(lattice.size(e_z));
  for (int i = 0; i <= iterations; i++) {
    
    reduction_vec = 0;
    reduction_vec.allreduce(false);

    onsites(ALL) {
      double val = f[X];
      reduction_vec[X.z()] += val;
    }
  }
  for (int i = 0; i < reduction_vec.size(); ++i) {
    hila::out0 << "reduction_vec[" << i << "] = " << reduction_vec[i] << std::endl;
  }

  hila::synchronize_threads();
  reduction_timer.stop();

  hila::finishrun();
  
}
