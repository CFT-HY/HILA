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


using mytype = SU<30,double>;
//using mytype = double;
int main(int argc, char **argv) {

  hila::initialize(argc, argv);

  hila::input par("parameters");

  CoordinateVector lsize = par.get("lattice size"); // reads NDIM numbers
  int iterations = par.get("iterations");
  
  par.close(); // file is closed also when par goes out of scope

  
  lattice.setup(lsize);

  hila::timer reduction_timer("Reduction timer");
  hila::seed_random(-1);
  Field<mytype> f,a;
  a.random();
  onsites(ALL) f[X] = X.z()+1;
  onsites(ALL) a[X] = a[X] + a[X];

  reduction_timer.start();
  ReductionVector<mytype> reduction_vec(lattice.size(e_z));
  for (int i = 0; i <= iterations; i++) {
    
    reduction_vec = 0;
    reduction_vec.allreduce(false);

    onsites(ALL) {
      mytype val = f[X];
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
