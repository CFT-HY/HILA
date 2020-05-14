#include "test.h"
#include "../datatypes/sun_vector.h"
#include "../datatypes/sun.h"
#include "../plumbing/algorithms/conjugate_gradients.h"

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

  field<SU_vector<N>> a, b, Db, Ddaggera, DdaggerDb;
  field<SU_vector<N>> sol;
  field<SU<N>> U[NDIM];

  onsites(ALL){
    a[X].gaussian();
    b[X].gaussian();
    sol[X] = 0;
  }

  foralldir(d) {
    onsites(ALL){
      U[d][X].random();
    }
  }

  // Check conjugate of the staggered Dirac operator
  double diffre = 0, diffim = 0;
  using dirac = dirac_staggered< SU_vector<N>, SU<N>>;
  dirac D(0.1, U);
  D.apply(b, Db);
  D.dagger(a, Ddaggera);
  onsites(ALL){
    diffre += a[X].dot(Db[X]).re - Ddaggera[X].dot(b[X]).re;
    diffim += a[X].dot(Db[X]).im - Ddaggera[X].dot(b[X]).im;
  }

  assert(diffre*diffre < 1e-16 && "test dirac_stagggered_dagger");
  assert(diffim*diffim < 1e-16 && "test dirac_stagggered_dagger");
  
  // Now run CG on DdaggerDb and check the result is b
  D.dagger(Db, DdaggerDb);
  CG<field<SU_vector<N>>, dirac> inverse(D);
  inverse.apply(DdaggerDb, a);

  onsites(ALL){
    diffre += norm_squared(a[X]-b[X]);
  }
  assert(diffre*diffre < 1e-8 && "test CG");


  finishrun();
}
