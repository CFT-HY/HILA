#include "test.h"
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

  field<matrix<1, N, cmplx<double>> > a, b, Db, Ddaggera, DdaggerDb;
  field<matrix<1, N, cmplx<double>> > sol;
  field<matrix<N,N, cmplx<double>> > U[NDIM];

  onsites(ALL){
    for(int i=0; i<N; i++){
      a[X].random();
      b[X].random();
      sol[X].c[0][i] = 0;
    }
  }

  foralldir(d) {
    onsites(ALL){
      U[d][X].random();
    }
  }

  // Check conjugate of the Dirac operator
  double diffre = 0, diffim = 0;
  dirac_stagggered(U, 0.1, b, Db);
  dirac_stagggered_dagger(U, 0.1, a, Ddaggera);
  onsites(ALL){
    diffre += (a[X]*Db[X]).re - (Ddaggera[X]*b[X]).re;
    diffim += (a[X]*Db[X]).im - (Ddaggera[X]*b[X]).im;
  }

  assert(diffre*diffre < 1e-16 && "test dirac_stagggered_dagger");
  assert(diffim*diffim < 1e-16 && "test dirac_stagggered_dagger");
  
  // Now run CG on DdaggerDb and check the result is b
  dirac_stagggered_dagger(U, 0.1, Db, DdaggerDb);
  CG_engine<staggered_dirac> engine;
  engine.solve(U, 0.1, DdaggerDb, a);

  onsites(ALL){
    diffre += norm_squared(a[X]-b[X]);
  }
  printf(" %g \n", diffre);
  assert(diffre*diffre < 1e-8 && "test CG");


  finishrun();
}
