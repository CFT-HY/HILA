#include "test.h"


#include "../plumbing/fermion/staggered.h"
#include "../plumbing/fermion/wilson.h"
#include "../plumbing/algorithms/conjugate_gradients.h"


#define N 3




void test_gamma_matrices(){
  Wilson_vector<N> w1, w2, w3;
  half_Wilson_vector<N> h1;
  SU<N> U; U.random();
  w1.gaussian();

#if NDIM == 4
  w2 = w1-gamma5*(gamma5*w1);
  assert(w2.norm_sq() < 0.0001 && "g5*g5 = 1");
#endif

  foralldir(d){
    w2 = w1-gamma_matrix[d]*(gamma_matrix[d]*w1);
    assert(w2.norm_sq() < 0.0001 && "gamma_d*gamma_d = 1");

    w2 = w1 + gamma_matrix[d]*w1;
    h1 = half_Wilson_vector(w1,d);
    double diff = w2.norm_sq() - h1.norm_sq();
    assert(diff*diff < 0.0001 && "half_Wilson_vector projection norm to direction XUP");

    w3 = (U*h1).expand(d) - U*w2;
    assert(w3.norm_sq() < 0.0001 && "half_wilson_vector expand");


    w2 = w1 - gamma_matrix[d]*w1;
    h1 = half_Wilson_vector(w1,opp_dir(d));
    diff = w2.norm_sq() - h1.norm_sq();
    assert(diff*diff < 0.0001 && "half_Wilson_vector projection norm to direction XUP");

    w3 = (U*h1).expand(opp_dir(d)) - U*w2;
    assert(w3.norm_sq() < 0.0001 && "half_wilson_vector expand");

  }

}




int main(int argc, char **argv){

  #if NDIM==1
  lattice->setup( 64, argc, argv );
  #elif NDIM==2
  lattice->setup( 32, 8, argc, argv );
  #elif NDIM==3
  lattice->setup( 16, 8, 8, argc, argv );
  #elif NDIM==4
  lattice->setup( 16, 8, 8, 8, argc, argv );
  #endif
  seed_random(2);

  test_gamma_matrices();

  field<SU<N>> U[NDIM];
  foralldir(d) {
    onsites(ALL){
      U[d][X].random();
    }
  }

  // Check conjugate of the staggered Dirac operator
  {
    using dirac = dirac_staggered< SU_vector<N>, SU<N>>;
    dirac D(0.1, U);
    field<SU_vector<N>> a, b, Db, Ddaggera, DdaggerDb;
    field<SU_vector<N>> sol;
    onsites(ALL){
      a[X].gaussian();
      b[X].gaussian();
      sol[X] = 0;
    }

    double diffre = 0, diffim = 0;
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
  }

  // Check conjugate of the wilson Dirac operator
  {
    using dirac = dirac_wilson<N>;
    dirac D(0.1, U);
    field<Wilson_vector<N>> a, b, Db, Ddaggera, DdaggerDb;
    field<Wilson_vector<N>> sol;
    onsites(ALL){
      a[X].gaussian();
      b[X].gaussian();
      sol[X] = 0;
    }

    double diffre = 0, diffim = 0;
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
    CG<field<Wilson_vector<N>>, dirac> inverse(D);
    inverse.apply(DdaggerDb, a);

    onsites(ALL){
      diffre += norm_squared(a[X]-b[X]);
    }
    assert(diffre*diffre < 1e-8 && "test CG");
  }


  finishrun();
}


