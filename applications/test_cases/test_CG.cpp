#include "test.h"


#include "dirac/staggered.h"
#include "dirac/wilson.h"
#include "dirac/Hasenbusch.h"
#include "dirac/conjugate_gradient.h"


#define N 3




void test_gamma_matrices(){
  Wilson_vector<N, double> w1, w2, w3;
  SU<N> U; U.random();
  w1.gaussian();

#if NDIM == 4
  w2 = w1-gamma5*(gamma5*w1);
  assert(w2.norm_sq() < 0.0001 && "g5*g5 = 1");
#endif

  foralldir(d){
    half_Wilson_vector<N, double> h1;
    w2 = w1-gamma_matrix[d]*(gamma_matrix[d]*w1);
    assert(w2.norm_sq() < 0.0001 && "gamma_d*gamma_d = 1");

    w2 = w1 + gamma_matrix[d]*w1;
    h1 = half_Wilson_vector(w1,d,1);
    double diff = w2.norm_sq() - 2*h1.norm_sq();
    assert(diff*diff < 0.0001 && "half_Wilson_vector projection +1 norm");

    w3 = (U*h1).expand(d,1) - U*w2;
    assert(w3.norm_sq() < 0.0001 && "half_wilson_vector expand");


    w2 = w1 - gamma_matrix[d]*w1;
    h1 = half_Wilson_vector(w1,d,-1);
    diff = w2.norm_sq() - 2*h1.norm_sq();
    assert(diff*diff < 0.0001 && "half_Wilson_vector projection -1 norm");

    w3 = (U*h1).expand(d,-1) - U*w2;
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
  lattice->setup( 8, 8, 8, 8, argc, argv );
  #endif
  seed_random(2);


  test_gamma_matrices();

  field<SU<N>> U[NDIM];
  foralldir(d) {
    onsites(ALL){
      if(disable_avx[X]==0){};
      U[d][X].random();
    }
  }


  // Check conjugate of the staggered Dirac operator
  {
    output0 << "Checking with dirac_staggered\n";
    using dirac = dirac_staggered<SU<N>>;
    dirac D(0.1, U);
    field<SU_vector<N, double>> a, b, Db, Ddaggera, DdaggerDb;
    onsites(ALL){
      if(disable_avx[X]==0){};
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
  
    assert(diffre*diffre < 1e-16 && "test dirac_staggered_dagger");
    assert(diffim*diffim < 1e-16 && "test dirac_staggered_dagger");
    
    // Now run CG on Ddaggera. b=1/D a -> Db = a 
    CG<dirac> inverse(D);
    b[ALL] = 0;
    inverse.apply(Ddaggera, b);
    D.apply(b, Db);

    diffre = 0;
    onsites(ALL){
      diffre += norm_squared(a[X]-Db[X]);
    }
    assert(diffre*diffre < 1e-16 && "test D (DdgD)^-1 Ddg");

    // Run CG on a, multiply by DdaggerD and check the same
    onsites(ALL){
      if(disable_avx[X]==0){};
      a[X].gaussian();
    }
    b[ALL] = 0;
    inverse.apply(a, b);
    D.apply(b, Db);
    D.dagger(Db, DdaggerDb);

    diffre = 0;
    onsites(ALL){
      diffre += norm_squared(a[X]-DdaggerDb[X]);
    }
    assert(diffre*diffre < 1e-16 && "test DdgD (DdgD)^-1");

    // The other way around, multiply first
    onsites(ALL){
      if(disable_avx[X]==0){};
      b[X].gaussian();
    }
    D.apply(b, Db);
    D.dagger(Db, DdaggerDb);
    a[ALL] = 0;
    inverse.apply(DdaggerDb, a);

    diffre = 0;
    onsites(ALL){
      diffre += norm_squared(a[X]-b[X]);
    }
    assert(diffre*diffre < 1e-16 && "test (DdgD)^-1 DdgD");
  }

  // Check conjugate of the wilson Dirac operator
  {
    output0 << "Checking with Dirac_Wilson\n";
    using dirac = Dirac_Wilson<SU<N, double>>;
    dirac D(0.05, U);
    field<Wilson_vector<N, double>> a, b, Db, Ddaggera, DdaggerDb;
    onsites(ALL){
      if(disable_avx[X]==0){};
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
  
    assert(diffre*diffre < 1e-16 && "test dirac_staggered_dagger");
    assert(diffim*diffim < 1e-16 && "test dirac_staggered_dagger");
  
    // Now run CG on Ddaggera. b=1/D a -> Db = a 
    CG<dirac> inverse(D);
    b[ALL] = 0;
    inverse.apply(Ddaggera, b);
    D.apply(b, Db);

    diffre = 0;
    onsites(ALL){
      diffre += norm_squared(a[X]-Db[X]);
    }
    assert(diffre*diffre < 1e-8 && "test CG");
  }

  // Check conjugate of the even-odd preconditioned staggered Dirac operator
  {
    output0 << "Checking with dirac_staggered_evenodd\n";
    dirac_staggered_evenodd D(5.0, U);
    field<SU_vector<N, double>> a, b, Db, Ddaggera, DdaggerDb;
    onsites(ALL){
      if(disable_avx[X]==0){};
      a[X].gaussian();
      b[X].gaussian();
    }

    double diffre = 0, diffim = 0;
    D.apply(b, Db);
    D.dagger(a, Ddaggera);
    onsites(EVEN){
      diffre += a[X].dot(Db[X]).re - Ddaggera[X].dot(b[X]).re;
      diffim += a[X].dot(Db[X]).im - Ddaggera[X].dot(b[X]).im;
    }

    assert(diffre*diffre < 1e-16 && "test dirac_staggered_dagger");
    assert(diffim*diffim < 1e-16 && "test dirac_staggered_dagger");
    
    // Now run CG on Ddaggera. b=1/D a -> Db = a 
    CG inverse(D);
    b[ALL] = 0;
    inverse.apply(Ddaggera, b);
    D.apply(b, Db);

    onsites(EVEN){
      diffre += norm_squared(a[X]-Db[X]);
    }
    assert(diffre*diffre < 1e-8 && "test CG");
  }

  // Check conjugate of the even-odd preconditioned wilson Dirac operator
  {
    output0 << "Checking with Dirac_Wilson_evenodd\n";
    using dirac = Dirac_Wilson_evenodd<SU<N>>;
    dirac D(0.12, U);
    field<Wilson_vector<N, double>> a, b, Db, Ddaggera, DdaggerDb;
    #if NDIM > 3
      a.set_boundary_condition(TUP, boundary_condition_t::ANTIPERIODIC);
      a.set_boundary_condition(TDOWN, boundary_condition_t::ANTIPERIODIC);
      b.copy_boundary_condition(a);
      Db.copy_boundary_condition(a);
      Ddaggera.copy_boundary_condition(a);
      DdaggerDb.copy_boundary_condition(a);
    #endif

    a[ODD] = 0; b[ODD] = 0;
    onsites(EVEN){
      if(disable_avx[X]==0){};
      a[X].gaussian();
      b[X].gaussian();
    }

    double diffre = 0, diffim = 0;
    D.apply(b, Db);
    D.dagger(a, Ddaggera);
    onsites(EVEN){
      diffre += a[X].dot(Db[X]).re - Ddaggera[X].dot(b[X]).re;
      diffim += a[X].dot(Db[X]).im - Ddaggera[X].dot(b[X]).im;
    }
  
    assert(diffre*diffre < 1e-16 && "test Dirac_Wilson_dagger");
    assert(diffim*diffim < 1e-16 && "test Dirac_Wilson_dagger");
  
    // Now run CG on Ddaggera. b=1/D a -> Db = a 
    CG<dirac> inverse(D);
    onsites(EVEN){
      if(disable_avx[X]==0){};
      a[X].gaussian();
    }
    b[ALL] = 0;
    D.dagger(a, Ddaggera);
    inverse.apply(Ddaggera, b);
    D.apply(b, Db);

    diffre = 0;
    onsites(EVEN){
      diffre += norm_squared(a[X]-Db[X]);
    }
    output0 << "D inv Dg " << diffre << "\n";
    assert(diffre*diffre < 1e-16 && "test CG");

    // Run CG on a, multiply by DdaggerD and check the same
    onsites(EVEN){
      if(disable_avx[X]==0){};
      a[X].gaussian();
    }
    b[ALL] = 0;
    inverse.apply(a, b);
    D.apply(b, Db);
    D.dagger(Db, DdaggerDb);

    diffre = 0;
    onsites(EVEN){
      diffre += norm_squared(a[X]-DdaggerDb[X]);
    }
    assert(diffre*diffre < 1e-16 && "test DdgD (DdgD)^-1");

    // The other way around, multiply first
    onsites(EVEN){
      if(disable_avx[X]==0){};
      b[X].gaussian();
    }
    D.apply(b, Db);
    D.dagger(Db, DdaggerDb);
    a[ALL] = 0;
    inverse.apply(DdaggerDb, a);

    diffre = 0;
    onsites(EVEN){
      diffre += norm_squared(a[X]-b[X]);
    }
    assert(diffre*diffre < 1e-16 && "test (DdgD)^-1 DdgD");
  }

  // The the Hasenbusch operator
  {
    output0 << "Checking with Hasenbusch_operator\n";
    using dirac_base = Dirac_Wilson_evenodd<SU<N>>;
    dirac_base Dbase(0.1, U);
    using dirac = Hasenbusch_operator<dirac_base>;
    dirac D(Dbase, 0.1);
    field<Wilson_vector<N, double>> a, b, Db, Ddaggera, DdaggerDb;
    field<Wilson_vector<N, double>> sol;

    a[ODD] = 0; b[ODD] = 0;
    onsites(EVEN){
      if(disable_avx[X]==0){};
      a[X].gaussian();
      b[X].gaussian();
      sol[X] = 0;
    }

    double diffre = 0, diffim = 0;
    D.apply(b, Db);
    D.dagger(a, Ddaggera);
    onsites(EVEN){
      diffre += a[X].dot(Db[X]).re - Ddaggera[X].dot(b[X]).re;
      diffim += a[X].dot(Db[X]).im - Ddaggera[X].dot(b[X]).im;
    }
  
    assert(diffre*diffre < 1e-16 && "test Dirac_Wilson_dagger");
    assert(diffim*diffim < 1e-16 && "test Dirac_Wilson_dagger");
  
    // Now run CG on Ddaggera. b=1/D a -> Db = a 
    CG<dirac> inverse(D);
    b[ALL] = 0;
    inverse.apply(Ddaggera, b);
    D.apply(b, Db);

    onsites(EVEN){
      diffre += norm_squared(a[X]-Db[X]);
    }
    assert(diffre*diffre < 1e-16 && "test CG");

    // Run CG on a, multiply by DdaggerD and check the same
    onsites(EVEN){
      if(disable_avx[X]==0){};
      a[X].gaussian();
    }
    b[ALL] = 0;
    inverse.apply(a, b);
    D.apply(b, Db);
    D.dagger(Db, DdaggerDb);

    diffre = 0;
    onsites(EVEN){
      diffre += norm_squared(a[X]-DdaggerDb[X]);
    }
    assert(diffre*diffre < 1e-16 && "test DdgD (DdgD)^-1");

    // The other way around, multiply first
    onsites(EVEN){
      if(disable_avx[X]==0){};
      b[X].gaussian();
    }
    D.apply(b, Db);
    D.dagger(Db, DdaggerDb);
    a[ALL] = 0;
    inverse.apply(DdaggerDb, a);

    diffre = 0;
    onsites(EVEN){
      diffre += norm_squared(a[X]-b[X]);
    }
    assert(diffre*diffre < 1e-16 && "test (DdgD)^-1 DdgD");
  }


  finishrun();
}


