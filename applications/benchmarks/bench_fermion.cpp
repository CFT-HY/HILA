#include "bench.h"
#include "plumbing/coordinates.h"
#include "dirac/conjugate_gradient.h"

#define N 3

#ifndef SEED
#define SEED 100
#endif

const int latsize[4] = { 8, 8, 8, 8 };

int main(int argc, char **argv){
    int n_runs=1;
    double msecs;
    struct timeval start, end;
    double timing;
    double sum;
    float fsum;

    // Runs lattice->setup 
    #if NDIM==2
    printf("Not running fermion benchmarks for NDIM = 1\n");
    return 0;
    #elif NDIM==2
    lattice->setup( latsize[0], latsize[1], argc, argv );
    #elif NDIM==3
    lattice->setup( latsize[0], latsize[1], latsize[2], argc, argv );
    #elif NDIM==4
    lattice->setup( latsize[0], latsize[1], latsize[2], latsize[3], argc, argv );
    #endif
    seed_random(SEED);


    // Define a gauge matrix
    field<SU<N,double>> U[NDIM];
    field<SU_vector<N, double>> sunvec1, sunvec2;

    foralldir(d) {
      onsites(ALL){
        U[d][X].random();
      }
    }
    onsites(ALL){
      if(disable_avx[X]==0){};
      sunvec1[X].gaussian();
      sunvec2[X].gaussian();
    }


    // Time staggered Dirac operator
    timing = 0;

    using sunvec = SU_vector<N, double>;
    using sunmat = SU<N, double>;
    using dirac_stg = dirac_staggered<sunmat>;
    dirac_stg D_staggered(0.1, U);
    D_staggered.apply(sunvec1, sunvec2);

    // synchronize();
    for(n_runs=1; timing < mintime; ){
      n_runs*=2;
      gettimeofday(&start, NULL);
      for( int i=0; i<n_runs; i++){
        sunvec1.mark_changed(ALL); // Ensure communication is included
        D_staggered.apply(sunvec1, sunvec2);
      }
      // synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
      broadcast(timing);
    }
    timing = timing / (double)n_runs;
    output0 << "Dirac staggered: " << timing << "ms \n";

    // Conjugate gradient step 
    CG<dirac_stg> stg_inverse(D_staggered, 1.0);
    timing = 0;
    for(n_runs=1; timing < mintime; ){
      n_runs*=2;
      
      for( int i=0; i<n_runs; i++){
        sunvec1[ALL]=0;
        stg_inverse.apply(sunvec2, sunvec1);
      }

      // synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
      broadcast(timing);

    }

    timing = timing / (double)n_runs;
    output0 << "Staggered CG: " << timing << "ms / iteration\n";



    field<Wilson_vector<N, double>> wvec1, wvec2;
    onsites(ALL){
      if(disable_avx[X]==0){};
      wvec1[X].gaussian();
      wvec2[X].gaussian();
    }
    // Time staggered Dirac operator
    timing = 0;
    //printf("node %d, dirac_staggered 0\n", mynode());
    using Dirac_Wilson = Dirac_Wilson_evenodd<sunmat>;
    Dirac_Wilson D_wilson(0.05, U);
    D_wilson.apply(wvec1, wvec2);

    // synchronize();
    for(n_runs=1; timing < mintime; ){
      n_runs*=2;
      gettimeofday(&start, NULL);
      for( int i=0; i<n_runs; i++){
        wvec1.mark_changed(ALL); // Ensure communication is included
        D_wilson.apply(wvec1, wvec2);
      }
      // synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
      broadcast(timing);
    }
    timing = timing / (double)n_runs;
    output0 << "Dirac Wilson: " << timing << "ms \n";

    // Conjugate gradient step (set accuracy=1 to run only 1 step)
    CG<Dirac_Wilson> w_inverse(D_wilson, 1.0);
    timing = 0;
    for(n_runs=1; timing < mintime; ){

      n_runs*=2;

      for( int i=0; i<n_runs; i++){
        wvec1[ALL]=0;
        w_inverse.apply(wvec2, wvec1);
      }

      // synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
      broadcast(timing);

    }

    timing = timing / (double)n_runs;
    output0 << "Dirac Wilson CG: " << timing << "ms / iteration\n";

    finishrun();
}



