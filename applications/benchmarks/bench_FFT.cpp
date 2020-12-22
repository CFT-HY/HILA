#include "bench.h"
#include "plumbing/timers.h"

#ifndef SEED
#define SEED 100
#endif

int main(int argc, char **argv){
  int n_runs=1;
  struct timeval start, end;
  double timing;

  // Runs lattice->setup 
  bench_setup(argc, argv);
  seed_random(SEED);

  using T = Matrix<2,2,cmplx<double>>;
  using Tf = Matrix<2,2,cmplx<float>>;

  field<T> d, d2;

  // Generate a random field
  onsites(ALL){
    d[X].random();
  }

  // Run once to make sure everything is set up
  FFT_field(d, d2);
  
  timing = 0;
  for(n_runs=1; timing < mintime; ){
    n_runs*=2;
    gettimeofday(&start, NULL);
    for( int i=0; i<n_runs; i++){
      FFT_field(d, d2);
    }
    // // synchronize();
    gettimeofday(&end, NULL);
    timing = timediff(start, end);
    broadcast(timing);
  }
  timing = timing / (double)n_runs;
  output0 << "FFT double precision : "<< timing << " ms \n";

  // Generate a random field
  field<Tf> f, f2;

  // Generate a random field
  onsites(ALL){
    f[X].random();
  }

  timing = 0;
  for(n_runs=1; timing < mintime; ){
    n_runs*=2;
    gettimeofday(&start, NULL);
    for( int i=0; i<n_runs; i++){
      FFT_field(f, f2);
    }
    // // synchronize();
    gettimeofday(&end, NULL);
    timing = timediff(start, end);
    broadcast(timing);
  }
  timing = timing / (double)n_runs;
  output0 << "FFT single precision : "<< timing << " ms \n";


  hila::finishrun();
}
