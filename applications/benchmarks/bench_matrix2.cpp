

#include "bench.h"
#include "plumbing/timers.h"

#ifndef MSIZE
#define MSIZE 3
#endif

#ifndef SEED
#define SEED 100
#endif

#define MADD(x) (MSIZE + x)

///////////////////////////////////////
// benchmark conjugate operations for 
// increasing matrix sizes starting 
// from msize.  
///////////////////////////////////////

int main(int argc, char **argv){
    int n_runs=1;
    double msecs;
    struct timeval start, end;
    double timing;
    double sum;

    // Runs lattice->setup 
    bench_setup(argc, argv);
    // Seed rng generator
    seed_random(SEED);

    timer timer1;
    
    field<matrix<MADD(0),MADD(0), cmplx<double>> > matrix1;
    field<matrix<MADD(1),MADD(1), cmplx<double>> > matrix2;
    field<matrix<MADD(3),MADD(3), cmplx<double>> > matrix3;
    field<matrix<MADD(6),MADD(6), cmplx<double>> > matrix4;

    onsites(ALL){
      matrix1[X].random();
      matrix2[X].random();
      matrix3[X].random();
      matrix4[X].random();
    }

    // Time conj(matrix) * matrix * conj(matrix) 
    timing = 0;
    for(n_runs=1; timing < mintime; ){
      n_runs*=2;
      gettimeofday(&start, NULL);
      for( int i=0; i<n_runs; i++){
          matrix1[ALL] = conj(matrix1[X])*matrix1[X]*conj(matrix1[X]);
      }
      synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
    }
    timing = timing / (double)n_runs;
    output0 << "matrix size " << (int) MSIZE << "*"  << (int) MSIZE << " : "<< timing << " ms \n";

    timer1.start();
    
    // Time conj(matrix) * matrix * conj(matrix) 
    timing = 0;
    for(n_runs=1; timing < mintime; ){
      n_runs*=2;
      gettimeofday(&start, NULL);
      for( int i=0; i<n_runs; i++){
          matrix2[ALL] = conj(matrix2[X])*matrix2[X]*conj(matrix2[X]);
      }
      synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
    }
    timing = timing / (double)n_runs;
    output0 << "matrix size " << (int) MSIZE + 1 << "*"  << (int) MSIZE + 1 << " : "<< timing << " ms \n";

    timer1.end();
    timer1.report("Timer 1");
    
    // Time conj(matrix) * matrix * conj(matrix) 
    timing = 0;
    for(n_runs=1; timing < mintime; ){
      n_runs*=2;
      gettimeofday(&start, NULL);
      for( int i=0; i<n_runs; i++){
          matrix3[ALL] = conj(matrix3[X])*matrix3[X]*conj(matrix3[X]);
      }
      synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
    }
    timing = timing / (double)n_runs;
    output0 << "matrix size " << (int) MSIZE + 3  << "*"  << (int) MSIZE + 3  << " : "<< timing << " ms \n";

    timing = 0;
    for(n_runs=1; timing < mintime; ){
      n_runs*=2;
      gettimeofday(&start, NULL);
      for( int i=0; i<n_runs; i++){
          matrix4[ALL] = conj(matrix4[X])*matrix4[X]*conj(matrix4[X]);
      }
      synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
    }
    timing = timing / (double)n_runs;
    output0 << "matrix size " << (int) MSIZE + 6 << "*"  << (int) MSIZE + 6 << " : "<< timing << " ms \n";


    //------------------------------------------------

    // Time conj(matrix) * matrix * conj(matrix) 
    timing = 0;
    for(n_runs=1; timing < mintime; ){
      n_runs*=2;
      gettimeofday(&start, NULL);
      for( int i=0; i<n_runs; i++){
          matrix1[ALL] = matrix1[X].conjugate()*matrix1[X]*matrix1[X].conjugate();
      }
      synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
    }
    timing = timing / (double)n_runs;
    output0 << "matrix size " << (int) MSIZE << "*"  << (int) MSIZE << " : "<< timing << " ms \n";

    // Time matrix) * .conjugate()atrix * matrix) 
    timing = 0;
    for(n_runs=1; timing < mintime; ){
      n_runs*=2;
      gettimeofday(&start, NULL);
      for( int i=0; i<n_runs; i++){
          matrix2[ALL] = matrix2[X].conjugate()*matrix2[X]*matrix2[X].conjugate();
      }
      synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
    }
    timing = timing / (double)n_runs;
    output0 << "matrix size " << (int) MSIZE + 1 << "*"  << (int) MSIZE + 1 << " : "<< timing << " ms \n";

    // Time matrix) * .conjugate()atrix * matrix) 
    timing = 0;
    for(n_runs=1; timing < mintime; ){
      n_runs*=2;
      gettimeofday(&start, NULL);
      for( int i=0; i<n_runs; i++){
          matrix3[ALL] = matrix3[X].conjugate()*matrix3[X]*matrix3[X].conjugate();
      }
      synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
    }
    timing = timing / (double)n_runs;
    output0 << "matrix size " << (int) MSIZE + 3  << "*"  << (int) MSIZE + 3  << " : "<< timing << " ms \n";

    timing = 0;
    for(n_runs=1; timing < mintime; ){
      n_runs*=2;
      gettimeofday(&start, NULL);
      for( int i=0; i<n_runs; i++){
          matrix4[ALL] = matrix4[X].conjugate()*matrix4[X]*matrix4[X].conjugate();
      }
      synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
    }
    timing = timing / (double)n_runs;
    output0 << "matrix size " << (int) MSIZE + 6 << "*"  << (int) MSIZE + 6 << " : "<< timing << " ms \n";

    finishrun();
}



