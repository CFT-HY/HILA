#include "bench.h"

#ifndef MSIZE
#define MSIZE 3
#endif

#ifndef SEED
#define SEED 100
#endif 

constexpr int mintime = CLOCKS_PER_SEC;

int main(int argc, char **argv){
    int n_runs=1;
    double msecs;
    clock_t init, end;
    double timing;
    double sum;

    // Runs lattice->setup 
    bench_setup(argc, argv);
    // Seed rng generator
    seed_random(SEED);

    field<matrix<MSIZE,MSIZE, cmplx<double>> > matrix1;
    field<matrix<MSIZE,MSIZE, cmplx<double>> > matrix2;
    field<matrix<MSIZE,MSIZE, cmplx<double>> > matrix3;
    field<matrix<MSIZE,MSIZE, cmplx<double>> > matrix4;

    onsites(ALL){
      matrix<MSIZE,MSIZE, cmplx<double> > rand1, rand2, rand3;
      rand1.random();
      rand2.random();
      rand3.random();
      matrix1[X] = rand1;
      matrix2[X] = rand2;
      matrix3[X] = rand3;
    }

    // Time conj(matrix) * matrix * conj(matrix) 
    output0 << "matrix size: " << (int) MSIZE << "*"  << (int) MSIZE << "\n";
    init = end = 0;
    for(n_runs=1; (end-init) < mintime; n_runs*=2){
      init = clock();
      for( int i=0; i<n_runs; i++){
          matrix4[ALL] = conj(matrix1[X])*matrix2[X]*conj(matrix3[X]);
      }
      synchronize();
      end = clock();
    }
    timing = (end - init) *1000.0 / ((double) CLOCKS_PER_SEC) / (double)n_runs;
    output0 << "conj(matrix) * matrix * conj(matrix): " << timing << "ms \n";

    // Time matrix.conjugate() * matrix * matrix.conjugate() 
    init = end = 0;
    for(n_runs=1; (end-init) < mintime; n_runs*=2){
      init = clock();
      for( int i=0; i<n_runs; i++){
          matrix4[ALL] = (matrix1[X].conjugate())*matrix2[X]*(matrix3[X].conjugate());
      }
      synchronize();
      end = clock();
    }

    timing = (end - init) *1000.0 / ((double) CLOCKS_PER_SEC) / (double)n_runs;
    output0 << "matrix.conjugate() * matrix * matrix.conjugate(): " << timing << "ms \n";

    finishrun();
}



