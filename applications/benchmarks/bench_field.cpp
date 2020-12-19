#include "bench.h"
#include "plumbing/coordinates.h"
#include "dirac/conjugate_gradient.h"

#define N 3

#ifndef SEED
#define SEED 100
#endif

int latsize[4] = {32,32,32,32};

int main(int argc, char **argv){
    int n_runs=1;
    double msecs;
    struct timeval start, end;
    double timing;
    double sum;
    float fsum;

    hila::initialize(argc, argv);

    lattice->setup(latsize);

    seed_random(SEED);


    field<double> dfield1, dfield2, dfield3;
    field<float> ffield1, ffield2, ffield3;
    onsites(ALL){
      if(disable_avx[X]==0){};
      dfield1[X] = hila_random();
      dfield2[X] = hila_random();
      dfield3[X] = hila_random();
    }
    onsites(ALL){
      if(disable_avx[X]==0){};
      ffield1[X] = hila_random();
      ffield2[X] = hila_random();
      ffield3[X] = hila_random();
    }

    // Benchmark simple scalar field operation (Memory bandwith)
    timing = 0;
    for(n_runs=1; timing < mintime; ){
      n_runs*=2;
      gettimeofday(&start, NULL);
      for( int i=0; i<n_runs; i++){
          dfield1[ALL] = dfield2[X]*dfield3[X];
      }
      // // synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
      broadcast(timing);

    }
    timing = timing / (double)n_runs;
    output0 << "Double multiply : "<< timing << " ms \n";

    timing = 0;
    for(n_runs=1; timing < mintime; ){
      n_runs*=2;
      gettimeofday(&start, NULL);
      for( int i=0; i<n_runs; i++){
          dfield1[ALL] = dfield2[X] + dfield3[X];
      }
      // // synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
      broadcast(timing);

    }
    timing = timing / (double)n_runs;
    output0 << "Double add : "<< timing << " ms \n";

    timing = 0;
    for(n_runs=1; timing < mintime; ){
      n_runs*=2;
      gettimeofday(&start, NULL);
      for( int i=0; i<n_runs; i++){
          ffield1[ALL] = ffield2[X]*ffield3[X];
      }
      // synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
      broadcast(timing);

    }
    timing = timing / (double)n_runs;
    output0 << "Float multiply : "<< timing << " ms \n";

    timing = 0;
    for(n_runs=1; timing < mintime; ){
      n_runs*=2;
      gettimeofday(&start, NULL);
      for( int i=0; i<n_runs; i++){
          ffield1[ALL] = ffield2[X] + ffield3[X];
      }
      // synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
      broadcast(timing);

    }
    timing = timing / (double)n_runs;
    output0 << "Float add : "<< timing << " ms \n";


    field<SquareMatrix<N, cmplx<double>> > matrix1;
    field<SquareMatrix<N, cmplx<double>> > matrix2;
    field<SquareMatrix<N, cmplx<double>> > matrix3;
    field<vector<N, cmplx<double>> > vector1;
    field<vector<N, cmplx<double>> > vector2;
    field<SquareMatrix<N, cmplx<float>> > fmatrix1;
    field<SquareMatrix<N, cmplx<float>> > fmatrix2;
    field<SquareMatrix<N, cmplx<float>> > fmatrix3;
    field<vector<N, cmplx<float>> > fvector1;
    field<vector<N, cmplx<float>> > fvector2;

    // Generate random values
    onsites(ALL){
      if(disable_avx[X]==0){};
      matrix1[X].random();
      matrix2[X].random();
      matrix3[X].random();
      vector1[X].random();
      vector2[X].random();
    }

    onsites(ALL){
      if(disable_avx[X]==0){};
      fmatrix1[X].random();
      fmatrix2[X].random();
      fmatrix3[X].random();
      fvector1[X].random();
      fvector2[X].random();
    }

    // Interesting case of using the same memory three times
    timing = 0;
    for(n_runs=1; timing < mintime; ){
      n_runs*=2;
      gettimeofday(&start, NULL);
      for( int i=0; i<n_runs; i++){
          matrix1[ALL] = matrix1[X]*matrix1[X];
      }
      // synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
      broadcast(timing);

    }
    timing = timing / (double)n_runs;
    output0 << "Matrix1 = Matrix1 * Matrix1 : "<< timing << " ms \n";


    // Time MATRIX * MATRIX
    timing = 0;
    for(n_runs=1; timing < mintime;){
      n_runs*=2;
      gettimeofday(&start, NULL);
      for( int i=0; i<n_runs; i++){
        matrix3[ALL] = matrix1[X]*matrix2[X];
      }
      // synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
      broadcast(timing);

    }
    timing = timing / (double)n_runs;
    output0 << "Matrix * Matrix: " << timing << "ms \n";


    // Time MATRIX * MATRIX
    timing = 0;
    for(n_runs=1; timing < mintime;){
      n_runs*=2;
      gettimeofday(&start, NULL);
      for( int i=0; i<n_runs; i++){
          fmatrix3[ALL] = fmatrix1[X]*fmatrix2[X];
      }
      // synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
      broadcast(timing);

    }
    timing = timing / (double)n_runs;
    output0 << "Single Precision Matrix * Matrix: " << timing << "ms \n";

    // Time VECTOR * MATRIX
    timing = 0;
    for(n_runs=1; timing < mintime; ){
      n_runs*=2;
      gettimeofday(&start, NULL);
      for( int i=0; i<n_runs; i++){
        vector2[ALL] = vector1[X]*matrix1[X];
      }
      // synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
      broadcast(timing);

    }
    timing = timing / (double)n_runs;
    output0 << "Vector * Matrix: " << timing << " ms \n";

    // Time VECTOR * MATRIX
    timing = 0;
    for(n_runs=1; timing < mintime;){
      n_runs*=2;
      gettimeofday(&start, NULL);
      for( int i=0; i<n_runs; i++){
          fvector2[ALL] = fvector1[X]*fmatrix1[X];
      }
      // synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
      broadcast(timing);
      //output0 << "timing " << timing << '\n';
    }
    timing = timing / (double)n_runs;
    output0 << "Single Precision Vector * Matrix: " << timing << " ms \n";



    // Time VECTOR NORM
    timing = sum = 0;
    onsites(ALL){ // Warm up. Why does this affect the time?
      sum += vector1[X].norm_sq();
    }
    for(n_runs=1; timing < mintime;){
      n_runs*=2;
      gettimeofday(&start, NULL);
      
      sum=0;
      for( int i=0; i<n_runs; i++){
        onsites(ALL){
          sum += vector1[X].norm_sq();
        }
      }
      // synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
      broadcast(timing);
    }
    timing = timing / (double)n_runs;
    output0 << "Vector square sum: " << timing << " ms \n";

    // Time FLOAT VECTOR NORM
    timing = 0;
    for(n_runs=1; timing < mintime;){
      n_runs*=2;
      gettimeofday(&start, NULL);
      fsum=0;
      for( int i=0; i<n_runs; i++){
        onsites(ALL){
          fsum += fvector1[X].norm_sq();
        }
      }
      // synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
      broadcast(timing);
    }
    timing = timing / (double)n_runs;
    output0 << "Single Precision vector square sum: " << timing << " ms \n";

    // Time COMMUNICATION of a MATRIX
    timing = 0;
    for(n_runs=1; timing < mintime; ){
      n_runs*=2;
      gettimeofday(&start, NULL);

      for( int i=0; i<n_runs; i++){
        matrix1.mark_changed(ALL);
        for(int dir=0; dir<NDIRS; dir++){
          matrix1.fetch((direction)dir,ALL);
        }
      }
      
      // synchronize();
      gettimeofday(&end, NULL);
      timing = timediff(start, end);
      broadcast(timing);

    }
    timing = timing / 2 / NDIRS / (double)n_runs;
    output0 << "Matrix nearest neighbour communication: " << timing << " ms \n";

    hila::finishrun();
}



