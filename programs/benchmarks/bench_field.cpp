#include "bench.h"

#define N 3
constexpr int mintime = CLOCKS_PER_SEC;

#ifndef SEED
#define SEED 100
#endif





int main(int argc, char **argv){
    int n_runs=1;
    double msecs;
    clock_t init, end;
    double timing;
    double sum;
    float fsum;

    // Runs lattice->setup 
    bench_setup(argc, argv);
    seed_random(SEED);


    field<double> dfield1, dfield2, dfield3;
    field<float> ffield1, ffield2, ffield3;
    onsites(ALL){
      dfield1[X] = hila_random();
      dfield2[X] = hila_random();
      dfield3[X] = hila_random();
    }
    onsites(ALL){
      ffield1[X] = hila_random();
      ffield2[X] = hila_random();
      ffield3[X] = hila_random();
    }

    // Benchmark simple scalar field operation (Memory bandwith)
    init = end = 0;
    for(n_runs=1; (end-init) < mintime; ){
      n_runs*=2;
      init = clock();
      for( int i=0; i<n_runs; i++){
          dfield1[ALL] = dfield2[X]*dfield3[X];
      }
      synchronize();
      end = clock();
    }
    timing = (end - init) *1000.0 / ((double) CLOCKS_PER_SEC) / (double)n_runs;
    output0 << "Double multiply : "<< timing << " ms \n";

    init = end = 0;
    for(n_runs=1; (end-init) < mintime; ){
      n_runs*=2;
      init = clock();
      for( int i=0; i<n_runs; i++){
          dfield1[ALL] = dfield2[X] + dfield3[X];
      }
      synchronize();
      end = clock();
    }
    timing = (end - init) *1000.0 / ((double) CLOCKS_PER_SEC) / (double)n_runs;
    output0 << "Double add : "<< timing << " ms \n";

    init = end = 0;
    for(n_runs=1; (end-init) < mintime; ){
      n_runs*=2;
      init = clock();
      for( int i=0; i<n_runs; i++){
          ffield1[ALL] = ffield2[X]*ffield3[X];
      }
      synchronize();
      end = clock();
    }
    timing = (end - init) *1000.0 / ((double) CLOCKS_PER_SEC) / (double)n_runs;
    output0 << "Float multiply : "<< timing << " ms \n";

    init = end = 0;
    for(n_runs=1; (end-init) < mintime; ){
      n_runs*=2;
      init = clock();
      for( int i=0; i<n_runs; i++){
          ffield1[ALL] = ffield2[X] + ffield3[X];
      }
      synchronize();
      end = clock();
    }
    timing = (end - init) *1000.0 / ((double) CLOCKS_PER_SEC) / (double)n_runs;
    output0 << "Float add : "<< timing << " ms \n";


    field<matrix<N,N, cmplx<double>> > matrix1;
    field<matrix<N,N, cmplx<double>> > matrix2;
    field<matrix<N,N, cmplx<double>> > matrix3;
    field<matrix<1,N, cmplx<double>> > vector1;
    field<matrix<1,N, cmplx<double>> > vector2;
    field<matrix<N,N, cmplx<float>> > fmatrix1;
    field<matrix<N,N, cmplx<float>> > fmatrix2;
    field<matrix<N,N, cmplx<float>> > fmatrix3;
    field<matrix<1,N, cmplx<float>> > fvector1;
    field<matrix<1,N, cmplx<float>> > fvector2;


    // NOTE: This is because of the failure of transformer to recognize
    // the member function call as changing the object!
    matrix1.mark_changed(ALL);
    matrix2.mark_changed(ALL);
    matrix3.mark_changed(ALL);
    vector1.mark_changed(ALL);
    vector2.mark_changed(ALL);
    fmatrix1.mark_changed(ALL);
    fmatrix2.mark_changed(ALL);
    fmatrix3.mark_changed(ALL);
    fvector1.mark_changed(ALL);
    fvector2.mark_changed(ALL);

    // Generate random values
    onsites(ALL){
      matrix1[X].random();
      matrix2[X].random();
      matrix3[X].random();
      vector1[X].random();
      vector2[X].random();
    }

    onsites(ALL){
      fmatrix1[X].random();
      fmatrix2[X].random();
      fmatrix3[X].random();
      fvector1[X].random();
      fvector2[X].random();
    }

    // Interesting case of using the same memory three times
    init = end = 0;
    for(n_runs=1; (end-init) < mintime; ){
      n_runs*=2;
      init = clock();
      for( int i=0; i<n_runs; i++){
          matrix1[ALL] = matrix1[X]*matrix1[X];
      }
      synchronize();
      end = clock();
    }
    timing = (end - init) *1000.0 / ((double) CLOCKS_PER_SEC) / (double)n_runs;
    output0 << "Matrix1 = Matrix1 * Matrix1 : "<< timing << " ms \n";


    // Time MATRIX * MATRIX
    init = end = 0;
    for(n_runs=1; (end-init) < mintime;){
      n_runs*=2;
      init = clock();
      for( int i=0; i<n_runs; i++){
        matrix3[ALL] = matrix1[X]*matrix2[X];
      }
      synchronize();
      end = clock();
    }
    timing = (end - init) *1000.0 / ((double) CLOCKS_PER_SEC) / (double)n_runs;
    output0 << "Matrix * Matrix: " << timing << "ms \n";


    // Time MATRIX * MATRIX
    init = end = 0;
    for(n_runs=1; (end-init) < mintime;){
      n_runs*=2;
      init = clock();
      for( int i=0; i<n_runs; i++){
          fmatrix3[ALL] = fmatrix1[X]*fmatrix2[X];
      }
      synchronize();
      end = clock();
    }
    timing = (end - init) *1000.0 / ((double)CLOCKS_PER_SEC) / (double)n_runs;
    output0 << "Single Precision Matrix * Matrix: " << timing << "ms \n";

    // Time VECTOR * MATRIX
    init = end = 0;
    for(n_runs=1; (end-init) < mintime; ){
      n_runs*=2;
      init = clock();
      for( int i=0; i<n_runs; i++){
        vector2[ALL] = vector1[X]*matrix1[X];
      }
      synchronize();
      end = clock();
    }
    timing = (end - init) *1000.0 / ((double)CLOCKS_PER_SEC) / (double)n_runs;
    output0 << "Vector * Matrix: " << timing << " ms \n";

    // Time VECTOR * MATRIX
    init = end = 0;
    for(n_runs=1; (end-init) < mintime;){
      n_runs*=2;
      init = clock();
      for( int i=0; i<n_runs; i++){
          fvector2[ALL] = fvector1[X]*fmatrix1[X];
      }
      synchronize();
      end = clock();
    }
    timing = (end - init) *1000.0 / ((double)CLOCKS_PER_SEC) / (double)n_runs;
    output0 << "Single Precision Vector * Matrix: " << timing << " ms \n";



    // Time VECTOR NORM
    init = end = 0;
    for(n_runs=1; (end-init) < mintime; ){
      n_runs*=2;
      init = clock();
      
      sum=0;
      for( int i=0; i<n_runs; i++){
        onsites(ALL){
          sum += norm_squared(vector1[X]);
        }
      }
      volatile double volatile_sum = sum;
      synchronize();
      end = clock();
    }
    timing = (end - init) *1000.0 / ((double)CLOCKS_PER_SEC) / (double)n_runs;
    output0 << "Vector square sum: " << timing << " ms \n";

    // Time FLOAT VECTOR NORM
    init = end = 0;
    for(n_runs=1; (end-init) < mintime; n_runs*=2){
      init = clock();
      fsum=0;
      for( int i=0; i<n_runs; i++){
        onsites(ALL){
          for(int j=0; j<N; j++){
            fvector1[X].c[0][j]=1;
          }
        }
        onsites(ALL){
          fsum += fvector1[X].norm_sq();
        }
        volatile double volatile_sum = fsum;
      }
      synchronize();
      end = clock();
    }
    timing = (end - init) *1000.0 / ((double)CLOCKS_PER_SEC) / (double)n_runs;
    output0 << "Single Precision vector square sum: " << timing << " ms \n";

    // Time COMMUNICATION of a MATRIX
    init = end = 0;
    for(n_runs=1; (end-init) < mintime; ){
      n_runs*=2;
      init = clock();

      for( int i=0; i<n_runs; i++){
        matrix1.mark_changed(ALL);
        for(int dir=0; dir<NDIRS; dir++){
          matrix1.wait_move((direction)dir,ALL);
        }
      }
      
      synchronize();
      end = clock();
    }
    timing = (end - init) *1000.0 / ((double)CLOCKS_PER_SEC) / 2 / NDIRS / (double)n_runs;
    output0 << "Matrix nearest neighbour communication: " << timing << " ms \n";

    //printf("node %d, create gauge\n", mynode());

    // Define a gauge matrix
    field<matrix<N,N, cmplx<double>> > U[NDIM];
    foralldir(d) {
      U[d].mark_changed(ALL);
      onsites(ALL){
        U[d][X].random();
        vector1[X].random();
        vector2[X].random();
      }
    }


    // Time staggered Dirac operator
    init = end = 0;
    //printf("node %d, dirac_stagggered 0\n", mynode());
    dirac_stagggered(U, 0.1, vector1, vector2);
    synchronize();
    for(n_runs=1; (end-init) < mintime; ){
      n_runs*=2;
      init = clock();
      for( int i=0; i<n_runs; i++){
        //printf("node %d, dirac_stagggered %d\n", mynode(), i);
        vector1.mark_changed(ALL); // Assure communication is included
        dirac_stagggered(U, 0.1, vector1, vector2);
      }
      synchronize();
      end = clock();
    }
    timing = (end - init) *1000.0 / ((double)CLOCKS_PER_SEC) / (double)n_runs;
    output0 << "Dirac: " << timing << "ms \n";


    // Time staggered Dirac operator with direction loop expanded
    #if (NDIM==4) 
    init = end = 0;
    dirac_stagggered_4dim(U, 0.1, vector1, vector2);
    synchronize();
    for(n_runs=1; (end-init) < mintime; ){
      n_runs*=2;
      init = clock();
      for( int i=0; i<n_runs; i++){
        vector1.mark_changed(ALL); // Assure communication is included
        dirac_stagggered_4dim(U, 0.1, vector1, vector2);
      }
      synchronize();
      end = clock();
    }
    timing = (end - init) *1000.0 / ((double)CLOCKS_PER_SEC) / (double)n_runs;
    output0 << "Dirac one loop: " << timing << "ms \n";
    #endif
    

    // Conjugate gradient step 
    init = end = 0;
    field<matrix<1,N, cmplx<double>> > r, rnew, p, Dp;
    for(n_runs=1; (end-init) < mintime; ){
      n_runs*=2;

      double pDDp = 0, rr = 0, rrnew = 0;
      double alpha, beta;

      onsites(ALL){
        r[X] = vector1[X];
        p[X] = vector1[X];
        for(int i=0; i<N; i++){
           vector2[X].c[0][i] = 0;
        }
      }

      init = clock();
      for( int i=0; i<n_runs; i++){
        
        dirac_stagggered(U, 0.1, p, Dp);

        rr=pDDp=0;
        onsites(ALL){
          rr += r[X].norm_sq();
          pDDp += Dp[X].norm_sq();
        }

        alpha = rr / pDDp;

        rrnew = 0;
        onsites(ALL){
          vector2[X] = r[X] + alpha*p[X];
          r[X] = r[X] - alpha*Dp[X];
          rrnew += r[X].norm_sq();
        }

        beta = rrnew/rr;
        p[ALL] = r[X] + beta*p[X];
      }
      synchronize();
      end = clock();
    }

    timing = (end - init) *1000.0 / ((double)CLOCKS_PER_SEC) / (double)n_runs;
    output0 << "CG: " << timing << "ms / iteration\n";


    finishrun();
}



