#include "bench.h"

#define N 3
constexpr int mintime = CLOCKS_PER_SEC;

int main(int argc, char **argv){
    int n_runs=1;
    double msecs;
    clock_t init, end;
    double timing;
    double sum;


    // Runs lattice->setup 
    bench_setup(argc, argv);

    field<matrix<N,N, cmplx<double>> > matrix1;
    field<matrix<N,N, cmplx<double>> > matrix2;
    field<matrix<N,N, cmplx<double>> > matrix3;
    field<matrix<1,N, cmplx<double>> > vector1;
    field<matrix<1,N, cmplx<double>> > vector2;

    matrix1[ALL] = 1;
    matrix2[ALL] = 1;

    // Time MATRIX * MATRIX
    init = end = 0;
    for(n_runs=1; (end-init) < mintime; n_runs*=2){
      init = clock();
      for( int i=0; i<n_runs; i++){
          matrix3[ALL] = matrix1[X]*matrix2[X];
      }
      synchronize();
      end = clock();
    }
    timing = (end - init) *1000.0 / ((double) CLOCKS_PER_SEC) / (double)n_runs;
    output0 << "Matrix * Matrix: " << timing << "ms \n";


    matrix1[ALL] = 1; 
    onsites(ALL){
        for(int i=0; i<N; i++){
            vector1[X].c[0][1]=1;
        }
    }

    // Time VECTOR * MATRIX
    init = end = 0;
    for(n_runs=1; (end-init) < mintime; n_runs*=2){
      init = clock();
      for( int i=0; i<n_runs; i++){
          vector1[ALL] = vector1[X]*matrix1[X];
      }
      synchronize();
      end = clock();
    }
    timing = (end - init) *1000.0 / ((double)CLOCKS_PER_SEC) / (double)n_runs;
    output0 << "Vector * Matrix: " << timing << " ms \n";


    // Time VECTOR * MATRIX
    init = end = 0;
    for(n_runs=1; (end-init) < mintime; n_runs*=2){
      init = clock();
      
      sum=0;
      for( int i=0; i<n_runs; i++){
        onsites(ALL){
          sum += norm_sq(vector1[X]);
        }
      }
      synchronize();
      end = clock();
    }
    timing = (end - init) *1000.0 / ((double)CLOCKS_PER_SEC) / (double)n_runs;
    output0 << "Vector square sum: " << timing << " ms \n";

    //printf("node %d, create gauge\n", mynode());

    // Define a gauge matrix
    field<matrix<N,N, cmplx<double>> > U[NDIM];
    foralldir(d) U[d][ALL] = 1;


    // Time naive Dirac operator
    init = end = 0;
    //printf("node %d, dirac_stagggered 0\n", mynode());
    dirac_stagggered(U, 0.1, vector1, vector2);
    synchronize();
    for(n_runs=1; (end-init) < mintime; n_runs*=2){
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


    // Time naive Dirac operator with direction loop expanded
    #if (NDIM==4) 
    init = end = 0;
    dirac_stagggered_alldim(U, 0.1, vector1, vector2);
    synchronize();
    for(n_runs=1; (end-init) < mintime; n_runs*=2){
      init = clock();
      for( int i=0; i<n_runs; i++){
        vector1.mark_changed(ALL); // Assure communication is included
        dirac_stagggered_alldim(U, 0.1, vector1, vector2);
      }
      synchronize();
      end = clock();
    }
    timing = (end - init) *1000.0 / ((double)CLOCKS_PER_SEC) / (double)n_runs;
    output0 << "Dirac separate loops: " << timing << "ms \n";
    #endif
    

    // Conjugate gradient step 
    init = end = 0;
    field<matrix<1,N, cmplx<double>> > r, rnew, p, Dp;
    for(n_runs=1; (end-init) < mintime; n_runs*=2){
      init = clock();

      for( int i=0; i<n_runs; i++){
        double pDDp = 0, rr = 0, rrnew = 0;
        double alpha, beta;

        onsites(ALL){
          r[X] = vector1[X];
          p[X] = vector1[X];
          for(int i=0; i<N; i++){
             vector2[X].c[0][i] = 0;
          }
        }
            
        dirac_stagggered(U, 0.1, p, Dp);

        onsites(ALL){
            rr += norm_sq(r[X]);
            pDDp += norm_sq(Dp[X]);
        }

        alpha = rr / pDDp;

        onsites(ALL){
          vector2[X] = r[X] + alpha*p[X];
          r[X] = r[X] - alpha*Dp[X];
          rrnew += norm_sq(r[X]);
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



