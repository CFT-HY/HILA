#include "test.h"
#include "../plumbing/algorithms/conjugate_gradients.h"

#define N 3

int main(int argc, char **argv){

    test_setup(argc, argv);

    field<matrix<1, N, cmplx<double>> > b;
    field<matrix<1, N, cmplx<double>> > sol;
    field<matrix<N,N, cmplx<double>> > U[NDIM];

    onsites(ALL){
      b[X].random();
      for(int i=0; i<N; i++)
        sol[X].c[1][i] = 0;
    }

    foralldir(d) {
      onsites(ALL){
        U[d][X].random();
      }
    }

    CG_engine<staggered_dirac> engine;
    engine.solve(U, 1.0, b, sol);

    finishrun();
}
