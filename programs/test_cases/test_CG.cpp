#include "test.h"
#include "../plumbing/algorithms/conjugate_gradients.h"

#define N 3

int main(int argc, char **argv){

    test_setup(argc, argv);

    field<matrix<1, N, cmplx<double>> > b;
    field<matrix<1, N, cmplx<double>> > sol;
    field<matrix<N,N, cmplx<double>> > U[NDIM];

    b.mark_changed(ALL);
    sol.mark_changed(ALL);

    onsites(ALL){
        b[X].random();
        sol[X].random();
    }

    foralldir(d) {
      U[d].mark_changed(ALL);
      onsites(ALL){
        U[d][X].random();
      }
    }

    CG_engine<staggered_dirac> engine;
    engine.solve(U, 0.1, b, sol);

    finishrun();
}
