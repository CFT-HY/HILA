#include "bench.h"
#include "plumbing/coordinates.h"
#include "dirac/conjugate_gradient.h"

#define N 3

#ifndef SEED
#define SEED 100
#endif

const CoordinateVector latsize = {32, 32, 32, 32};

int main(int argc, char **argv) {
    int n_runs = 1;
    double msecs;
    struct timeval start, end;
    double timing;
    double sum;
    float fsum;

    hila::initialize(argc, argv);

    lattice.setup(latsize);

    hila::seed_random(SEED);

    // Define a gauge matrix
    Field<SU<N, double>> U[NDIM];
    Field<SU_vector<N, double>> sunvec1, sunvec2;

    foralldir(d) {
        onsites(ALL) { U[d][X].random(); }
    }
    onsites(ALL) {
        sunvec1[X].gaussian_random();
        sunvec2[X].gaussian_random();
    }

    // Time staggered Dirac operator
    timing = 0;

    using sunvec = SU_vector<N, double>;
    using sunmat = SU<N, double>;
    using dirac_stg = dirac_staggered<sunmat>;
    dirac_stg D_staggered(0.1, U);
    D_staggered.apply(sunvec1, sunvec2);

    // synchronize();
    for (n_runs = 1; timing < mintime;) {
        n_runs *= 2;
        gettimeofday(&start, NULL);
        for (int i = 0; i < n_runs; i++) {
            sunvec1.mark_changed(ALL); // Ensure communication is included
            D_staggered.apply(sunvec1, sunvec2);
        }
        // synchronize();
        gettimeofday(&end, NULL);
        timing = timediff(start, end);
        hila::broadcast(timing);
    }
    timing = timing / (double)n_runs;
    hila::out0 << "Dirac staggered: " << timing << "ms \n";

    // Conjugate gradient step
    CG<dirac_stg> stg_inverse(D_staggered, 1e-5, 1);
    timing = 0;
    sunvec1[ALL] = 0;
    stg_inverse.apply(sunvec2, sunvec1);
    for (n_runs = 1; timing < mintime;) {
        n_runs *= 2;

        gettimeofday(&start, NULL);
        for (int i = 0; i < n_runs; i++) {
            sunvec1[ALL] = 0;
            stg_inverse.apply(sunvec2, sunvec1);
        }

        // synchronize();
        gettimeofday(&end, NULL);
        timing = timediff(start, end);
        hila::broadcast(timing);
    }

    timing = timing / (double)n_runs;
    hila::out0 << "Staggered CG: " << timing << "ms / iteration\n";

    Field<Wilson_vector<N, double>> wvec1, wvec2;
    onsites(ALL) {
        wvec1[X].gaussian_random();
        wvec2[X].gaussian_random();
    }
    // Time staggered Dirac operator
    timing = 0;
    // printf("node %d, dirac_staggered 0\n", hila::myrank());
    using Dirac_Wilson = Dirac_Wilson_evenodd<sunmat>;
    Dirac_Wilson D_wilson(0.05, U);
    D_wilson.apply(wvec1, wvec2);

    // synchronize();
    for (n_runs = 1; timing < mintime;) {
        n_runs *= 2;
        gettimeofday(&start, NULL);
        for (int i = 0; i < n_runs; i++) {
            wvec1.mark_changed(ALL); // Ensure communication is included
            D_wilson.apply(wvec1, wvec2);
        }
        // synchronize();
        gettimeofday(&end, NULL);
        timing = timediff(start, end);
        hila::broadcast(timing);
    }
    timing = timing / (double)n_runs;
    hila::out0 << "Dirac Wilson: " << timing << "ms \n";

    // Conjugate gradient step (set accuracy=1 to run only 1 step)
    CG<Dirac_Wilson> w_inverse(D_wilson, 1e-12, 5);
    timing = 0;
    wvec1[ALL] = 0;
    w_inverse.apply(wvec2, wvec1);
    for (n_runs = 1; timing < mintime;) {

        n_runs *= 2;

        gettimeofday(&start, NULL);
        for (int i = 0; i < n_runs; i++) {
            wvec1[ALL] = 0;
            w_inverse.apply(wvec2, wvec1);
        }

        // synchronize();
        gettimeofday(&end, NULL);
        timing = timediff(start, end);
        hila::broadcast(timing);
    }

    timing = timing / (double)n_runs;
    hila::out0 << "Dirac Wilson CG: " << timing << "ms / iteration\n";

    hila::finishrun();
}
