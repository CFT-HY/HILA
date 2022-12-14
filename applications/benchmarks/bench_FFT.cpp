#include "bench.h"
#include "plumbing/timing.h"
#include "plumbing/fft.h"
//#include "plumbing/FFT.h"

#ifndef SEED
#define SEED 100
#endif

const CoordinateVector latsize = {32, 32, 32, 32};

int main(int argc, char **argv) {
    int n_runs = 1;
    struct timeval start, end;
    double timing;

    hila::initialize(argc, argv);

    lattice.setup(latsize);

    hila::seed_random(SEED);

    using T = Matrix<2, 2, Complex<double>>;
    using Tf = Matrix<2, 2, Complex<float>>;

    Field<T> d, d2;

    // Generate a random field
    onsites(ALL) { d[X].random(); }

    // Run once to make sure everything is set up
    FFT_field(d, d2);

    timing = 0;
    for (n_runs = 1; timing < mintime;) {
        n_runs *= 2;
        gettimeofday(&start, NULL);
        for (int i = 0; i < n_runs; i++) {
            FFT_field(d, d2);
        }
        // // synchronize();
        gettimeofday(&end, NULL);
        timing = timediff(start, end);
        hila::broadcast(timing);
    }
    timing = timing / (double)n_runs;
    hila::out0 << "FFT double precision : " << timing << " ms \n";

    // Generate a random field
    Field<Tf> f, f2;

    // Generate a random field
    onsites(ALL) { f[X].random(); }

    timing = 0;
    for (n_runs = 1; timing < mintime;) {
        n_runs *= 2;
        gettimeofday(&start, NULL);
        for (int i = 0; i < n_runs; i++) {
            FFT_field(f, f2);
        }
        // // synchronize();
        gettimeofday(&end, NULL);
        timing = timediff(start, end);
        hila::broadcast(timing);
    }
    timing = timing / (double)n_runs;
    hila::out0 << "FFT single precision : " << timing << " ms \n";

    hila::finishrun();
}
