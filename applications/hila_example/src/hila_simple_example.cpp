#include "hila.h"
static_assert(NDIM == 3, "NDIM must be 3");

int main(int argc, char *argv[]) {

    hila::initialize(argc, argv);

    // set up 32^3 lattice
    lattice.setup({32, 32, 32});

    // Random numbers are used here
    hila::seed_random(32345);

    Field<Complex<double>> f;
    Field<double> g = 0;

    // make f Gaussian random distributed
    onsites(ALL) f[X].gaussian_random();

    // calculate sum of 2nd derivatives of f in to g
    foralldir(d) {
        g[ALL] += abs(f[X + d] - 2 * f[X] + f[X - d]);
    }

    // get average of g
    double average = 0;
    onsites(ALL) {
        average += g[X];
    }

    average = average / lattice.volume();
    hila::out0 << "Average of g is " << average << '\n';

    // make a clean exit
    hila::finishrun();
}