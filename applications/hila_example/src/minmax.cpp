#include "hila.h"

#include <random>

static_assert(NDIM == 3, "NDIM must be 3 here");


using MyType = Complex<double>;

int main(int argc, char *argv[]) {

    // initialize system
    hila::initialize(argc, argv);

    // set up 32^3 lattice
    lattice->setup({32, 32, 32});

    // Random numbers are used here - use time to seed
    hila::seed_random(0);

    // lattice field
    Field<double> g = 2.0;

    // make f Gaussian random distributed
    g[{2,3,4}] = 1.0;

    double val;
    CoordinateVector loc;
    val = g.min(loc);
    output0 << "Max value " << val << '\n';
    output0 << "Location:" << loc << '\n';

    hila::finishrun();
    return 0;
}
