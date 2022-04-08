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
    Field<MyType> f;
    Field<double> g = 0;

    // make f Gaussian random distributed
    onsites(ALL) f[X].gaussian_random();

    foralldir(d) {
        g[ALL] += abs(f[X+d] - 2*f[X] + f[X-d]);
    }

    double val;
    CoordinateVector loc;
    val = g.max(loc);
    output0 << "Max value " << val << '\n';
    output0 << "Location:" << loc << '\n';

    hila::finishrun();
    return 0;
}
