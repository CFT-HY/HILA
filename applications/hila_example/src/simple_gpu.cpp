#include "hila.h"

static_assert(NDIM == 3, "NDIM must be 3 here");

using MyType = Complex<float>;

int main(int argc, char *argv[]) {

    // initialize system
    hila::initialize(argc, argv);

    // set up the lattice
    lattice->setup({128, 128, 128});

    // Random numbers are used here - use time to seed
    hila::seed_random(0);

    // lattice field
    Field<MyType> f = 1;
    Field<MyType> g;
    // make f Gaussian random distributed
    onsites(ALL) g[X] = f[X] + f[X+e_x];
    
    hila::finishrun();
    return 0;
}