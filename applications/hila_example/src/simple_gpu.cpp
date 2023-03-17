#include "hila.h"

static_assert(NDIM == 3, "NDIM must be 3 here");

using MyType = double;

int main(int argc, char *argv[]) {

    // initialize system
    hila::initialize(argc, argv);

    // set up the lattice
    lattice.setup({128, 128, 128});

    // Random numbers are used here - use time to seed
    hila::seed_random(0);
    double reduce = 0;
    // lattice field
    Field<MyType> f = 1;
    // make f Gaussian random distributed
    onsites(ALL) reduce += f[X];
    hila::out0 << reduce/lattice.volume()<< '\n';

    f.generate_random_field();
    onsites(ALL) reduce += f[X];
    hila::out0 << reduce/lattice.volume() << '\n';

    hila::finishrun();
    return 0;
}