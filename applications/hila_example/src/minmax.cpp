#include "hila.h"

#include <random>

static_assert(NDIM == 3, "NDIM must be 3 here");


using MyType = Complex<double>;

int main(int argc, char *argv[]) {

    // initialize system
    hila::initialize(argc, argv);

    // set up 32^3 lattice
    lattice->setup({256, 256, 256});

    // Random numbers are used here - use time to seed
    hila::seed_random(0);

    // lattice field
    Field<double> g = 2.0;

    // make f Gaussian random distributed
    g[{1, 1, 1}] = 2.5;
    g[{1, 0, 1}] = 1.0;

    double val1;
    double val2;
    CoordinateVector loc1;
    CoordinateVector loc2;
    for (auto i = 0; i < 1000; i++)
    {
        val1 = g.min(loc1);
        val2 = g.max(loc2);
    }
    
    // val = g.max(loc);
    output0 << "Min value " << val1 << " at location: " << loc1 << '\n';
    output0 << "Max value " << val2 << " at location: " << loc2 << '\n';

    hila::finishrun();
    return 0;
}
