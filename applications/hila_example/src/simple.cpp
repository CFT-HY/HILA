#include "hila.h"

static_assert(NDIM == 3, "NDIM must be 3 here");

int main(int argc, char * argv[]) {

    // initialize system
    hila::initialize(argc,argv);

    // set up 32^3 lattice
    lattice->setup({32,32,32});

    // Random numbers are used here
    hila::seed_random(32345);

    // two lattice fields, set g=0
    Field<Complex<double>> f;
    Field<double> g = 0;

    // make f Gaussian random distributed
    onsites(ALL) f[X].gaussian();

    // calculate sum of 2nd derivatives of f to g
    foralldir(d) {
        g[ALL] += abs(f[X+d] - 2*f[X] + f[X-d]);
    }

    // get average value of g
    double ave = 0;
    onsites(ALL) {
        ave += g[X];
    }

    output0 << "Average of g is " << ave/lattice->volume() << '\n';

    // make a clean exit
    hila::finishrun();    
}
