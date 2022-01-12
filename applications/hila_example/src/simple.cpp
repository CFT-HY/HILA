#include "hila.h"

static_assert(NDIM == 3, "NDIM must be 3 here");

using MyType = Complex<double>;

int main(int argc, char *argv[]) {

    // initialize system
    hila::initialize(argc, argv);

    // set up 32^3 lattice
    lattice->setup({32, 32, 32});

    // Random numbers are used here
    hila::seed_random(32345);

    // lattice field
    Field<MyType> f;

    // make f Gaussian random distributed
    onsites(ALL) f[X].gaussian_random();

    // Measure hopping term and f^2
    Complex<double> hopping = 0;
    double fsqr = 0;

    onsites(ALL) {
        foralldir(d) {
            hopping += f[X] * f[X + d].conj();
        }
        fsqr += f[X].squarenorm();
    }

    output0 << "Average f^2 : " << fsqr / lattice->volume() << '\n';
    output0 << "Average hopping term " << hopping / (NDIM*lattice->volume()) << '\n';

    hila::finishrun();
    return 0;
}
