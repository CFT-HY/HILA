#include "hila.h"
#include <random>

static_assert(NDIM == 3, "NDIM must be 3 here");

using MyType = SquareMatrix<4,Complex<float>>;

int main(int argc, char *argv[]) {

    // initialize system
    hila::initialize(argc, argv);

    // set up the lattice
    lattice->setup({128, 128, 128});

    // Random numbers are used here - use time to seed
    hila::seed_random(0);

    // lattice field
    Field<MyType> f;
    Field<float> g;
    g = 0;
    g.product();
    // make f Gaussian random distributed
    onsites(ALL) f[X].gaussian_random();
    // Measure hopping term and f^2
    MyType hopping = 0;
    hopping = 0;
    double fsqr = 0;

    onsites(ALL) {
        foralldir(d) {
            hopping += f[X] * f[X + d].dagger();
        }
        fsqr += f[X].squarenorm();
    }
    output0 << "Average f^2 : " << fsqr / lattice->volume() << '\n';
    output0 << "Average hopping term " << hopping / (NDIM*lattice->volume()) << '\n';

    // auto sd = spectraldensity(f,64);
    // for (int i=0; i<sd.size(); i++) output0 << i << ' ' << sd[i] << '\n';
    
    hila::finishrun();
    return 0;
}
