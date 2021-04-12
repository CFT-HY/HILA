#include "plumbing/defs.h"
#include "plumbing/field.h"
#include "datatypes/cmplx.h"
#include "datatypes/matrix.h"
#include "plumbing/fft.h"
#include "plumbing/input.h"

using vtype = Complex<double>;

static_assert(NDIM == 3, "NDIM must be 3 in this program");

int main(int argc, char **argv) {

    hila::initialize(argc, argv);

    hila::input par("parameters");
    int nx = par.get("nx");
    int ny = par.get("ny");
    int nz = par.get("nz");
    int loops = par.get("loops");
    int seed = par.get("random seed");

    par.close();

    lattice->setup({nx, ny, nz});

    hila::seed_random(seed);

    Field<vtype> f, g;

    f[ALL] = hila::gaussian_ran();

    static hila::timer cmplx_timer("cmplx fft");

    for (int i = 0; i < loops; i++) {
        cmplx_timer.start();
        FFT_field(f, g);
        cmplx_timer.stop();
    }

    // Field<Vector<6,vtype>> v = 0,v1;

    // onsites(ALL) {
    //   v[X].e(1).gaussian();
    // }

    // hila::timer vector_timer("vector fft");
    // for (int i=0; i<loops; i++) {
    //   vector_timer.start();
    //   FFT_field(v,v1);
    //   vector_timer.stop();
    // }

    hila::finishrun();
}
