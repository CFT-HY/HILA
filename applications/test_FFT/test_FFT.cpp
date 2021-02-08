
#include "plumbing/defs.h"
#include "datatypes/cmplx.h"
#include "datatypes/matrix.h"
#include "plumbing/field.h"
#include "plumbing/param_input.h"
#include "plumbing/FFT_new.h"

using vtype = Cmplx<double>;

static_assert( NDIM == 3, "NDIM must be 3");


int main(int argc, char **argv){

    hila::initialize(argc,argv);

    input par("parameters");
    int nx = par.get("nx");
    int ny = par.get("ny");
    int nz = par.get("nz");
    int loops = par.get("loops");
    int seed = par.get("random seed");
    
    lattice->setup(nx,ny,nz);

    
    par.close();

    
    seed_random(seed);

    Field<vtype> f,g;
    
    f[ALL] = gaussian_ran();

    static timer loop_timer("fft loops"); 
    
    for (int i=0; i<loops; i++) {
      loop_timer.start();
      FFT_field(f,g);
      loop_timer.stop();
      g = f;
    }

    hila::finishrun();
}



