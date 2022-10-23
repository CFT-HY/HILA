#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <assert.h>

#include "plumbing/defs.h"
#include "datatypes/cmplx.h"
#include "datatypes/sun.h"
#if NDIM == 4
#include "datatypes/wilson_vector.h"
#endif
#include "plumbing/field.h"

// // Define the lattice global variable
// lattice_struct my_lattice;
// lattice_struct * lattice = & my_lattice;

const CoordinateVector nd = {32, 32, 32, 32};

inline void checkLatticeSetup() {
    for (int dir = 0; dir < NDIRS; dir++) {
        // check that neighbor arrays are allocated
        assert(lattice.neighb[dir] != nullptr);
#ifdef CUDA
        assert(lattice.backend_lattice->d_neighb[dir] != nullptr);
#endif
    }
    for (int dir = 0; dir < NDIM; dir++) {
        assert(lattice.size(dir) == nd[dir]);
    }
}

inline void test_setup(int &argc, char **argv) {
    hila::initialize(argc, argv);
    lattice.setup(nd);

    checkLatticeSetup();

    hila::seed_random(1);
}
