#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <assert.h>

#include "plumbing/defs.h"
#include "datatypes/cmplx.h"
// #include "datatypes/sun.h"
#if NDIM == 4
// #include "datatypes/wilson_vector.h"
#endif
#include "plumbing/field.h"

// // Define the lattice global variable
// lattice_struct my_lattice;
// lattice_struct * lattice = & my_lattice;

const CoordinateVector nd = {32, 32, 32, 32};


inline void test_setup(int &argc, char **argv) {
    hila::initialize(argc, argv);
    lattice.setup(nd);

    hila::seed_random(1);
}
