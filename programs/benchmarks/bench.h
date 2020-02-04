#include <iostream>
#include <string>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "../plumbing/defs.h"
#include "../datatypes/general_matrix.h"
#include "../datatypes/sun.h"
#include "../plumbing/field.h"
#include "../plumbing/dirac.h"

// Direct output to stdout
std::ostream &hila::output = std::cout;

// Define the lattice global variable
lattice_struct my_lattice;
lattice_struct * lattice = & my_lattice;

const int nd[4] = { 32, 32, 32, 32 };

inline void bench_setup(int &argc, char **argv){
    #if NDIM==1
    lattice->setup( nd[0], argc, argv );
    #elif NDIM==2
    lattice->setup( nd[0], nd[1], argc, argv );
    #elif NDIM==3
    lattice->setup( nd[0], nd[1], nd[2], argc, argv );
    #elif NDIM==4
    lattice->setup( nd[0], nd[1], nd[2], nd[3], argc, argv );
    #endif
}
