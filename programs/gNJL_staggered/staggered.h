#ifndef STAGGERED_H_
#define STAGGERED_H_

#include <iostream>
#include <string>
#include <math.h>

#define PI 3.14159265358979323846
#define NDIM 4

// Include the lattice field definition
#include "../plumbing/defs.h"
#include "../plumbing/inputs.h"
#include "../datatypes/sun_vector.h"
#include "../datatypes/sun.h"
#include "../plumbing/field.h"
#include "dirac_staggered_NJL.h"
#include "../plumbing/algorithms/hmc.h"
#include "../plumbing/gauge_field.h"
#include "../plumbing/fermion_field.h"


const int N=2;

using SUN = SU<N,double>;
using NMAT = matrix<N,N,cmplx<double>>;
using VEC = SU_vector<N,double>;

// Direct output to stdout
std::ostream &hila::output = std::cout;

// Define the lattice global variable
lattice_struct my_lattice;
lattice_struct * lattice = & my_lattice;



// Define some parameters for the simulation
extern const int nd[4] = { 8, 8, 8, 8 };

extern int n_measurements;
extern int n_updates_per_measurement;
extern long seed;


#endif //STAGGERED_H_