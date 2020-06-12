#ifndef STAGGERED_H_
#define STAGGERED_H_

#include <iostream>
#include <string>
#include <math.h>

#define PI 3.14159265358979323846
#define NDIM 4

// Include the lattice field definition
#include "../../libraries/plumbing/defs.h"
#include "../../libraries/plumbing/inputs.h"
#include "../../libraries/datatypes/vector.h"
#include "../../libraries/datatypes/sun.h"
#include "../../libraries/datatypes/representations.h"
#include "../../libraries/datatypes/wilson_vector.h"
#include "../../libraries/plumbing/field.h"
#include "../../libraries/dirac/wilson.h"
#include "../../libraries/hmc/update.h"
#include "../../libraries/hmc/gauge_field.h"
#include "../../libraries/hmc/fermion_field.h"


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
