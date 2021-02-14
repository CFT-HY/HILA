#ifndef STAGGERED_H_
#define STAGGERED_H_

#include <iostream>
#include <math.h>
#include <string>

#define PI 3.14159265358979323846
#define NDIM 4

// Include the lattice field definition
#include "plumbing/defs.h"
#include "datatypes/matrix.h"
#include "datatypes/sun.h"
#include "datatypes/representations.h"
#include "plumbing/field.h"
#include "hmc/hmc.h"
#include "hmc/gauge_field.h"
#include "dirac/staggered.h"
#include "hmc/fermion_field.h"
#include "plumbing/input.h"

const int N = 3;

using SUN = SU<N, double>;
using NMAT = Matrix<N, N, Cmplx<double>>;
using VEC = SU_vector<N, double>;

// Define some parameters for the simulation
extern const CoordinateVector nd{8, 8, 8, 8};

extern int n_measurements;
extern int n_updates_per_measurement;
extern long seed;

#endif // STAGGERED_H_
