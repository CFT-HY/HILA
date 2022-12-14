#ifndef WILSON_H_
#define WILSON_H_

#include <iostream>
#include <string>
#include <math.h>

#define PI 3.14159265358979323846
#define NDIM 4

// Include the lattice field definition
#include "plumbing/defs.h"
#include "datatypes/matrix.h"
#include "datatypes/sun.h"
#include "datatypes/representations.h"
#include "datatypes/wilson_vector.h"
#include "plumbing/field.h"
#include "hmc/hmc.h"
#include "hmc/gauge_field.h"
#include "dirac/wilson.h"
#include "hmc/fermion_field.h"
#include "plumbing/input.h"

const int N = 2;

using SUN = SU<N, double>;
using NMAT = Matrix<N, N, Complex<double>>;
using VEC = SU_vector<N, double>;

// Define some parameters for the simulation
extern const CoordinateVector nd = {8, 8, 8, 8};

extern int n_measurements;
extern int n_updates_per_measurement;
extern long seed;

#endif // STAGGERED_H_
