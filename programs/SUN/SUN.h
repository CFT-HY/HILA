#ifndef SUN_H_
#define SUN_H_

#include <iostream>
#include <string>
#include <math.h>

#define PI 3.14159265358979323846

#define NDIM 4
const int N=2;

// Include the lattice field definition
#include "../plumbing/field.h"
#include "../datatypes/general_matrix.h"


// Define some parameters for the simulation
extern double beta;
extern int n_measurements;
extern int n_updates_per_measurement;
extern long seed;
extern int NX, NY, NZ, NT;
extern int VOLUME;

loop_callable double monte(
  matrix<N,N,cmplx<double>> &U, 
  matrix<N,N,cmplx<double>> &staple,
  double beta);

loop_callable void KennedyPendleton(
  matrix<2,2,cmplx<double>> &U,
  matrix<2,2,cmplx<double>> &staple
);


#endif //SUN_H_