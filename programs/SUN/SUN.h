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
extern "C"
{
    #include "mersenne.h"
};


// Define some parameters for the simulation
extern double beta;
extern int n_measurements;
extern int n_updates_per_measurement;
extern long seed;
extern int NX, NY, NZ, NT;
extern int VOLUME;


double monte(
  matrix<N,N,cmplx<double>> &U, 
  matrix<N,N,cmplx<double>> &staple,
  double beta);



#endif //SUN_H_