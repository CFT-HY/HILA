#ifndef SUN_H_
#define SUN_H_

#include <sstream>
#include<iostream>
#include <string>
#include <math.h>

#define PI 3.14159265358979323846

#define NDIM 4
const int N=2;

// Include the lattice field definition
#include "plumbing/defs.h"
#include "datatypes/matrix.h"
#include "plumbing/field.h"


// Define some parameters for the simulation
extern double beta;
extern int n_measurements;
extern int n_updates_per_measurement;
extern long seed;
extern int NX, NY, NZ, NT;
extern int VOLUME;

//#pragma hila loop_function  //TODO
double monte(
  Matrix<N,N,Cmplx<double>> &U, 
  Matrix<N,N,Cmplx<double>> &staple,
  double beta);

//#pragma hila loop_function  //TODO
void KennedyPendleton(
  Matrix<2,2,Cmplx<double>> &U,
  Matrix<2,2,Cmplx<double>> &staple
);


#endif //SUN_H_