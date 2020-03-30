/*********************************************
 * Simulates staggered fermions with a gauge *
 * interaction                               *
 *********************************************/


#include "staggered.h"




int main(int argc, char **argv){

  lattice->setup( nd[0], nd[1], nd[2], nd[3], argc, argv );

  field<SU<N,double>> gauge;
  field<matrix<N,N,cmplx<double>>> gauge_m;
  field<SU_vector<N,double>> vector;

  // Starting with a unit configuration
  gauge_m[ALL] = 1.0;

  return 0;
}