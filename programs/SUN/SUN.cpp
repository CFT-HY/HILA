#include <iostream>
#include <string>
#include <math.h>

#define NDIM 4

// Include the lattice field definition
#include "../plumbing/field.h"
#include "../datatypes/general_matrix.h"
extern "C"
{
    #include "mersenne.h"
};

// Direct output to stdout
std::ostream &hila::output = std::cout;
std::ostream &output = std::cout;

// Define the lattice global variable
lattice_struct my_lattice;
lattice_struct * lattice = & my_lattice;


// Define some parameters for the simulation
double beta = 0.4;
int n_measurements=100;
int n_updates_per_measurement=1;
long seed = 123456;
int NX=8, NY=8, NZ=8, NT=8;
int VOLUME = NX*NY*NZ*NT;

const int N=2;

int main()
{
  // Basic setup
  lattice->setup( NX, NY, NZ, NT );
  // Define a field
  field<matrix<N,N,double>> U[NDIM];

  seed_mersenne( seed );

  /* "Warm up" the rng generator */
  for( int i=0; i<543210; i++ ) mersenne();
  
  // Set to 1
  foralldir(d) {
    printf(" %d\n", d);
    U[d][ALL] = 1;
  }

  // Run update-measure loop
  for( int i=0; i<n_measurements; i++ ){

    // Run a number of updates, starting with EVEN sites
    // and alternating between parities
    parity p = EVEN;
    for(int j=0; j<2*n_updates_per_measurement; j++){
      p=opp_parity(p);
    }

    // Measure plauqette
    double Plaq=0;
    foralldir(d1) foralldir(d2) if(d1 != d2){
      direction dir1 = (direction)d1, dir2 = (direction)d2;
      onsites(ALL){
        matrix<N,N,double> temp;
        temp =  U[dir1][X] * U[dir2][X+dir1];
        temp *= U[dir2][X+dir1].conjugate();
        temp *= U[dir1][X+dir2].conjugate();
        Plaq += 1-temp.trace()/N;
      }
    }
    printf("Plaquette %f\n", Plaq/VOLUME);
  }
  

  return 0;
}
