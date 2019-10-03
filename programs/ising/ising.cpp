#include <iostream>
#include <string>
#include <math.h>

#define NDIM 2

// Include the lattice field definition
#include "../plumbing/field.h"
#include "../datatypes/scalar.h"
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

double beta = 0.1;
int n_measurements=100;
int n_updates_per_measurement=10;
long seed = 123456;

int main()
{
  // Basic setup
  lattice->setup( 8, 8 );

  // Define a field
  field<scalar<double>> spin;

  seed_mersenne( seed );

  /* "Warm up" the rng generator */
  for( int i=0; i<543210; i++ ) mersenne();
  
  // Set to 1
  spin[ALL] = 1;

  // Run update-measure loop
  for( int i=0; i<n_measurements; i++ ){

    // Run a number of updates, starting with EVEN sites
    // and alternating between parities
    parity p = EVEN;
    for(int i=0; i<n_updates_per_measurement; i++){
      onsites(p){
        double S, neighbour_sum=0;
        neighbour_sum += spin[X+XUP];
        neighbour_sum += spin[X+XDOWN];
        neighbour_sum += spin[X+YUP];
        neighbour_sum += spin[X+YDOWN];

        S = 2*neighbour_sum*(spin[X]);
        if( mersenne() < exp(-beta*S) ){
          spin[X] = -spin[X];
        }
      }
      p=opp_parity(p);
    }

    // Measure magnetisation
    double M=0;
    onsites(ALL){
      M += spin[X];
    }

    printf("Magnetisation %f\n", M);
  }
  
  return 0;
}
