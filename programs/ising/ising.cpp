#include <iostream>
#include <string>
#include <math.h>

#define NDIM 2

// Include the lattice field definition
#include "../plumbing/field.h"

// Direct output to stdout
std::ostream &hila::output = std::cout;

// Define the lattice global variable
lattice_struct my_lattice;
lattice_struct * lattice = & my_lattice;


// Define some parameters for the simulation
double beta = 0.1;
int n_measurements=100;
int n_updates_per_measurement=10;
long seed = 123456;
int NX=64, NY=64;
int VOLUME = NX*NY;

int main()
{
  // Basic setup
  lattice->setup( NX, NY );
  // Define a field
  field<double> spin;

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
    for(int j=0; j<2*n_updates_per_measurement; j++){

      // A temporary field for the local change in action
      onsites(p){
        double deltaS;
        deltaS = 2.0*spin[X]*( spin[X+XUP] + spin[X+XDOWN]
                             + spin[X+YUP] + spin[X+YDOWN] );

        if( mersenne() < exp(-beta*deltaS) ){
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
    printf("Magnetisation %f\n", M/VOLUME);
  }
  

  return 0;
}
