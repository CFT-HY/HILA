#include <sstream>
#include<iostream>
#include <string>
#include <math.h>

// Include the lattice field definition
#include "plumbing/field.h"


// Define some parameters for the simulation
double beta = 0.1;
int n_measurements=100;
int n_updates_per_measurement=10;
long seed = 123456;
int NX=64, NY=64;
int VOLUME = NX*NY;

int main(int argc, char **argv)
{
  // Basic setup
  const int nd[2] = { NX, NY };
  hila::initialize(argc,argv);
  lattice->setup(nd);

  // Define a field
  field<double> spin;

  seed_random( seed );
  
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
        deltaS = 2.0*spin[X]*( spin[X+e_x] + spin[X-e_x]
                             + spin[X+e_y] + spin[X-e_y] );

        if( hila_random() < exp(-beta*deltaS) ){
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
    output0 << "Magnetisation " << M/VOLUME << "\n";
  }
  
  hila::finishrun();
  return 0;
}
