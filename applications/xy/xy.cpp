#include <sstream>
#include<iostream>
#include <string>
#include <math.h>

#define NDIM 2

// Include the lattice field definition
#include "plumbing/field.h"

// Direct output to stdout
std::ostream &hila::output=std::cout;

// Define the lattice global variable
lattice_struct my_lattice;
lattice_struct * lattice=&my_lattice;


// Define some parameters for the simulation
double beta=0.1;
int n_measurements=100;
int n_updates_per_measurement=10;
long seed=123456;
int NX=64,NY=64;
int VOLUME=NX*NY;

int main(int argc, char **argv){
	// Basic setup
	lattice->setup(NX,NY, argc, argv);
	// Define a field
	field<double> spin;

	seed_random(seed);

	/* "Warm up" the rng generator */
	for(int i=0; i<543210; i++) hila_random();

	// Set to 1
	spin[ALL]=0;

	// Run update-measure loop
	for(int i=0; i<n_measurements; i++) {

		// Run a number of updates, starting with EVEN sites
		// and alternating between parities
		parity p=EVEN;
		for(int j=0; j<2*n_updates_per_measurement; j++) {
			onsites(p) {
				double deltaS;
				double tspin=spin[X];
				double tnspin=tspin+M_PI*(1.-2.*hila_random());
				deltaS=cos(spin[X+XUP]-tspin)+cos(spin[X+XDOWN]-tspin)+cos(spin[X+YUP]-tspin)+cos(spin[X+YDOWN]-tspin);
				deltaS-=cos(spin[X+XUP]-tnspin)+cos(spin[X+XDOWN]-tnspin)+cos(spin[X+YUP]-tnspin)+cos(spin[X+YDOWN]-tnspin);

				if(deltaS<0 || hila_random()<exp(-beta*deltaS)) {
					if(tnspin<-M_PI) {
						tnspin+=2.0*M_PI;
					} else if(tnspin>M_PI) {
						tnspin-=2.0*M_PI;
					}
					spin[X]=tnspin;
				}
			}
			if(hila_random()<0.5) {
				p=opp_parity(p);
			}
		}

		// Measure magnetisation
		double M=0;
		onsites(ALL) {
			M+=cos(spin[X]);
		}
		output0 << "Magnetisation " << M/VOLUME << "\n";
	}

	finishrun();
	return 0;
}
