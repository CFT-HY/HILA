#include <iostream>
#include <string>
#include <math.h>
#include <assert.h>

#include "../plumbing/defs.h"
#include "../datatypes/cmplx.h"
#include "../datatypes/general_matrix.h"
#include "../datatypes/sun.h"
#include "../plumbing/field.h"
#include "../plumbing/inputs.h"

// Direct output to stdout
std::ostream &hila::output = std::cout;

// Define the lattice global variable
lattice_struct my_lattice;
lattice_struct * lattice = & my_lattice;

//define lattice dimensions 
int size[3] = {30, 30, 30};

int main(int argc, char ** argv){

	input a;  
	
	lattice->setup(size, argc, argv);
	seed_random(20);

	field<cmplx<double>> phi;
	field<cmplx<double>> pi;

	onsites(ALL){
		phi[X] = hila_random();
		pi[X] = hila_random(); 
	}
	return 0;
}
