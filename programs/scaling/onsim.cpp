
#define _USE_MATH_DEFINES
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

int main(int argc, char ** argv){

	input parameters = input(); 
	parameters.add_essential("l"); //define parameters that must be in parameter file
	parameters.add_essential("sigma");
	parameters.add_essential("seed");
	parameters.add_essential("param1");
	parameters.import("sim_params.txt"); //import sim_params
  
	int l = parameters.get("l");
	int seed = parameters.get("seed");
	
	double param1 = parameters.get("param1");
	double sigma = parameters.get("sigma");

	int box_dimensions[3] = {l,l,l}; 

	lattice->setup(box_dimensions, argc, argv);
	seed_random(seed);

	//define fields
	field<cmplx<double>> phi;
	field<cmplx<double>> pi;

	//initialize vaccuum state phi, set pi to zero
 
	onsites(ALL){
		double theta, r;
		r = param1*sigma;
		theta = hila_random()*2*M_PI;
		cmplx<double> val;
		phi[X] = val.polar(r, theta);
		pi[X] = 0; 
	}

	return 0;
}
