
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
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

inline double scalefactor(double t, double t_end){
	return t/t_end;
}

int main(int argc, char ** argv){

	int l = 10;
	int seed = 20;
	int steps = 30;
	double param1 = 1.0;
	double sigma = 0.1;

  	int dim = 3;
  	dx = 1.0;
  	dt = 0.1;
  	tStart = 0.0;
  	tEnd = 32.0;
 
	int box_dimensions[3] = {l,2*l,4*l}; 

	lattice->setup(box_dimensions, argc, argv);
	seed_random(seed);

	//define fields

	field<cmplx<double>> phi;
	field<cmplx<double>> pi;
	field<cmplx<double>> e; //energy 

	//initialize vaccuum state phi, set pi to zero
 
	onsites(ALL){
		double theta, r;
		r = param1*sigma;
		theta = hila_random()*2*M_PI;
		cmplx<double> val;
		phi[X] = val.polar(r, theta);
		pi[X] = 0; 
	}

	//evolve fields

	double t = 0;
	for (int i = 0; i < steps; i++){

		double a = scaleFactor(t);  
      		double aHalfPlus = scaleFactor(t+dt/2.0);
      		double aHalfMinus = scaleFactor(t-dt/2.0);

		double aadt_aadxdx = pow( a / aHalfPlus , Real(2) ) * dt / (dx*dx);
  		double aadt2D_aadxdx = aadt_aadxdx * Real(2) * dim;
 		double aaaaldt_aa = pow( a, Real(4) ) * lambda * dt / pow(aHalfPlus, 2);
  		double daa_aa = ( pow(aHalfPlus, 2) - pow(aHalfMinus, 2) ) / pow(aHalfPlus, 2);
  		double ss = sigma*sigma;
		double*  deltaPi= new double[nf]; //delta-pi at each site
    		double mod; //at each site 
	}

	return 0;
}
