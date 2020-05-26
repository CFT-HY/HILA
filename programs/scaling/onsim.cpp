#define _USE_MATH_DEFINES
#include <sstream>
#include<iostream>
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

inline double scaleFactor(double t, double t_end){
	return t/t_end;
}

int main(int argc, char ** argv){

	input parameters = input();
	parameters.import("sim_params.txt");

	int l = parameters.get("l");
	int seed = parameters.get("seed");
	double param1 = parameters.get("param1");
	double sigma = parameters.get("sigma");
  	double dx = parameters.get("dx");
  	double dt = parameters.get("dt");
  	double tStart = parameters.get("tStart");
  	double tEnd = parameters.get("tEnd");
	double lambda = parameters.get("lambda");

	//parameters.close();
 
	int box_dimensions[3] = {l,2*l,4*l}; 

	lattice->setup(box_dimensions, argc, argv);
	seed_random(seed);

	//define fields

	field<cmplx<double>> phi;
	field<cmplx<double>> pi;
	field<cmplx<double>> deltaPi; 
	field<cmplx<double>> e;  

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

	for (double t = tStart; t < tEnd; t += dt){

		double a = scaleFactor(t, tEnd);  
      		double aHalfPlus = scaleFactor(t+dt/2.0, tEnd);
      		double aHalfMinus = scaleFactor(t-dt/2.0, tEnd);

		double aadt_aadxdx = pow( a / aHalfPlus , 2.0) * dt / (dx*dx);
  		double aadt2D_aadxdx = aadt_aadxdx * 2.0 * NDIM;
 		double aaaaldt_aa = pow( a, 4.0 ) * lambda * dt / pow(aHalfPlus, 2.0);
  		double daa_aa = ( pow(aHalfPlus, 2.0) - pow(aHalfMinus, 2.0) ) / pow(aHalfPlus, 2.0);
  		double ss = sigma*sigma;
		
		onsites(ALL){ 
			cmplx<double> mod = phi[X].conj()*phi[X];
			deltaPi[X] = -1.0*(aadt2D_aadxdx + aaaaldt_aa*(mod - ss))*phi[X];
		}

		direction d;
		foralldir(d){
			deltaPi[ALL] += aadt_aadxdx*phi[X + d]; 
		}

		pi[ALL] = pi[X] + deltaPi[X];
	}
	return 0;
}
