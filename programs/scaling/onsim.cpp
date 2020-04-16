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

inline double scaleFactor(double t, double t_end){
	return t/t_end;
}

class scaling_sim {

	public:
		scaling_sim() = default;
		void allocate(const std::string fname, int argc, char ** argv);
		void initialize();
		void write_energies();
		void next();
		inline double scaleFactor(double t);

		field<cmplx<double>> phi;
		field<cmplx<double>> pi;
		field<cmplx<double>> deltaPi; 

		double t;

		struct {
			int l;
			int seed;
			int smoothing;
			double param1;
			double sigma;
			double dx;
			double dt;
			double tStart;
			double tEnd;
			double lambda;
		} config;
};


void scaling_sim::allocate(const std::string fname, int argc, char ** argv){
	input parameters = input();
	parameters.import(fname);
	config.l = parameters.get("N");
	config.seed = parameters.get("seed");
	config.param1 = parameters.get("param1");
	config.sigma = parameters.get("sigma");
  	config.dx = parameters.get("dx");
	double ratio = parameters.get("dtdxRatio");
  	config.dt = config.dx*ratio;
  	config.tStart = parameters.get("tStart");
  	config.tEnd = parameters.get("tEnd");
	config.lambda = parameters.get("lambda");
	config.smoothing = parameters.get("smooth");

	int box_dimensions[3] = {config.l,config.l,config.l}; 

	lattice->setup(box_dimensions, argc, argv);
	seed_random(config.seed);
}

inline double scaling_sim::scaleFactor(double t){
	return t/config.tEnd;
}

void scaling_sim::initialize(){

	//initialize vaccuum state
	onsites(ALL){
		double theta, r;
		r = config.param1*config.sigma;
		theta = hila_random()*2*M_PI;
		cmplx<double> val;
		phi[X] = val.polar(r, theta);
		pi[X] = 0; 
	}

	//smoothing iterations
	for (int iter = 0; iter < config.smoothing; iter++){
		direction d;
		pi[ALL] = 6.0*phi[X];
		foralldir(d){
			pi[ALL] = pi[X] + phi[X + d];
		}
		onsites(ALL){
			cmplx<double> norm = pi[X].conj()*pi[X];
			if (norm.re == 0) norm = cmplx<double>(1.0, 0.0);
			phi[X] = pi[X]/norm;
			pi[X] = cmplx<double>(0.0, 0.0);
		}
	}

	

}

void scaling_sim::write_energies(){
	double a = scaleFactor(t);
	double ss = config.sigma*config.sigma;

	//non-weighted energies
    	double sumPi = 0.0;
   	double sumDiPhi = 0.0;
    	double sumV = 0.0;
	double sumPhiDiPhi = 0.0;
	double sumPhiPi = 0.0; 

	//weighted energies
        double w_sumPi = 0.0;
        double w_sumDiPhi = 0.0;
        double w_sumV = 0.0;
        double w_sumPhiDiPhi = 0.0; 
        double w_sumPhiPi = 0.0;

	onsites(ALL){
		double v = 0;
                cmplx<double> norm = phi[X].conj()*phi[X]; //calculate phi norm
                sumV += 0.25*config.lambda*a*a*pow((norm - ss).re, 2.0); //reduce potential
		sumPi += (pi[X].conj()*pi[X]).re; 
	}

	direction d;
	foralldir(d){
		onsites(ALL){
			cmplx<double> diff_phi = (phi[X + d] - phi[X])/config.dx;  
			double diPhi = (diff_phi.conj()*diff_phi).re;
			sumDiPhi += diPhi; //reduce diPhi
		}
	}

	if (mynode() == 0){
		double vol = (double) config.l*config.l*config.l;
		hila::output << sumV/vol << '\t' << sumDiPhi/vol << '\t' << sumPi/vol << '\n';
	}
}

void scaling_sim::next(){
	double a = scaleFactor(t);
	double aHalfPlus = scaleFactor(t + config.dt/2.0); 
	double aHalfMinus = scaleFactor(t - config.dt/2.0);
	double aadt_aadxdx = pow( a / aHalfPlus , 2.0) * config.dt / (config.dx*config.dx);
  	double aadt2D_aadxdx = aadt_aadxdx * 2.0 * NDIM;
 	double aaaaldt_aa = pow( a, 4.0 ) * config.lambda * config.dt / pow(aHalfPlus, 2.0);
  	double daa_aa = ( pow(aHalfPlus, 2.0) - pow(aHalfMinus, 2.0) ) / pow(aHalfPlus, 2.0);
  	double ss = config.sigma*config.sigma;

	phi[ALL] = phi[X] + config.dt*pi[X];

	onsites(ALL){
		cmplx<double> norm = phi[X].conj()*phi[X]; //calculate phi norm
		deltaPi[X] = -1.0*(aadt2D_aadxdx + aaaaldt_aa*(norm - ss))*phi[X]; 	
	}

	direction d;
	foralldir(d){
		deltaPi[ALL] += aadt_aadxdx*phi[X + d]; 
	}

	pi[ALL] = pi[X] - daa_aa*pi[X];
	pi[ALL] = pi[X] + deltaPi[X]; 

	t += config.dt;
}

int main(int argc, char ** argv){
	scaling_sim sim;
	sim.allocate("sim_params.txt", argc, argv);
	sim.initialize();
	while(sim.t < sim.config.tEnd){
		sim.next();
		sim.write_energies();
	}
	return 0;
}
