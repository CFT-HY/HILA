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

class filewriter {

};

class scaling_sim {

	public:

		scaling_sim() = default;
		void allocate(const std::string fname, int argc, char ** argv);
		void initialize();
		void write_moduli();
		void write_energies();
		void next();
		inline double scaleFactor(double t);

		field<cmplx<double>> phi;
		field<cmplx<double>> pi;
		field<cmplx<double>> deltaPi; 

		double t;

		struct {
			int l;
			int m;
			int seed;
			int smoothing;
			int initalCondition;
			double initialModulus;
			double epsilon;
			double sigma;
			double dx;
			double dt;
			double tStart;
			double tStats;
			double nOutputs;
			double tEnd;
			double lambda;
			std::fstream stream;
		} config;
};

void scaling_sim::allocate(const std::string fname, int argc, char ** argv){
	input parameters = input();
	parameters.import(fname);
	config.l = parameters.get("N");
	config.m = parameters.get("m");
	config.epsilon = parameters.get("epsilon");
	config.seed = parameters.get("seed");
	config.initialModulus = parameters.get("initialModulus");
	config.sigma = parameters.get("sigma");
  	config.dx = parameters.get("dx");
  	config.tStart = parameters.get("tStart");
  	config.tEnd = parameters.get("tEnd");
	config.lambda = parameters.get("lambda");
	config.smoothing = parameters.get("smooth");
	config.initalCondition = parameters.get("initialCondition");
	config.tStats = parameters.get("tStats");
	config.nOutputs = parameters.get("numberStatsOutputs");
	double ratio = parameters.get("dtdxRatio");
  	config.dt = config.dx*ratio;
	t = config.tStart;

	int box_dimensions[3] = {config.l,config.l,config.l}; 

	lattice->setup(box_dimensions, argc, argv);
	seed_random(config.seed);
}

inline double scaling_sim::scaleFactor(double t){
	return t/config.tEnd;
}

void scaling_sim::initialize(){

	int m =config.m;
	double epsilon = config.epsilon;
	double s = config.sigma;
	int N = config.l;

	switch(config.initalCondition){

		case 2 : {
			onsites(ALL){
				pi[X] = cmplx<double>(0.0, 0.0);
				phi[X].re = config.sigma;
				phi[X].im = config.sigma;
			}

			if (mynode() == 0){
					hila::output << "Field real and imaginary components set to sigma = "<< config.sigma << '\n';
			}

			break;
		}

		case 1 : {
			onsites(ALL){
				element<coordinate_vector> coord = coordinates(X); //coordinates of the current lattice site
				pi[X] = cmplx<double>(0.0, 0.0);
				phi[X].re = s*sqrt(1-epsilon*epsilon*sin(2.0*M_PI*coord[0]*m/N)*sin(2.0*M_PI*coord[0]*m/N));
				phi[X].im = s*epsilon*sin(2.0*M_PI*coord[0]*m/N);
			}

			if (mynode() == 0){
					hila::output << "Axion wave generated with amplitude: "<< config.epsilon << '\n';
			}

			break;
		}

		default : { 
			
			onsites(ALL){
				double theta, r;
				r = config.initialModulus*s;
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

			if (mynode() == 0){
					hila::output << "Field elements randomly generated \n";
			}

			break;
		}
	}
}

void scaling_sim::write_moduli(){

	double a = scaleFactor(t);

	double phimod = 0.0;
	double pimod = 0.0;

	onsites(ALL){
		double p_r = 0.0, p_i = 0.0;
		cmplx<double> norm_1 = phi[X].conj()*phi[X];
		cmplx<double> norm_2 = pi[X].conj()*pi[X]; 
		p_r = norm_1.re;
		p_i = norm_2.re;
		phimod += sqrt(p_r);
		pimod += sqrt(p_i);
	}

	if (mynode() == 0){
		double vol = (double) (config.l*config.l*config.l);
		config.stream << t << "," << a << "," << config.lambda << "," << phimod/vol << "," << pimod/vol << ",";
	}
}

void scaling_sim::write_energies(){

	double a = scaleFactor(t);
	double ss = config.sigma*config.sigma;

	//non-weighted energies
    double sumPi = 0.0; // 
   	double sumDiPhi = 0.0; //
    double sumV = 0.0; //
	double sumPhiDiPhi = 0.0;
	double sumPhiPi = 0.0; 

	//weighted energies
    double w_sumPi = 0.0;
    double w_sumDiPhi = 0.0;
    double w_sumV = 0.0;
    double w_sumPhiDiPhi = 0.0; 
    double w_sumPhiPi = 0.0;

	onsites(ALL){
        double phinorm = (phi[X].conj()*phi[X]).re; 
		double v = 0.25*config.lambda*a*a*pow((phinorm - ss), 2.0);
		double pinorm = (pi[X].conj()*pi[X]).re;
		double pPi = (phi[X].conj()*pi[X]).re;

        sumV += v;
		w_sumV += v*v;

		sumPi += 0.5*pinorm; 
		w_sumPi += 0.5*pinorm*v;

		sumPhiPi += 0.5*pPi*pPi;
		w_sumPhiPi += 0.5*pPi*pPi*v;
	}

	direction d;
	foralldir(d){
		onsites(ALL){
			cmplx<double> norm = phi[X].conj()*phi[X];
			double v = 0.25*config.lambda*a*a*pow((norm.re - ss), 2.0);
			cmplx<double> diff_phi = (phi[X + d] - phi[X])/config.dx;
			double pDphi = 0.5*(diff_phi.conj()*phi[X]).re;

			sumDiPhi += 0.5*(diff_phi.conj()*diff_phi).re; 
			sumPhiDiPhi += pDphi*pDphi/norm.re;

			w_sumDiPhi += 0.5*(diff_phi.conj()*diff_phi).re*v; 
			w_sumPhiDiPhi += pDphi*pDphi/norm.re*v;
		}
	}

	if (mynode() == 0){
		double vol = (double) config.l*config.l*config.l;
		config.stream << sumPi/vol << "," << w_sumPi/vol << ",";
		config.stream << sumDiPhi/vol << "," << w_sumDiPhi/vol << ","; 
		config.stream << sumPhiPi/vol << "," << w_sumPhiPi/vol << ","; 
		config.stream << sumPhiDiPhi/vol << "," << w_sumPhiDiPhi/vol << ","; 
		config.stream << sumV/vol << "," << w_sumV/vol << "\n"; 
	}
}

void scaling_sim::next(){

	double a = scaleFactor(t);
	double aHalfPlus = scaleFactor(t + config.dt/2.0); 
	double aHalfMinus = scaleFactor(t - config.dt/2.0);

	double aadt_aadxdx = pow( a / aHalfPlus , 2.0) * config.dt / (config.dx*config.dx);
  	double aadt2D_aadxdx = aadt_aadxdx * 2.0 * 3.0;
 	double aaaaldt_aa = pow( a, 4.0 ) * config.lambda * config.dt / pow(aHalfPlus, 2.0);
  	double daa_aa = ( pow(aHalfPlus, 2.0) - pow(aHalfMinus, 2.0) ) / pow(aHalfPlus, 2.0);
  	double ss = config.sigma*config.sigma;

	phi[ALL] = phi[X] + config.dt*pi[X];

	onsites(ALL){
		cmplx<double> norm = phi[X].conj()*phi[X]; //calculate phi norm
		deltaPi[X] = phi[X]*(aaaaldt_aa*(ss - norm.re) - aadt2D_aadxdx);
	}

	direction d;
	forALLdir(d){
		onsites(ALL){
			deltaPi[X] = deltaPi[X] + aadt_aadxdx*phi[X + d];
		}
	}

	pi[ALL] = pi[X] - daa_aa*pi[X];
	pi[ALL] = pi[X] + deltaPi[X]; 

	t += config.dt;

}

int main(int argc, char ** argv){
	scaling_sim sim;

	std::string input_fname(argv[1]);
	std::string output_fname(argv[2]);

	sim.allocate(input_fname, argc, argv);
	sim.initialize();

	int steps = (sim.config.tEnd - sim.config.tStats)/(sim.config.dt * sim.config.nOutputs); //number of steps between printing stats 
	int stat_counter = 0;

	if (mynode() == 0){
		sim.config.stream.open(output_fname, std::ios::out);
	}

	while (sim.t < sim.config.tEnd){
		if (sim.t >= sim.config.tStats){
			if (stat_counter%steps==0){
				synchronize();
				sim.write_moduli();
				sim.write_energies();
			}
			stat_counter++;
		}
		sim.next();
	}

	if (mynode() == 0){
		sim.config.stream.close();
	}

	finishrun();
	return 0;
}
