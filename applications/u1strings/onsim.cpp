#define _USE_MATH_DEFINES
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
//#include <math.h>
#include <assert.h>

#include "plumbing/hila.h"
#include "plumbing/fft.h"

using real_t = float;

// inline double scaleFactor(double t, double t_end) {
//     return t / t_end;
// }

/// Container for simulation parameters and methods
class scaling_sim {

  public:
    scaling_sim() = default;
    const std::string allocate(const std::string &fname, int argc, char **argv);
    void initialize();
    void write_moduli();
    void write_energies();
    void write_windings();
    void next();
    void updateCosmology(); 
    inline real_t scaleFactor(real_t t);

    Field<Complex<real_t>> phi;
    Field<Complex<real_t>> pi;

    real_t t;
    real_t a;
    real_t aHalfPlus;
    real_t aHalfMinus;
    real_t lambda;

    real_t adif;
    real_t lambdadif;
    real_t acg;
    real_t lambdacg;

    struct config {
        int l;
        int m;
        int seed;
        int smoothing;
        int initalCondition;
        real_t initialModulus;
        real_t PhiLength;
        real_t epsilon;
        real_t sigma;
        real_t dx;
        real_t dt;
	real_t era;
        real_t tStart;
        real_t tdif;
        real_t difFac;
        real_t tcg;
	real_t s1;
	real_t s2;
   	real_t tStats;
        real_t nOutputs;
        real_t tEnd;
        real_t lambda0;
        std::fstream stream;
    } config;
};

const std::string scaling_sim::allocate(const std::string &fname, int argc,
                                        char **argv) {
    hila::initialize(argc, argv);

    hila::input parameters(fname);
    config.l = parameters.get("N");
    config.m = parameters.get("m");
    config.epsilon = parameters.get("epsilon");
    config.seed = parameters.get("seed");
    config.initialModulus = parameters.get("initialModulus");
    config.sigma = parameters.get("sigma");
    config.dx = parameters.get("dx");
    config.era = parameters.get("era");
    config.tStart = parameters.get("tStart");
    config.tEnd = parameters.get("tEnd");
    config.lambda0 = parameters.get("lambda");
    config.tdif = parameters.get("tdif");
    config.difFac = parameters.get("difFac");
    config.tcg = parameters.get("tcg");
    config.s1 = parameters.get("s1");
    config.s2 = parameters.get("s2");
    config.smoothing = parameters.get("smooth");
    config.initalCondition = parameters.get("initialCondition");
    config.PhiLength = parameters.get("PhiLength");
    config.tStats = parameters.get("tStats");
    config.nOutputs = parameters.get("numberStatsOutputs");
    real_t ratio = parameters.get("dtdxRatio");
    const std::string output_file = parameters.get("output_file");
    config.dt = config.dx * ratio;
    t = config.tStart;

    if (config.tcg<config.tdif)
    {
	output0 << "Core growth time has to be bigger or equal to Diffusion time"<< '\n';
	output0 << "Setting Core growth time equal to Diffusion time"<< '\n';
	config.tcg=config.tdif;		
    }
 
    acg = pow(config.tcg/config.tEnd, config.era);
    adif = acg * pow (config.tdif/config.tcg, config.era);

    lambdacg = config.lambda0 * pow(acg, -2.0*(1.0-config.s2));
    lambdadif = lambdacg * pow( adif/acg, -2.0*(1.0-config.s1));
     
    CoordinateVector box_dimensions = {config.l, config.l, config.l};
    lattice->setup(box_dimensions);
    hila::seed_random(config.seed);

    return output_file;
}

inline real_t scaling_sim::scaleFactor(real_t t) {
    return pow(t / config.tEnd, config.era);
}

void scaling_sim::updateCosmology() { 


    if (t <= config.tdif)
    {
        a = adif;
        aHalfPlus = adif;
        aHalfMinus = adif;
    	lambda = lambdadif;
    }
    else if (t > config.tdif && t<= config.tcg)
    {
	a = acg * pow(t/config.tcg, config.era);
    	aHalfPlus = acg * pow((t+0.5*config.dt)/config.tcg, config.era);
    	aHalfMinus = acg * pow((t-0.5*config.dt)/config.tcg, config.era);
	lambda = lambdacg * pow (a/acg, -2.0*(1.0-config.s1));
    }    
    else
    {
	a = scaleFactor(t);
	aHalfPlus = scaleFactor(t+0.5*config.dt);
	aHalfMinus = scaleFactor(t-0.5*config.dt);
	lambda = lambdacg * pow(a/acg, -2.0*(1.0-config.s2));	
    }
}

void scaling_sim::initialize() {

    int m = config.m;
    real_t epsilon = config.epsilon;
    real_t s = config.sigma;
    int N = config.l;
    real_t dx = config.dx;
    
    switch (config.initalCondition) {

    case 2: {
        pi = 0;
        phi = Complex<real_t>(config.sigma, config.sigma);

        output0 << "Field real and imaginary components set to sigma = " << config.sigma
                << '\n';

        break;
    }

    case 1: {
        pi = 0;
        onsites (ALL) {
            auto xcoord = X.coordinate(e_x);
            phi[X].re =
                s * sqrt(1 - epsilon * epsilon * sin(2.0 * M_PI * xcoord * m / N) *
                                 sin(2.0 * M_PI * xcoord * m / N));
            phi[X].im = s * epsilon * sin(2.0 * M_PI * xcoord * m / N);
        }

        output0 << "Axion wave generated with amplitude: " << config.epsilon << '\n';

        break;
    }

    case 3: {
      auto kphi = phi;
      
      onsites (ALL) {
	real_t constant = pow(config.initialModulus,2.0)*pow(2.0*M_PI,1.5)*pow(config.PhiLength,3.0)/(2.0*N*N*N*dx*dx*dx);
	real_t kSqu;
	real_t std;
	kSqu = 0.0;
	auto k = X.coordinates();

	foralldir (d) {
	  kSqu += pow( sin(M_PI*k.e(d)/N), 2.0);
	}
	kSqu *= pow(2.0/dx, 2.0);

	if (kSqu > 0.0) {
	  std = sqrt(0.5*constant*exp(-0.5*kSqu*config.PhiLength*config.PhiLength));
	  kphi[X].re = hila::gaussrand()*std;
	  kphi[X].im = hila::gaussrand()*std;
	}
	else {
	  kphi[X].re = 0.0;
          kphi[X].im = 0.0;
	}	
      }

      FFT_field(kphi, phi, fft_direction::back);

      pi[ALL] = 0;

      output0 << "k space generation \n";    

      break;

    }
      
    default: {

        // #pragma hila ast_dump
        onsites (ALL) {
            real_t theta, r;
            r = config.initialModulus * s;
            theta = hila::random() * 2 * M_PI;
            phi[X].polar(r, theta);
            pi[X] = 0;
        }

        // smoothing iterations
        for (int iter = 0; iter < config.smoothing; iter++) {
            pi[ALL] = 6.0 * phi[X];
            foralldir (d) { pi[ALL] += phi[X + d] + phi[X - d]; }
            onsites (ALL) {
                phi[X] = pi[X] / pi[X].abs();
                pi[X] = 0;
            }
        }

        output0 << "Field elements randomly generated \n";

        break;
    }
    }
}

void scaling_sim::write_moduli() {

   // real_t a = scaleFactor(t);

    double phimod = 0.0;
    double pimod = 0.0;

    hila::set_allreduce(false);
    onsites (ALL) {
        phimod += phi[X].abs();
        pimod += pi[X].abs();
    }

    if (hila::myrank() == 0) {
        config.stream << t << " " << a << " " << sqrt(lambda/2.0) << " "
                      << phimod / lattice->volume() << " " << pimod / lattice->volume()
                      << " ";
    }
}

void scaling_sim::write_energies() {

   // double a = scaleFactor(t);
    double ss = config.sigma * config.sigma;

    // non-weighted energies
    double sumPi = 0.0;    //
    double sumDiPhi = 0.0; //
    double sumV = 0.0;     //
    double sumPhiDiPhi = 0.0;
    double sumPhiPi = 0.0;

    // weighted energies
    double w_sumPi = 0.0;
    double w_sumDiPhi = 0.0;
    double w_sumV = 0.0;
    double w_sumPhiDiPhi = 0.0;
    double w_sumPhiPi = 0.0;

    double phi2 = 0.0;

    hila::synchronize();

    hila::set_allreduce(false);
    onsites (ALL) {
        double phinorm = phi[X].squarenorm();
        double v = 0.25 * lambda * a * a * pow((phinorm - ss), 2.0);
        double pinorm = pi[X].squarenorm();
        double pPi = (phi[X].conj()*pi[X]).re;

        sumV += v;
        w_sumV += v * v;

        sumPi += 0.5 * pinorm;
        w_sumPi += 0.5 * pinorm * v;

        sumPhiPi += 0.5 * pPi * pPi;
        w_sumPhiPi += 0.5 * pPi * pPi * v;

        phi2 += phinorm;
    }


    hila::set_allreduce(false);
    onsites (ALL) {
        auto norm = phi[X].squarenorm();
        real_t v2 = 0.25 * lambda * a * a * pow((norm - ss), 2.0);
        auto diff_phi = (phi[X + e_x] - phi[X - e_x] + phi[X + e_y] - phi[X - e_y] + 
                         phi[X + e_z] - phi[X - e_z]) / (2 * config.dx);
        real_t pDphi = 0.5 * (diff_phi.conj() * phi[X]).re;
        real_t diff_phi_norm2 = diff_phi.squarenorm();

        sumDiPhi += 0.5 * diff_phi_norm2; 
        sumPhiDiPhi += pDphi * pDphi / norm;
        w_sumDiPhi += 0.5 * diff_phi_norm2 * v2; 
        w_sumPhiDiPhi += pDphi * pDphi / norm * v2;
    }
    

    if (hila::myrank() == 0) {
        double vol = (double)config.l * config.l * config.l;
        config.stream << sumPi / vol << " " << w_sumPi / vol << " ";
        config.stream << sumDiPhi / vol << " " << w_sumDiPhi / vol << " ";
        config.stream << sumPhiPi / vol << " " << w_sumPhiPi / vol << " ";
        config.stream << sumPhiDiPhi / vol << " " << w_sumPhiDiPhi / vol << " ";
        config.stream << sumV / vol << " " << w_sumV / vol << " ";
	config.stream << phi2 / vol << " "; 
    }
}

void scaling_sim::write_windings()
{
#ifdef OLD_WINDING

    real_t length = 0.0;
    
    foralldir(d1) foralldir(d2) if (d1 < d2) {
   	onsites(ALL) {
       		float plaq = (phi[X] * phi [X+d1].conj()).arg()  
            	+ (phi[X+d1] * phi[X+d1+d2].conj()).arg()
            	+ (phi[X+d1+d2] * phi[X+d2].conj()).arg()
            	+ (phi[X+d2] * phi[X].conj()).arg();

		length += abs(plaq) * config.dx /(2.0 * M_PI);
	}
    }

    if (hila::myrank() == 0) 
    {
        config.stream << length * config.dx/(2.0 * M_PI) << "\n";
    }

#else

    Reduction<real_t> length(0);
    length.allreduce(false).delayed(true);

    Field<real_t> twist[NDIM];
    foralldir(d) 
        twist[d][ALL] = (phi[X] * phi[X+d].conj()).arg();

    foralldir(d1) foralldir(d2) if (d1 < d2) {
        onsites(ALL) {
            real_t plaq = twist[d1][X] + twist[d2][X+d1] - twist[d1][X+d2] - twist[d2][X];

            length += abs(plaq);            
        }
    }

    auto v = length.value() * config.dx/(2.0 * M_PI);

    if (hila::myrank() == 0) 
    {
        config.stream << v << "\n";
    }


#endif

}

void scaling_sim::next() {

    //real_t a = scaleFactor(t);
    //real_t aHalfPlus = scaleFactor(t + config.dt / 2.0);
    //real_t aHalfMinus = scaleFactor(t - config.dt / 2.0);

    real_t aadt_aadxdx = pow(a / aHalfPlus, 2.0) * config.dt / (config.dx * config.dx);
    real_t aadt2D_aadxdx = aadt_aadxdx * 2.0 * 3.0;
    real_t aaaaldt_aa = pow(a, 4.0) * lambda * config.dt / pow(aHalfPlus, 2.0);
    real_t daa_aa = (pow(aHalfPlus, 2.0) - pow(aHalfMinus, 2.0)) / pow(aHalfPlus, 2.0);
    real_t ss = config.sigma * config.sigma;

    static hila::timer next_timer("timestep");

    Field<Complex<real_t>> deltaPi;

    next_timer.start();

    onsites (ALL) {
        phi[X] += config.dt * pi[X];
        deltaPi[X] = phi[X] * (aaaaldt_aa * (ss - phi[X].squarenorm()) - aadt2D_aadxdx);
    }

    foralldir(d) {
         phi.start_fetch(d); 
         phi.start_fetch(-d);
    }

    foralldir (d) {
          onsites (ALL) {
              deltaPi[X] += aadt_aadxdx * (phi[X + d] + phi[X - d]);
          }
    }

    // onsites (ALL) {
    //     deltaPi[X] += aadt_aadxdx * (phi[X + e_x] + phi[X - e_x] + 
    //                                  phi[X + e_y] + phi[X - e_y] + 
    //                                  phi[X + e_z] + phi[X - e_z]);
    // }

    //pi[ALL] = pi[X] - daa_aa * pi[X] + deltaPi[X];

    if (t < config.tdif)
    {
        pi[ALL] = deltaPi[X]/(config.difFac*config.dt);
        t += config.dt/config.difFac;
    }
    else
    {
        pi[ALL] = pi[X] - daa_aa * pi[X] + deltaPi[X];
        t += config.dt;
    }

    next_timer.stop();

    //t += config.dt;
}

int main(int argc, char **argv) {
    scaling_sim sim;
    const std::string output_fname = sim.allocate("sim_params.txt", argc, argv);
    sim.initialize();

    int steps =
        (sim.config.tEnd - sim.config.tStats) /
        (sim.config.dt * sim.config.nOutputs); // number of steps between printing stats
    if (steps == 0)
        steps = 1;
        
    int stat_counter = 0;

    if (hila::myrank() == 0) {
        sim.config.stream.open(output_fname, std::ios::out);
    }


    // on gpu the simulation timer is fake, because there's no sync here.  
    // BUt we want to avoid unnecessary sync anyway.
    static hila::timer run_timer("Simulation time"), meas_timer("Measurements");
    run_timer.start();
    
    //auto tildephi = sim.phi;
    while (sim.t < sim.config.tEnd) {
	sim.updateCosmology();
        if (sim.t >= sim.config.tStats) {
            if (stat_counter % steps == 0) {
                meas_timer.start();
                sim.write_moduli();
                sim.write_energies();
                sim.write_windings();
		meas_timer.stop();
            }
            stat_counter++;
        }
        sim.next();
    }
    run_timer.stop();

    if (hila::myrank() == 0) {
        sim.config.stream.close();
    }

    hila::finishrun();
    return 0;
}
