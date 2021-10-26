#define _USE_MATH_DEFINES
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
//#include <math.h>
#include <assert.h>

#include "plumbing/hila.h"

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
    void next();
    inline real_t scaleFactor(real_t t);

    Field<Complex<real_t>> phi;
    Field<Complex<real_t>> pi;
    Field<Complex<real_t>> deltaPi;

    real_t t;

    struct config {
        int l;
        int m;
        int seed;
        int smoothing;
        int initalCondition;
        real_t initialModulus;
        real_t epsilon;
        real_t sigma;
        real_t dx;
        real_t dt;
        real_t tStart;
        real_t tStats;
        real_t nOutputs;
        real_t tEnd;
        real_t lambda;
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
    config.tStart = parameters.get("tStart");
    config.tEnd = parameters.get("tEnd");
    config.lambda = parameters.get("lambda");
    config.smoothing = parameters.get("smooth");
    config.initalCondition = parameters.get("initialCondition");
    config.tStats = parameters.get("tStats");
    config.nOutputs = parameters.get("numberStatsOutputs");
    real_t ratio = parameters.get("dtdxRatio");
    const std::string output_file = parameters.get("output_file");
    config.dt = config.dx * ratio;
    t = config.tStart;

    CoordinateVector box_dimensions = {config.l, config.l, config.l};
    lattice->setup(box_dimensions);
    hila::seed_random(config.seed);

    return output_file;
}

inline real_t scaling_sim::scaleFactor(real_t t) {
    return t / config.tEnd;
}

void scaling_sim::initialize() {

    int m = config.m;
    real_t epsilon = config.epsilon;
    real_t s = config.sigma;
    int N = config.l;

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

    real_t a = scaleFactor(t);

    double phimod = 0.0;
    double pimod = 0.0;

    onsites (ALL) {
        phimod += phi[X].abs();
        pimod += pi[X].abs();
    }

    if (hila::myrank() == 0) {
        config.stream << t << "," << a << "," << config.lambda << ","
                      << phimod / lattice->volume() << "," << pimod / lattice->volume()
                      << ",";
    }
}

void scaling_sim::write_energies() {

    double a = scaleFactor(t);
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

    onsites (ALL) {
        double phinorm = phi[X].squarenorm();
        double v = 0.25 * config.lambda * a * a * pow((phinorm - ss), 2.0);
        double pinorm = pi[X].squarenorm();
        double pPi = (phi[X].conj()*pi[X]).re;

        sumV += v;
        w_sumV += v * v;

        sumPi += 0.5 * pinorm;
        w_sumPi += 0.5 * pinorm * v;

        sumPhiPi += 0.5 * pPi * pPi;
        w_sumPhiPi += 0.5 * pPi * pPi * v;
    }

    static hila::timer vreduction("reduction");

    vreduction.start();

    onsites (ALL) {
            auto norm = phi[X].squarenorm();
            real_t v = 0.25 * config.lambda * a * a * pow((norm - ss), 2.0);
            auto diff_phi = (phi[X + e_x] - phi[X - e_x] + phi[X + e_y] - phi[X - e_y] + 
                             phi[X + e_z] - phi[X - e_z]) / (2 * config.dx);
            real_t pDphi = 0.5 * (diff_phi.conj() * phi[X]).re;
            real_t diff_phi_norm2 = diff_phi.squarenorm();

            sumDiPhi += 0.5 * diff_phi_norm2; 
            // red[0] += 0.5 * diff_phi_norm2;
            sumPhiDiPhi += pDphi * pDphi / norm;
            // red[1] += pDphi * pDphi / norm;
            w_sumDiPhi += 0.5 * diff_phi_norm2 * v; 
            // red[2] += 0.5 * diff_phi_norm2 * v;
            w_sumPhiDiPhi += pDphi * pDphi / norm * v;
            // red[3] += pDphi * pDphi / norm * v;

        
    }

    // sumDiPhi = red[0];
    // sumPhiDiPhi = red[1];
    // w_sumDiPhi = red[2];
    // w_sumPhiDiPhi = red[3];
    
    vreduction.stop();

    if (hila::myrank() == 0) {
        double vol = (double)config.l * config.l * config.l;
        config.stream << sumPi / vol << "," << w_sumPi / vol << ",";
        config.stream << sumDiPhi / vol << "," << w_sumDiPhi / vol << ",";
        config.stream << sumPhiPi / vol << "," << w_sumPhiPi / vol << ",";
        config.stream << sumPhiDiPhi / vol << "," << w_sumPhiDiPhi / vol << ",";
        config.stream << sumV / vol << "," << w_sumV / vol << "\n";
    }
}

void scaling_sim::next() {

    real_t a = scaleFactor(t);
    real_t aHalfPlus = scaleFactor(t + config.dt / 2.0);
    real_t aHalfMinus = scaleFactor(t - config.dt / 2.0);

    real_t aadt_aadxdx = pow(a / aHalfPlus, 2.0) * config.dt / (config.dx * config.dx);
    real_t aadt2D_aadxdx = aadt_aadxdx * 2.0 * 3.0;
    real_t aaaaldt_aa = pow(a, 4.0) * config.lambda * config.dt / pow(aHalfPlus, 2.0);
    real_t daa_aa = (pow(aHalfPlus, 2.0) - pow(aHalfMinus, 2.0)) / pow(aHalfPlus, 2.0);
    real_t ss = config.sigma * config.sigma;

    static hila::timer next_timer("timestep");

    next_timer.start();

    foralldir (d) {
              phi.start_fetch(d);
              phi.start_fetch(-d);
    }

    onsites (ALL) {
        phi[X] += config.dt * pi[X];
        deltaPi[X] = phi[X] * (aaaaldt_aa * (ss - phi[X].squarenorm()) - aadt2D_aadxdx);
    }

    // foralldir (d) {
    //      onsites (ALL) {
    //          deltaPi[X] += aadt_aadxdx * (phi[X + d] + phi[X - d]);
    //      }
    // }

    onsites (ALL) {
        deltaPi[X] += aadt_aadxdx * (phi[X + e_x] + phi[X - e_x] + 
                                     phi[X + e_y] + phi[X - e_y] + 
                                     phi[X + e_z] + phi[X - e_z]);
    }

    pi[ALL] = pi[X] - daa_aa * pi[X] + deltaPi[X];

    next_timer.stop();

    t += config.dt;
}

int main(int argc, char **argv) {
    scaling_sim sim;
    const std::string output_fname = sim.allocate("sim_params.txt", argc, argv);
    sim.initialize();

    int steps =
        (sim.config.tEnd - sim.config.tStats) /
        (sim.config.dt * sim.config.nOutputs); // number of steps between printing stats
    int stat_counter = 0;

    if (hila::myrank() == 0) {
        sim.config.stream.open(output_fname, std::ios::out);
    }


    // on gpu the simulation timer is fake, because there's no synch here.  
    // BUt we want to avoid unnecessary sync anyway.
    static hila::timer run_timer("Simulation time"), meas_timer("Measurements");
    run_timer.start();
    while (sim.t < sim.config.tEnd) {
        if (sim.t >= sim.config.tStats) {
            if (stat_counter % steps == 0) {
                meas_timer.start();
                sim.write_moduli();
                sim.write_energies();
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
