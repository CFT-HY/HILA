/*********************************************
 * Simulates staggered fermions with a gauge *
 * interaction                               *
 *********************************************/


#include "staggered.h"





class auxiliary_momentum_action {
  public:

    field<double> &sigma, &pi;
    field<double> &sigma_momentum, &pi_momentum;

    auxiliary_momentum_action(field<double> &s, field<double> &p, field<double> &sm, field<double> &pm) 
    : sigma(s), pi(s), sigma_momentum(sm), pi_momentum(pm){}

    auxiliary_momentum_action(auxiliary_momentum_action &ma)
    : sigma(ma.sigma), pi(ma.pi), sigma_momentum(ma.sigma_momentum), pi_momentum(ma.pi_momentum){}

    double action(){
      double a=0;
      onsites(ALL){
        a += sigma_momentum[X]*sigma_momentum[X]
           + pi_momentum[X]*pi_momentum[X];
      }
      return a;
    }

    /// Gaussian random momentum for each element
    void draw_gaussian_fields(){
      onsites(ALL){
        sigma_momentum[X] = gaussian_ran();
        pi_momentum[X] = gaussian_ran();
      }
    }

    // Integrator step: apply the momentum on the gauge field
    void step(double eps){
      sigma[ALL] = sigma[X] + eps*sigma_momentum[X];
      pi[ALL] = pi[X] + eps*pi_momentum[X];
    }

    // Called by hmc
    void back_up_fields(){}
    void restore_backup(){}
};


class auxiliary_action {
  public:
    field<double> &sigma, &pi;
    field<double> &sigma_momentum, &pi_momentum;
    field<double> sigma_backup, pi_backup;
    double G;

    auxiliary_action(field<double> &s, field<double> &p, field<double> &sm, field<double> &pm, double g)
    : sigma(s), pi(s), sigma_momentum(sm), pi_momentum(pm), G(g){}

    auxiliary_action(auxiliary_action &ma)
    : sigma(ma.sigma), pi(ma.pi), sigma_momentum(ma.sigma_momentum), pi_momentum(ma.pi_momentum), G(ma.G){}

    double action(){
      double a=0;
      onsites(ALL){
        a += 1.0/(4.0*G*G) * (sigma[X]*sigma[X] + pi[X]*pi[X]);
      }
      return a;
    }

    void draw_gaussian_fields(){}

    // Update the momentum with the auxiliary field
    void force_step(double eps){
      sigma_momentum[ALL] = sigma_momentum[X] - 2*eps*1.0/(4.0*G*G)*sigma[X];
      pi_momentum[ALL] = pi_momentum[X] - 2*eps*1.0/(4.0*G*G)*pi[X];
    }

    // Make a copy of fields updated in a trajectory
    void back_up_fields(){
      sigma_backup = sigma;
      pi_backup = pi;
    }

    // Restore the previous backup
    void restore_backup(){
      sigma = sigma_backup;
      pi = pi_backup;
    }
};




int main(int argc, char **argv){

  input parameters = input();
  parameters.import("parameters");
  double beta = parameters.get("beta");
  double gamma = parameters.get("gamma");
  double mass = parameters.get("mass");
  int seed = parameters.get("seed");
	double hmc_steps = parameters.get("hmc_steps");
	double traj_length = parameters.get("traj_length");

  lattice->setup( nd[0], nd[1], nd[2], nd[3], argc, argv );
  seed_random(seed);


  // Define gauge field and momentum field
  field<SUN> gauge[NDIM];
  field<SUN> momentum[NDIM];

  // Define gauge and momentum action terms
  gauge_momentum_action gma(gauge, momentum);
  gauge_action ga(gauge, momentum, beta);

  // Initialize the gauge field
  ga.set_unity();
  
  // Define auxiliary fields and momentum
  field<double> sigma, pi, sigma_mom, pi_mom;

  // And action
  auxiliary_momentum_action ama(sigma, pi, sigma_mom, pi_mom);
  auxiliary_action aa(sigma, pi, sigma_mom, pi_mom, gamma);

  // Initialize
  sigma[ALL] = 0;
  pi[ALL] = 0;

  // Define a Dirac operator (2 flavors)
  dirac_staggered_NJL<N> D(mass, gauge, sigma, pi);
  fermion_action fa(D, gauge, momentum);


  // Build two integrator levels. Gauge is on the lowest level and
  // the fermions are on higher level
  integrator integrator_level_1(ga+aa, gma+ama);
  integrator integrator_level_2(fa, integrator_level_1);
  

  // Run HMC using the integrator
  for(int step = 0; step < 100; step ++){
    // Run update
    update_hmc(integrator_level_2, hmc_steps, traj_length);
    
    // Measurements

    double plaq = plaquette(ga.gauge);
    output0 << "Plaq: " << plaq << "\n";

    double sigmasq = 0, sigma_ave = 0;
    onsites(ALL){
      sigma_ave += sigma[X];
      sigmasq += sigma[X]*sigma[X];
    }
    output0 << "Sigma: " << sigma_ave << "\n";
    output0 << "Sigma sq: " << sigmasq << "\n";
  }



  finishrun();

  return 0;
}
