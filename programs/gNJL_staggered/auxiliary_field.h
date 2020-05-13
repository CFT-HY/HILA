
class auxiliary_momentum_action {
  public:

    field<double> &sigma, &pi;
    field<double> &sigma_momentum, &pi_momentum;

    auxiliary_momentum_action(field<double> &s, field<double> &p, field<double> &sm, field<double> &pm) 
    : sigma(s), pi(p), sigma_momentum(sm), pi_momentum(pm){}

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
    : sigma(s), pi(p), sigma_momentum(sm), pi_momentum(pm), G(g){}

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
      sigma_momentum[ALL] = sigma_momentum[X] - eps*1.0/(4.0*G*G)*sigma[X];
      pi_momentum[ALL] = pi_momentum[X] - eps*1.0/(4.0*G*G)*pi[X];
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

