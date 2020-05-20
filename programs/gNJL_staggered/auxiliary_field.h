

class hypercube{
  private:
    int s[16];
  public:
  int size = 16;
  int & operator[](const int i) {
    return s[i];
  }
};

std::vector<hypercube> get_hypercubes(){
  static std::vector<hypercube> hypercubes;
  const int loop_begin = lattice->loop_begin(parity::all);
  const int loop_end   = lattice->loop_end(parity::all);

  if(hypercubes.size() == 0){
    for(int Index = loop_begin; Index < loop_end; ++Index){
      coordinate_vector l = lattice->coordinates(Index);
      bool is_hypercube_root = true;
      foralldir(d){
        is_hypercube_root *= (l[d]%2) == 0;
      }
      if(is_hypercube_root){
        hypercube h;
        int hi = 0;
        h[hi] = Index;
        foralldir(d1){
          int id1 = lattice->neighb[d1][Index];
          hi++; h[hi] = id1;
          foralldir(d2) if(d2>d1) {
            int id2 = lattice->neighb[d2][id1];
            hi++; h[hi] = id2;
            foralldir(d3) if(d3>d2) {
              int id3 = lattice->neighb[d3][id2];
              hi++; h[hi] = id3;
              foralldir(d4) if(d4>d3) {
                int id4 = lattice->neighb[d4][id3];
                hi++; h[hi] = id4;
              }
            }
          }
        }
        hypercubes.push_back(h);
      }
    }
  }

  return hypercubes;
}


void hypercube_sum(const field<double> &f1, field<double> &f2){

  f1.check_alloc();
  f2.check_alloc();
  for(hypercube cube : get_hypercubes()){
    double val = 0;
    for(int ci=0; ci < cube.size; ci++){
      val += f1.get_value_at(cube[ci]);
    }
    for(int ci=0; ci < cube.size; ci++){
      f2.set_value_at(val, cube[ci]);
    }
  }
  f2.mark_changed(ALL);
}


void momentum_gaussian_hypercube(field<double> &mom){

  mom.check_alloc();
  for(hypercube cube : get_hypercubes()){
    double val = gaussian_ran();
    mom.set_value_at(val, cube[0]);
    for(int ci=1; ci < cube.size; ci++){
      mom.set_value_at(0, cube[ci]);
    }
  }
  mom.mark_changed(ALL);
}


class auxiliary_momentum_action {
  public:

    field<double> &sigma, &pi;
    field<double> &sigma_momentum, &pi_momentum;
    double alpha=1;

    auxiliary_momentum_action(field<double> &s, field<double> &p, field<double> &sm, field<double> &pm, double a=1) 
    : sigma(s), pi(p), sigma_momentum(sm), pi_momentum(pm), alpha(a){}

    auxiliary_momentum_action(auxiliary_momentum_action &ma)
    : sigma(ma.sigma), pi(ma.pi), sigma_momentum(ma.sigma_momentum), pi_momentum(ma.pi_momentum), alpha(ma.alpha){}

    double action(){
      field<double> tpi, tsigma;
      hypercube_sum(sigma_momentum, tsigma);
      hypercube_sum(pi_momentum, tpi);
      double a=0;
      onsites(ALL){
        a += alpha*0.5*(tsigma[X]*tsigma[X] + tpi[X]*tpi[X])/16/16;
      }
      return a;
    }

    /// Gaussian random momentum for each element
    void draw_gaussian_fields(){
      momentum_gaussian_hypercube(sigma_momentum);
      momentum_gaussian_hypercube(pi_momentum);
    }

    // Integrator step: apply the momentum on the gauge field
    void step(double eps){
      field<double> tpi, tsigma;
      hypercube_sum(sigma_momentum, tsigma);
      hypercube_sum(pi_momentum, tpi);
      sigma[ALL] = sigma[X] - alpha*eps*tsigma[X]/16;
      pi[ALL] = pi[X] - alpha*eps*tpi[X]/16;
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
        a += 1.0/(4*G*G) * (sigma[X]*sigma[X] + pi[X]*pi[X]);
      }
      return a;
    }

    void draw_gaussian_fields(){}

    // Update the momentum with the auxiliary field
    void force_step(double eps){
      sigma_momentum[ALL] = sigma_momentum[X] + 2*eps/(4.0*G*G)*sigma[X];
      pi_momentum[ALL] = pi_momentum[X] + 2*eps/(4.0*G*G)*pi[X];
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


