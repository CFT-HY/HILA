#ifndef GAUGE_FIELD_H
#define GAUGE_FIELD_H


// A conveniance class for a gauge field.
// Contains an SU(N) matrix in each direction for the
// gauge field and for the momentum
template<int N,typename radix=double>
struct gauge_field {
  using gauge_type = SU<N,radix>;
  field<SU<N,double>> gauge[NDIM];
  field<SU<N,double>> momentum[NDIM];

  // Set the gauge field to unity
  void set_unity(){
    foralldir(dir){
      onsites(ALL){
        gauge[dir][X] = 1;
      }
    }
  }

  // Read the gauge field from a file
  void read_file(std::string filename){
    std::ifstream inputfile;
    inputfile.open(filename, std::ios::in | std::ios::binary);
    foralldir(dir){
      read_fields(inputfile, gauge[dir]);
    }
    inputfile.close();
  }

  // Write the gauge field to a file
  void write_file(std::string filename){
    std::ofstream outputfile;
    outputfile.open(filename, std::ios::out | std::ios::trunc | std::ios::binary);
    foralldir(dir){
      write_fields(outputfile, gauge[dir]);
    }
    outputfile.close();
  }
};




/// Calculate the Polyakov loop for a given gauge field.
template<int N>
double polyakov_loop(direction dir, field<SU<N>> (&gauge)[NDIM]){
  // This implementation uses the onsites() to cycle through the
  // NDIM-1 dimensional planes. This is probably not the most
  // efficient implementation.
  coordinate_vector vol = lattice->size();
  field<SU<N>> polyakov; polyakov[ALL] = 1;
  for(int t=0; t<vol[dir]; t++){
    onsites(ALL){
      if(X.coordinates()[dir]==(t+1)%vol[dir]){
        polyakov[X] = polyakov[X] * gauge[dir][X-dir];
      }
    }
  }

  double poly=0;
  onsites(ALL){
    if(X.coordinates()[dir]==0){
      poly += polyakov[X].trace().re;
    }
  }

  double v3 = lattice->volume()/vol[dir];
  return poly/(N*v3);
}





/// Calculate the sum of staples connected to links in direction dir 
template<typename SUN>
field<SUN> calc_staples(field<SUN> *U, direction dir)
{
  field<SUN> staple_sum;
  static field<SUN> down_staple;
  staple_sum[ALL] = 0;
  foralldir(dir2){
    //Calculate the down side staple.
    //This will be communicated up.
    down_staple[ALL] = U[dir2][X+dir].conjugate()
                     * U[dir][X].conjugate()
                     * U[dir2][X];
    // Forward staple
    staple_sum[ALL]  = staple_sum[X]
                     + U[dir2][X+dir]
                     * U[dir][X+dir2].conjugate()
                     * U[dir2][X].conjugate();
    // Add the down staple
    staple_sum[ALL] = staple_sum[X] + down_staple[X-dir2];
  }
  return staple_sum;
}


/// Apply the force of the gauge field on the momentum field 
template<typename SUN, typename MATRIX>
void gauge_force(field<SUN> *gauge, field<MATRIX> *momentum, double eps){
  foralldir(dir){
    field<SUN> staples = calc_staples(gauge, dir);
    onsites(ALL){
      element<MATRIX> force;
      force = gauge[dir][X]*staples[X];
      project_antihermitean(force);
      momentum[dir][X] = momentum[dir][X] - eps*force;
    }
  }
}




/// Measure the plaquette
template<int N, typename radix>
double plaquette_sum(field<SU<N,radix>> *U){
  double Plaq=0;
  foralldir(dir1) foralldir(dir2) if(dir2 < dir1){
    onsites(ALL){
      element<SU<N,radix>> temp;
      temp = U[dir1][X] * U[dir2][X+dir1]
           * U[dir1][X+dir2].conjugate()
           * U[dir2][X].conjugate();
      Plaq += 1-temp.trace().re/N;
    }
  }
  return Plaq;
}


template<typename SUN>
double plaquette(field<SUN> *gauge){
  return plaquette_sum(gauge)/(lattice->volume()*NDIM*(NDIM-1)/2);
}






// Action term for the momentum of a gauge field
// This is both an action term and an integrator. It can form the
// lowest level step to an integrator construct.
template<typename gauge_field>
class gauge_momentum_action {
  public:
    using gauge_mat = typename gauge_field::gauge_type;

    gauge_field &gauge;
    field<gauge_mat> gauge_copy[NDIM];

    gauge_momentum_action(gauge_field &g) 
    : gauge(g){}

    gauge_momentum_action(gauge_momentum_action &ma)
    : gauge(ma.gauge){}

    double action(){
      double sum = 0;
      foralldir(dir) {
        onsites(ALL){
          sum += gauge.momentum[dir][X].algebra_norm();
        }
      }
      return sum;
    }

    /// Gaussian random momentum for each element
    void draw_gaussian_fields(){
      foralldir(dir) {
        onsites(ALL){
          if(disable_avx[X]==0){};
          gauge.momentum[dir][X].gaussian_algebra();
        }
      }
    }

    // Integrator step: apply the momentum on the gauge field
    void step(double eps){
      foralldir(dir){
        onsites(ALL){
          element<gauge_mat> momexp = eps*gauge.momentum[dir][X];
          momexp.exp();
          gauge.gauge[dir][X] = momexp*gauge.gauge[dir][X];
        }
      }
    }

    // Called by hmc
    // Make a copy of fields updated in a trajectory
    void back_up_fields(){
      foralldir(dir) gauge_copy[dir] = gauge.gauge[dir];
    }

    // Restore the previous backup
    void restore_backup(){
      foralldir(dir) gauge.gauge[dir] = gauge_copy[dir];
    }
};



// Represents a sum of two momentum terms. Useful for adding them
// to an integrator on the same level.
template<typename momentum_action_1, typename momentum_action_2>
class momentum_action_sum {
  public:
    momentum_action_1 a1;
    momentum_action_2 a2;

    momentum_action_sum(momentum_action_1 _a1, momentum_action_2 _a2) 
    : a1(_a1), a2(_a2){}

    momentum_action_sum(momentum_action_sum &asum) : a1(asum.a1), a2(asum.a2){}

    //The gauge action
    double action(){
      return a1.action() + a2.action();
    }

    /// Gaussian random momentum for each element
    void draw_gaussian_fields(){
      a1.draw_gaussian_fields();
      a2.draw_gaussian_fields();
    }

    // Called by hmc
    void back_up_fields(){}
    void restore_backup(){}

    // Integrator step: apply the momentum on the gauge fields
    void step(double eps){
      a1.step(eps);
      a2.step(eps);
    }

};

// Sum operator for creating a momentum action_sum object
template<typename gauge_type, typename action2>
momentum_action_sum<gauge_momentum_action<gauge_type>, action2> operator+(gauge_momentum_action<gauge_type> a1, action2 a2){
  momentum_action_sum<gauge_momentum_action<gauge_type>, action2> sum(a1, a2);
  return sum;
}


template<typename gauge_field>
class gauge_action {
  public:
    using gauge_mat = typename gauge_field::gauge_type;
    static constexpr int N = gauge_mat::size;

    gauge_field &gauge;
    double beta;

    gauge_action(gauge_field &g, double b) 
    : gauge(g), beta(b){}

    gauge_action(gauge_action &ga)
    : gauge(ga.gauge), beta(ga.beta) {}

    //The gauge action
    double action(){
      double s = beta*plaquette_sum(gauge.gauge); 
      return s;
    }

    /// Gaussian random momentum for each element
    void draw_gaussian_fields(){}

    // Update the momentum with the gauge field
    void force_step(double eps){
      gauge_force(gauge.gauge, gauge.momentum, beta*eps/N);
    }

    // Draw a random gauge field
    void random(){
      foralldir(dir){
        onsites(ALL){
          gauge.gauge[dir][X].random();
        }
      }
    }


    // Called by HMC. The gauge field is copied by 
    // the momentum action.
    void back_up_fields(){}
    void restore_backup(){}
};


// Represents a sum of two acttion terms. Useful for adding them
// to an integrator on the same level.
template<typename action_type_1, typename action_type_2>
class action_sum {
  public:
    action_type_1 a1;
    action_type_2 a2;

    action_sum(action_type_1 _a1, action_type_2 _a2) 
    : a1(_a1), a2(_a2){}

    action_sum(action_sum &asum) : a1(asum.a1), a2(asum.a2){}

    //The gauge action
    double action(){
      return a1.action() + a2.action();
    }

    /// Gaussian random momentum for each element
    void draw_gaussian_fields(){
      a1.draw_gaussian_fields();
      a2.draw_gaussian_fields();
    }

    // Update the momentum with the gauge field
    void force_step(double eps){
      a1.force_step(eps);
      a2.force_step(eps);
    }

    // Set the gauge field to unity
    void set_unity(){
      a1.set_unity();
      a2.set_unity();
    }

    // Draw a random gauge field
    void random(){
      a1.random();
      a2.random();
    }


    // Make a copy of fields updated in a trajectory
    void back_up_fields(){
      a1.back_up_fields();
      a2.back_up_fields();
    }

    // Restore the previous backup
    void restore_backup(){
      a1.restore_backup();
      a2.restore_backup();
    }
};


// Sum operator for creating an action_sum object
template<typename gauge_field, typename action2>
action_sum<gauge_action<gauge_field>, action2> operator+(gauge_action<gauge_field> a1, action2 a2){
  action_sum<gauge_action<gauge_field>, action2> sum(a1, a2);
  return sum;
}



/// Define an integration step for a Molecular Dynamics
/// trajectory.
template<typename action_type, typename lower_integrator_type>
class integrator{
  public:
    action_type action_term;
    lower_integrator_type lower_integrator;

    integrator(action_type a, lower_integrator_type i)
    : action_term(a), lower_integrator(i) {}

    // The current total action of fields updated by this
    // integrator. This is kept constant up to order eps^3.
    double action(){
      return action_term.action() + lower_integrator.action();
    }

    // Refresh fields that can be drawn from a gaussian distribution
    // This is needed at the beginning of a trajectory
    void draw_gaussian_fields(){
      action_term.draw_gaussian_fields();
      lower_integrator.draw_gaussian_fields();
    }

    // Make a copy of fields updated in a trajectory
    void back_up_fields(){
      action_term.back_up_fields();
      lower_integrator.back_up_fields();
    }

    // Restore the previous backup
    void restore_backup(){
      action_term.restore_backup();
      lower_integrator.restore_backup();
    }


    // Update the momentum with the gauge field
    void force_step(double eps){
      action_term.force_step(eps);
    }

    // Update the gauge field with momentum
    void momentum_step(double eps){
      lower_integrator.step(eps);
    }

    // A single gauge update
    void step(double eps){
      O2_step(*this, eps);
    }

};







#endif