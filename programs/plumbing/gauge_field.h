#ifndef GAUGE_FIELD_H
#define GAUGE_FIELD_H


/// Project to the antihermitean part of a matrix
template<typename MATRIX>
void project_antihermitean(MATRIX &matrix){
  double tr = 0;
  for(int i=0; i<N; i++) {
    for(int j=0; j<i; j++) {
      double a = 0.5*(matrix.c[i][j].re - matrix.c[j][i].re);
      double b = 0.5*(matrix.c[i][j].im + matrix.c[j][i].im);
      matrix.c[i][j].re = a;
      matrix.c[j][i].re =-a;
      matrix.c[i][j].im = b;
      matrix.c[j][i].im = b;
    }
    tr += matrix.c[i][i].im;
    matrix.c[i][i].re = 0;
  }
  for(int i=0; i<N; i++) {
    matrix.c[i][i].im -= tr/N;
  }
}

/// Generate a random antihermitean matrix
template<typename MATRIX>
void gaussian_momentum(field<MATRIX> *momentum){
  foralldir(dir) {
    onsites(ALL){
      for(int i=0; i<N; i++) {
        for(int j=0; j<i; j++) {
          double a = gaussian_ran();
          double b = gaussian_ran();
          momentum[dir][X].c[i][j].re = a;
          momentum[dir][X].c[j][i].re =-a;
          momentum[dir][X].c[i][j].im = b;
          momentum[dir][X].c[j][i].im = b;
        }
      }

      for(int i=0; i<N; i++) {
        momentum[dir][X].c[i][i].re = 0;
        momentum[dir][X].c[i][i].im = 0;
      }
      for(int i=1; i<N; i++) {
        double a = gaussian_ran()*sqrt(2.0/(i*(i+1)));
        for(int j=0; j<i; j++)
          momentum[dir][X].c[j][j].im += a;
        momentum[dir][X].c[i][i].im -= i*a;
      }
    }
  }
}


// The momentum action
template<typename MATRIX>
double momentum_action(field<MATRIX> *momentum){
  double sum = 0;
  foralldir(dir) {
    onsites(ALL){
      double thissum = 0;
      for(int i=0; i<N; i++) {
        for(int j=0; j<i; j++) {
          thissum += momentum[dir][X].c[i][j].squarenorm();
        }
        double diag = momentum[dir][X].c[i][i].im;
        thissum += 0.5*diag*diag;
      }
      sum += thissum;
    }
  }
  return sum;
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


/// Apply the momentum on the gauge field
template<typename SUN, typename MATRIX>
void apply_momentum(field<SUN> *gauge, field<MATRIX> *momentum, double eps){
  foralldir(dir){
    onsites(ALL){
      element<SUN> momexp = eps*momentum[dir][X];
      momexp.exp();
      gauge[dir][X] = momexp*gauge[dir][X];
    }
  }
}




/// Measure the plaquette
template<typename SUN>
double plaquette_sum(field<SUN> *U){
  double Plaq=0;
  foralldir(dir1) foralldir(dir2) if(dir2 < dir1){
    onsites(ALL){
      element<SUN> temp;
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
  return plaquette_sum(gauge)/(lattice->volume()*NDIM*(NDIM-1));
}




template<int N, typename float_t>
class gauge_action{
  public:
    using SUN = SU<N, float_t>;
    using MATRIX = matrix<N,N,cmplx<float_t>>;

    field<SUN> *gauge;
    field<MATRIX> *momentum;
    double beta;

    gauge_action(field<SUN> *g, field<MATRIX> *m, double b){
      gauge = g; momentum = m; beta = b;
    }

    //The gauge action
    double action(){
      return beta*plaquette_sum(gauge) + momentum_action(momentum);
    }

    /// Gaussian random momentum for each element
    void generate_momentum(){
      gaussian_momentum(momentum);
    }

    // Update the momentum with the gauge field
    void force_step(double eps){
      gauge_force(gauge, momentum, beta*eps/N);
    }

    // Update the gauge field with momentum
    void momentum_step(double eps){
      apply_momentum(gauge, momentum, eps);
    }

    // A single gauge update
    void integrator_step(double eps){
      O2_step(*this, eps);
    }
};





#endif