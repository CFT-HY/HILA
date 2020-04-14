/*********************************************
 * Simulates staggered fermions with a gauge *
 * interaction                               *
 *********************************************/


#include "staggered.h"
#include "../plumbing/inputs.h"
#include "../plumbing/algorithms/hmc.h"





/// Gaussian random momentum for each element
void project_antihermitean(element<NMAT> &matrix){
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



/// Calculate the sum of staples connected to links in direction dir 
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



/// Measure the plaquette
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

double plaquette(field<SUN> *gauge){
  return plaquette_sum(gauge)/(lattice->volume()*NDIM*(NDIM-1));
}





// The momentum action
double momentum_action(field<NMAT> *momentum){
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

void gaussian_momentum_impl(field<NMAT> *momentum){
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


void force_step_impl(field<SUN> *gauge, field<NMAT> *momentum, double eps){
  foralldir(dir){
    field<SUN> staples = calc_staples(gauge, dir);
    onsites(ALL){
      element<NMAT> force;
      force = gauge[dir][X]*staples[X];
      project_antihermitean(force);
      momentum[dir][X] = momentum[dir][X] - eps*force;
    }
  }
}



class gauge_term{
  public:
    field<SUN> *gauge;
    field<NMAT> *momentum;
    double beta;

    gauge_term(field<SUN> *g, field<NMAT> *m, double b){
      gauge = g; momentum = m; beta = b;
    }

    //The gauge action
    double action(){
      return beta*plaquette_sum(gauge) + momentum_action(momentum);
    }

    /// Gaussian random momentum for each element
    void gaussian_momentum(){
      gaussian_momentum_impl(momentum);
    }

    // Update the momentum with the gauge field
    void force_step(double eps){
      foralldir(dir){
        field<SUN> staples = calc_staples(gauge, dir);
        onsites(ALL){
          element<NMAT> force;
          force = gauge[dir][X]*staples[X];
          project_antihermitean(force);
          momentum[dir][X] = momentum[dir][X] - beta*eps/N*force;
        }
      }
    }

    // Update the gauge field with momentum
    void momentum_step(double eps){
      foralldir(dir){
        onsites(ALL){
          element<SUN> momexp = eps*momentum[dir][X];
          momexp.exp();
          gauge[dir][X] = momexp*gauge[dir][X];
        }
      }
    }

    // A single gauge update
    void integrator_step(double eps){
      O2_step(*this, eps);
    }
};



class fermion_term{
  public:
    gauge_term gt;
    field<SUN> *gauge;
    field<NMAT> *momentum;

    fermion_term(gauge_term g, field<NMAT> *m) : gt(g) {
      gauge = g.gauge; momentum = m;
    }

    //The gauge action
    double action(){
      return gt.action();
    }

    /// Gaussian random momentum for each element
    void gaussian_momentum(){
      gt.gaussian_momentum();
    }

    // Update the momentum with the gauge field
    void force_step(double eps){
      // Fermion force here
     }

    // Update the gauge field with momentum
    void momentum_step(double eps){
      gt.integrator_step(eps);
    }

    // A single gauge update
    void integrator_step(double eps){
      O2_step(*this, eps);
    }
};


class full_action{
  public:
    fermion_term ft;
    field<SUN> *gauge;

    full_action(fermion_term f) : ft(f) {gauge = f.gauge;}

    //The gauge action
    double action(){
      return ft.action();
    }

    /// Gaussian random momentum for each element
    void gaussian_momentum(){
      ft.gaussian_momentum();
    }

    // Update the momentum with the gauge field
    void integrate(int steps, double dt){
      for(int step=0; step < steps; step++){
        ft.integrator_step(dt/steps);
      }
    }
};










int main(int argc, char **argv){

  input parameters = input();
  parameters.import("parameters");
  double beta = parameters.get("beta");
  int seed = parameters.get("seed");
	double hmc_steps = parameters.get("hmc_steps");
	double traj_length = parameters.get("traj_length");

  lattice->setup( nd[0], nd[1], nd[2], nd[3], argc, argv );
  seed_random(seed);


  field<SUN> gauge[NDIM];

  // Starting with a unit configuration
  foralldir(dir){
    onsites(ALL){
      gauge[dir][X].random();
    }
  }


  field<NMAT> momentum[NDIM];
  foralldir(dir){
    onsites(ALL){
      momentum[dir][X] = 0;
    }
  }


  // Test the force calculation by varying one gauge link
  // (Needs to be moved to tests)
  double eps = 1e-6;
  SUN g1 = gauge[0].get_value_at(50);
  SUN h = 1;
  h.c[0][1].re += eps;
  h.c[1][0].re -= eps;
  SUN g12 = h*g1;

  gauge_term gt(gauge, momentum, 1.0);

  if(mynode()==0)
    gauge[0].set_value_at(g1, 50);
  gauge[0].mark_changed(ALL);
  double s1 = plaquette_sum(gauge);

  if(mynode()==0)
    gauge[0].set_value_at(g12,50);
  gauge[0].mark_changed(ALL);
  double s2 = plaquette_sum(gauge);

  if(mynode()==0)
    gauge[0].set_value_at(g1, 50);
  gauge[0].mark_changed(ALL);

  force_step_impl(gauge, momentum, 1.0/N);
  NMAT f = momentum[0].get_value_at(50);
  double diff = 2*f.c[0][1].re + (s2-s1)/eps;
  if(mynode()==0) 
    assert( diff*diff < eps*eps*100 );


  // Check also the momentum action and derivative
  gt.gaussian_momentum();

  s1 = momentum_action(momentum);
  h = momentum[0].get_value_at(0);
  h.c[0][0].im += eps;
  if(mynode()==0)
    momentum[0].set_value_at(h, 0);
  s2 = momentum_action(momentum);

  diff = h.c[0][0].im - (s2-s1)/eps;
  if(mynode()==0) 
    assert( diff*diff < eps*eps*100 );




  // Now the actual simulation
  gt = gauge_term(gauge, momentum, beta);
  fermion_term ft = fermion_term(gt, momentum);
  full_action action = full_action(ft);
  for(int step = 0; step < 100; step ++){
    update_hmc(action, hmc_steps, traj_length);
    double plaq = plaquette(gauge);
    output0 << "Plaq: " << plaq << "\n";
  }



  finishrun();

  return 0;
}
