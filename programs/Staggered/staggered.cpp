/*********************************************
 * Simulates staggered fermions with a gauge *
 * interaction                               *
 *********************************************/


#include "staggered.h"
#include "../plumbing/inputs.h"






/// Gaussian random momentum for each element
void gaussian_momentum(field<NMAT> *momentum){
  foralldir(dir) {
    onsites(ALL){
      for(int i=0; i<N; i++) {
        for(int j=0; j<i; j++) {
          double a = gaussian_ran();
          double b = gaussian_ran();
          momentum[dir][X].c[i][j].re = a;
          momentum[dir][X].c[i][j].im = b;
          momentum[dir][X].c[j][i].re =-a;
          momentum[dir][X].c[j][i].im = b;
        }
      }

      for(int i=1; i<N; i++) {
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


/// Gaussian random momentum for each element
void project_antihermitean(element<NMAT> &matrix){
  double tr = 0;
  for(int i=0; i<N; i++) {
    for(int j=0; j<N; j++) {
      double a = 0.5*(matrix.c[i][j].re - matrix.c[j][i].re);
      double b = 0.5*(matrix.c[i][j].im + matrix.c[j][i].im);
      matrix.c[i][j].re = a;
      matrix.c[i][j].im = b;
      matrix.c[j][i].re =-a;
      matrix.c[j][i].im = b;
    }
    tr += matrix.c[i][i].im;
    matrix.c[i][i].re = 0;
  }
  for(int i=0; i<N; i++) {
    matrix.c[i][i].im -= tr/N;
  }
}


/// Calculate the action of the gauge momentum
double momentum_action(field<NMAT> *momentum){
  double sum = 0;
  foralldir(dir) {
    onsites(ALL){
      double thissum = 0;
      for(int i=0; i<N; i++) {
        for(int j=0; j<N; j++) {
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
double plaquette(field<SUN> *U){
  double Plaq=0;
  foralldir(dir1) foralldir(dir2) if(dir1 != dir2){
    onsites(ALL){
      element<SUN> temp;
      temp =  U[dir1][X] * U[dir2][X+dir1];
      temp = temp*U[dir1][X+dir2].conjugate();
      temp = temp*U[dir2][X].conjugate();
      Plaq += 1.0-temp.trace().re/N;
    }
  }
  return Plaq;
}

/// Calculate the action
double gauge_action(field<SUN> *gauge, double beta){
  return beta*plaquette(gauge);
}





// Update the momentum with the gauge field
void gauge_force(field<SUN> *gauge, field<NMAT> *momentum, double eps){
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


// Update the gauge field with momentum
void gauge_step(field<SUN> *gauge, field<NMAT> *momentum, double eps){
  foralldir(dir){
    onsites(ALL){
      element<SUN> momexp = eps*momentum[dir][X];
      momexp.exp();
      gauge[dir][X] = momexp*gauge[dir][X];
    }
  }
}



void leapfrog_step(field<SUN> *gauge, field<NMAT> *momentum, double beta, double eps){
  gauge_step(gauge, momentum, 0.5*eps);
  gauge_force(gauge, momentum, beta*eps/N);
  gauge_step(gauge, momentum, 0.5*eps);
}





void update_hmc(field<SUN> *gauge, double beta, int steps){
  field<VEC> vector;
  field<NMAT> momentum[NDIM];

  gaussian_momentum(momentum);
  double S_gauge = gauge_action(gauge, beta);
  double S_mom = momentum_action(momentum);
  double start_action = S_gauge + S_mom;
  output0 << "Action " << S_gauge << " " << S_mom << " "
          << start_action << "\n";

  for(int step=0; step < steps; step++){

    leapfrog_step(gauge, momentum, beta, 0.1/steps);

    S_gauge = gauge_action(gauge, beta);
    S_mom = momentum_action(momentum);
  }
    output0 << "Action " << S_gauge << " " << S_mom << " "
            << S_gauge + S_mom  << " " << start_action-S_gauge - S_mom
            << "\n";
}








int main(int argc, char **argv){

  input parameters = input();
  parameters.import("parameters");
  double beta = parameters.get("beta");
  int seed = parameters.get("seed");
	double hmc_steps = parameters.get("hmc_steps");

  lattice->setup( nd[0], nd[1], nd[2], nd[3], argc, argv );
  seed_random(seed);


  field<SUN> gauge[NDIM];

  // Starting with a unit configuration
  foralldir(dir){
    onsites(ALL){
      gauge[dir][X].random();
    }
  }

  double eps = 0.001;
  double s1 = gauge_action(gauge, beta);
  output0 << s1 << "\n";

  SUN g = gauge[0].get_value_at(0);
  g.c[0][1].re += eps;
  g.c[1][0].re -= eps;
  gauge[0].set_value_at(g, 0);

  double s2 = gauge_action(gauge, beta);
  output0 << s2 << " " << (s2-s1)/(eps*beta) << "\n";
  g.c[0][1].re -= eps;
  g.c[1][0].re += eps;
  gauge[0].set_value_at(g, 0);

  field<NMAT> momentum[NDIM];
  foralldir(dir){
    onsites(ALL){
      momentum[dir][X] = 0;
    }
  }

  gauge_force(gauge, momentum, 1);
  NMAT f = momentum[0].get_value_at(0);
  output0 << f.c[0][1].re << " " << f.c[1][0].re << "\n";
  output0 << f.c[0][1].re + (s2-s1)/(eps*beta) << "\n";




  for(int step = 0; step < 20; step ++){
    update_hmc(gauge, beta, hmc_steps);
  }



  finishrun();

  return 0;
}
