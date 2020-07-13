#ifndef SMEARING_H
#define SMEARING_H


#include "datatypes/sun.h"
#include "hmc/gauge_field.h"



// Calculate the exponential of Q and the matrix lambda=d/dQ (e^Q m0) 
template<typename sun>
void exp_and_derivative(sun &Q, sun &m0, sun &lambda, sun &eQ, int exp_steps){
  sun m1, qn;
  eQ=1; qn=Q;
  double n=1.0;
  m1 = m0;
  if(exp_steps==0){
    lambda = 0;
  } else {
    lambda = m1;
    eQ=eQ+Q;
  }
  // gamma matrix (morningstar paper, eq 74)
  for(int k=2; k<=exp_steps; k++){
    n=n*1.0/((double)k);
    m1 = m0*qn+Q*m1;
    qn = qn*Q;
    eQ = eQ+n*qn;
    // lambda = 1/(n+1)!*m1
    //output0 << " * " << 1.0/n << "\n";
    lambda = lambda + m1*n;
  }
}


// Calculate the derivative of with respect to the links a positive and negative staple and add to result
template<typename matrix, typename forcetype>
void staple_dir_derivative(field<matrix> *basegauge, field<matrix> (&Lambda)[NDIM], field<forcetype> *result, direction dir1, direction dir2){
  field<matrix> stapleder2, stapleder3; // Two derivatives that need to be communicated
  
  onsites(ALL){
    element<matrix> U1, U2, U3, U4, L, L2;
    U1 = basegauge[dir1][X];
    U2 = basegauge[dir2][X+dir1];
    U3 = basegauge[dir1][X+dir2];
    U4 = basegauge[dir2][X];
    L = Lambda[dir1][X];
    L2 = Lambda[dir1][X+dir2];

    // Up staple
    result[dir2][X] += (L*U2*U3.conjugate()).conjugate();
    stapleder2[X] = U3.conjugate()*U4.conjugate()*L;
    stapleder3[X] = (U4.conjugate()*L*U2).conjugate();

    // Down staple
    stapleder2[X] = stapleder2[X] + L2.conjugate()*U4.conjugate()*U1;
    result[dir1][X] += U2*L2.conjugate()*U4.conjugate();
    result[dir2][X] += L2*U2.conjugate()*U1.conjugate();
  }

  // Move derivatives up where necessary
  onsites(ALL){
    result[dir2][X] += stapleder2[X-dir1];
    result[dir1][X] += stapleder3[X-dir2];
  }
}






template<typename sun>
struct stout_smeared_field {
  using gauge_type = sun;
  using fund_type = sun;
  using basetype = typename sun::base_type;
  static constexpr int Nf = sun::size;
  static constexpr int N = sun::size;

  double c;
  int smear_steps = 1;
  int exp_steps = 10;

  gauge_field<sun> &fundamental;
  field<sun> gauge[NDIM];
  field<sun> **staples;
  field<sun> **smeared_fields;

  stout_smeared_field(gauge_field<fund_type>  &f, double coeff) 
    : fundamental(f), c(coeff){
      allocate();
    }
  stout_smeared_field(gauge_field<fund_type>  &f, double coeff, int nsteps)
    : fundamental(f), c(coeff), smear_steps(nsteps) {
      allocate();
    }
  stout_smeared_field(gauge_field<fund_type>  &f, double coeff, int nsteps, int expsteps)
    : fundamental(f), c(coeff), smear_steps(nsteps), exp_steps(expsteps) {
      allocate();
    }
  stout_smeared_field(stout_smeared_field &r)
    : fundamental(r.fundamental), c(r.c), smear_steps(r.smear_steps), exp_steps(r.exp_steps){
      allocate();
    }

  void allocate(){
    staples = (field<sun>**)malloc(smear_steps*sizeof(field<sun>*));
    smeared_fields = (field<sun>**)malloc(smear_steps*sizeof(field<sun>*));
    for(int step=0; step<smear_steps-1; step++){
      staples[step] = new field<sun>[NDIM];
      smeared_fields[step] = new field<sun>[NDIM];
    }
    staples[smear_steps-1]  = new field<sun>[NDIM];
    smeared_fields[smear_steps-1]  = &(gauge[0]);
  }


  // Represent the fields
  void refresh(){
    field<sun> *previous;
    previous = &fundamental.gauge[0];

    for(int step=0; step<smear_steps; step++){
      foralldir(dir){
        previous[dir].check_alloc();
      }
      foralldir(dir){
        staples[step][dir] = calc_staples(previous, dir);
        onsites(ALL){
          element<sun> Q;
          Q = -c*previous[dir][X]*staples[step][dir][X];
          project_antihermitean(Q);
          Q.exp(exp_steps);
          smeared_fields[step][dir][X] = previous[dir][X]*Q;
        }
      }
      previous = smeared_fields[step];
    }
  }

  void set_unity(){
    fundamental.set_unity();
    refresh();
  }

  void random(){
    fundamental.random();
    refresh();
  }


  void add_momentum(field<squarematrix<N,cmplx<basetype>>> *force){
    // Two storage fields for the current and previous levels of the force
    field<squarematrix<N,cmplx<basetype>>> storage1[NDIM];
    field<squarematrix<N,cmplx<basetype>>> storage2[NDIM];
    foralldir(dir){
      storage1[dir] = force[dir];
    }
    field<squarematrix<N,cmplx<basetype>>> *previous= &storage1[0];
    field<squarematrix<N,cmplx<basetype>>> *result= &storage2[0];

    // Another storage field, for the derivative of the exponential
    field<sun> Lambda[NDIM];

    for(int step=smear_steps-1; step>=0; step--){
      // Find the gauge field the current level is calculated from
      field<sun> *basegauge;
      if(step==0){
        basegauge = &fundamental.gauge[0];
      } else {
        basegauge = smeared_fields[step-1];
      } 

      // Take the derivative of the exponential
      foralldir(dir){
        result[dir][ALL] = 0;
        onsites(ALL){
          element<sun> m0, m1, qn, eQ, Q;
          Q = -c*basegauge[dir][X]*staples[step][dir][X];
          project_antihermitean(Q);

          m0 = previous[dir][X]*basegauge[dir][X];
          exp_and_derivative(Q, m0, Lambda[dir][X], eQ, exp_steps);
          
          project_antihermitean(Lambda[dir][X]);

          // First derivative term, R in R*exp(Q)*L
          result[dir][X] = eQ*previous[dir][X]; 

          // second derivative term, the first link in the plaquette
          result[dir][X] -= c*staples[step][dir][X]*Lambda[dir][X];
          
          // Now update Lambda to the derivative of the staple
          Lambda[dir][X] = -c*Lambda[dir][X]*basegauge[dir][X];
        }
      }

      // Take the derivetive with respect to the links in the staples
      foralldir(dir1) foralldir(dir2) if(dir1!=dir2){
        staple_dir_derivative(basegauge, Lambda, result, dir1, dir2);
      }

      // Swap previous and result for the next iteration
      field<squarematrix<N,cmplx<basetype>>> *tmp = previous;
      previous = result;
      result = tmp;
      basegauge = &smeared_fields[step][0];
    }

    // Since we swap at the end, the force is now in "previous"
    fundamental.add_momentum(previous);
  }

  void draw_momentum(){
    fundamental.draw_momentum();
  }
  void zero_momentum(){
    fundamental.zero_momentum();
  }

  void backup(){
    fundamental.backup();
  }

  // Restore the previous backup
  void restore_backup(){
    fundamental.restore_backup();
  }


  field<sun> & get_momentum(int dir){
    return fundamental.get_momentum(dir);
  }
  field<sun> & get_gauge(int dir){
    return fundamental.get_gauge(dir);
  }
};




#endif