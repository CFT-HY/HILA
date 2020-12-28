#ifndef SMEARING_H
#define SMEARING_H


#include "datatypes/sun.h"
#include "hmc/gauge_field.h"



/// Calculate the exponential of Q and the matrix lambda=d/dQ (e^Q m0) 
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


/// Calculate the derivative of with respect to the links a positive and negative staple and add to result
template<typename matrix, typename forcetype>
void staple_dir_derivative(field<matrix> &basegauge1, field<matrix> &basegauge2, field<matrix> &Lambda, field<forcetype> &result1, field<forcetype> &result2, direction dir1, direction dir2){
  field<matrix> stapleder2, stapleder3; // Two derivatives that need to be communicated
  
  onsites(ALL){
    element<matrix> U1, U2, U3, U4, L, L2;
    U1 = basegauge1[X];
    U2 = basegauge2[X+dir1];
    U3 = basegauge1[X+dir2];
    U4 = basegauge2[X];
    L = Lambda[X];
    L2 = Lambda[X+dir2];

    // Up staple
    result2[X] += (L*U2*U3.conjugate()).conjugate();
    stapleder2[X] = U3.conjugate()*U4.conjugate()*L;
    stapleder3[X] = (U4.conjugate()*L*U2).conjugate();

    // Down staple
    stapleder2[X] = stapleder2[X] + L2.conjugate()*U4.conjugate()*U1;
    result1[X] += U2*L2.conjugate()*U4.conjugate();
    result2[X] += L2*U2.conjugate()*U1.conjugate();
  }

  // Move derivatives up where necessary
  onsites(ALL){
    result2[X] += stapleder2[X-dir1];
    result1[X] += stapleder3[X-dir2];
  }
}



/// A stout smeared gauge field built of a standard gauge field.
/// The structure is the same as for the standard and represented
/// gauge fields.
/// refresh(): recalculates the smeared field from the underlying
///    gauge field.
/// add_momentum(): transforms a derivative with respect to this
///    field to a derivative with respect to the underlying gauge
///    field and add to the momentum of the gauge field.
template<typename sun>
class stout_smeared_field : public gauge_field_base<sun>{
  public:
  using gauge_type = sun;
  using fund_type = sun;
  using basetype = typename sun::base_type;
  static constexpr int Nf = sun::size;
  static constexpr int N = sun::size;

  double c;
  int smear_steps = 1;
  int exp_steps = 10;

  gauge_field<sun> &base_field;
  field<sun> **staples;
  field<sun> **smeared_fields;

  stout_smeared_field(gauge_field<fund_type>  &f, double coeff) 
    : base_field(f), c(coeff){
      gauge_field_base<sun>();
      allocate();
    }
  stout_smeared_field(gauge_field<fund_type>  &f, double coeff, int nsteps)
    : base_field(f), c(coeff), smear_steps(nsteps) {
      gauge_field_base<sun>();
      allocate();
    }
  stout_smeared_field(gauge_field<fund_type>  &f, double coeff, int nsteps, int expsteps)
    : base_field(f), c(coeff), smear_steps(nsteps), exp_steps(expsteps) {
      gauge_field_base<sun>();
      allocate();
    }
  stout_smeared_field(stout_smeared_field &r)
    : base_field(r.base_field), c(r.c), smear_steps(r.smear_steps), exp_steps(r.exp_steps){
      gauge_field_base<sun>();
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
    smeared_fields[smear_steps-1]  = &(this->gauge[0]);
  }

  ~stout_smeared_field() {
    for(int step=0; step<smear_steps-1; step++){
      delete staples[step];
      delete smeared_fields[step];
    }
    free(staples); free(smeared_fields);
  }


  // Represent the fields
  void refresh(){
    field<sun> *previous;
    previous = &base_field.gauge[0];

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
    base_field.set_unity();
    refresh();
  }

  void random(){
    base_field.random();
    refresh();
  }


  void add_momentum(field<SquareMatrix<N,cmplx<basetype>>> *force){
    // Two storage fields for the current and previous levels of the force
    field<SquareMatrix<N,cmplx<basetype>>> storage1[NDIM];
    field<SquareMatrix<N,cmplx<basetype>>> storage2[NDIM];
    foralldir(dir){
      storage1[dir] = force[dir];
    }
    field<SquareMatrix<N,cmplx<basetype>>> *previous= &storage1[0];
    field<SquareMatrix<N,cmplx<basetype>>> *result= &storage2[0];

    // Another storage field, for the derivative of the exponential
    field<sun> Lambda[NDIM];

    for(int step=smear_steps-1; step>=0; step--){
      // Find the gauge field the current level is calculated from
      field<sun> *basegauge;
      if(step==0){
        basegauge = &base_field.gauge[0];
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
        staple_dir_derivative(basegauge[dir1], basegauge[dir2], Lambda[dir1], result[dir1], result[dir2], dir1, dir2);
      }

      // Swap previous and result for the next iteration
      field<SquareMatrix<N,cmplx<basetype>>> *tmp = previous;
      previous = result;
      result = tmp;
      basegauge = &smeared_fields[step][0];
    }

    // Since we swap at the end, the force is now in "previous"
    base_field.add_momentum(previous);
  }

  void draw_momentum(){
    base_field.draw_momentum();
  }
  void zero_momentum(){
    base_field.zero_momentum();
  }

  void backup(){
    base_field.backup();
  }

  // Restore the previous backup
  void restore_backup(){
    base_field.restore_backup();
  }


  field<sun> & get_momentum(int dir){
    return base_field.get_momentum(dir);
  }
  field<sun> & get_gauge(int dir){
    return base_field.get_gauge(dir);
  }
};






#if NDIM==4




template<typename sun>
struct HEX_smeared_field : public gauge_field_base<sun> {
  using gauge_type = sun;
  using fund_type = sun;
  using basetype = typename sun::base_type;
  static constexpr int Nf = sun::size;
  static constexpr int N = sun::size;

  // SU2 default parameters
  double c1=0.13;
  double c2=0.1525;
  double c3=0.175;
  int exp_steps = 10;

  gauge_field<sun> &base_field;
  field<sun> staples3[NDIM][NDIM];
  field<sun> level3[NDIM][NDIM];
  field<sun> staples2[NDIM][NDIM];
  field<sun> level2[NDIM][NDIM];
  field<sun> staples1[NDIM];

  HEX_smeared_field(gauge_field<fund_type>  &f)
    : base_field(f) {
      gauge_field_base<sun>();
    }
  HEX_smeared_field(gauge_field<fund_type>  &f, int nsteps)
    : base_field(f), exp_steps(nsteps) {
      gauge_field_base<sun>();
    }
  HEX_smeared_field(gauge_field<fund_type>  &f, double _c1, double _c2, double _c3)
    : base_field(f), c1(_c1), c2(_c2), c3(_c3) {
      gauge_field_base<sun>();
    }
  HEX_smeared_field(gauge_field<fund_type>  &f, double _c1, double _c2, double _c3, int nsteps)
    : base_field(f), c1(_c1), c2(_c2), c3(_c3), exp_steps(nsteps) {
      gauge_field_base<sun>();
    }


  // Represent the fields
  void refresh(){
    field<sun> *previous;
    previous = &base_field.gauge[0];

    foralldir(mu){
      base_field.gauge[mu].check_alloc();
    }

    // Level 3, link to direction mu, staples only summed to
    // to direction nu
    foralldir(mu) foralldir(nu)  if(mu!=nu){
      staples3[nu][mu] = calc_staples(base_field.gauge, base_field.gauge, mu, nu);
      onsites(ALL){
        element<sun> Q;
        Q = -c3*base_field.gauge[mu][X]*staples3[nu][mu][X];
        project_antihermitean(Q);
        Q.exp(exp_steps);
        level3[nu][mu][X] = base_field.gauge[mu][X]*Q;
      }
    }

    // Level 2, link to direction mu, staples summed to direction
    // rho != mu, nu. label directions nu and rho by eta != mu, nu, rho
    foralldir(mu) foralldir(nu) if(mu!=nu) {
      staples2[nu][mu][ALL] = 0;
      foralldir(rho) if(rho!=mu) if(rho!=nu){
        direction eta;
        foralldir(e) if(e!=mu && e!=nu && e!=rho) eta = e;
        field<sun> stp = calc_staples(level3[eta], level3[eta], mu, rho);
        staples2[nu][mu][ALL] = staples2[nu][mu][X] + stp[X];
      }
      onsites(ALL){
        element<sun> Q;
        Q = -c2*base_field.gauge[mu][X]*staples2[nu][mu][X];
        project_antihermitean(Q);
        Q.exp(exp_steps);
        level2[nu][mu][X] = base_field.gauge[mu][X]*Q;
      }
    }

    // Level 1, link to direction mu, staples summed to directions nu
    // with direction nu excluded from lower levels
    foralldir(mu) {
      staples1[mu][ALL] = 0;
      foralldir(nu) if(mu!=nu) {
        field<sun> stp = calc_staples(level2[nu], level2[mu], mu, nu);
        staples1[mu][ALL] = staples1[mu][X] + stp[X];
      }
      onsites(ALL){
        element<sun> Q;
        Q = -c1*base_field.gauge[mu][X]*staples1[mu][X];
        project_antihermitean(Q);
        Q.exp(exp_steps);
        this->gauge[mu][X] = base_field.gauge[mu][X]*Q;
      }
    }

  }

  void set_unity(){
    base_field.set_unity();
    refresh();
  }

  void random(){
    base_field.random();
    refresh();
  }


  void add_momentum(field<SquareMatrix<N,cmplx<basetype>>> *force){
    field<sun> lambda1[NDIM];
    field<SquareMatrix<N,cmplx<basetype>>> result1[NDIM][NDIM];
    field<sun> lambda2[NDIM][NDIM];
    field<SquareMatrix<N,cmplx<basetype>>> result2[NDIM][NDIM];
    field<sun> lambda3[NDIM][NDIM];
    field<SquareMatrix<N,cmplx<basetype>>> result[NDIM];

    foralldir(mu) {
      result[mu] = 0;
      foralldir(nu) {
        result1[mu][nu] = 0;
        result2[mu][nu] = 0;
      }
    }
    
    // Level1 exponential
    foralldir(mu){
      onsites(ALL){
        element<sun> m0, m1, qn, eQ, Q;
        Q = -c1*base_field.gauge[mu][X]*staples1[mu][X];
        project_antihermitean(Q);

        m0 = force[mu][X]*base_field.gauge[mu][X];
        exp_and_derivative(Q, m0, lambda1[mu][X], eQ, exp_steps);
        project_antihermitean(lambda1[mu][X]);

        // First derivative term, R in R*exp(Q)*L
        result[mu][X] = eQ*force[mu][X];

        // second derivative term, the first link in the plaquette
        result[mu][X] -= c1*staples1[mu][X]*lambda1[mu][X];
        
        // Now update Lambda to the derivative with respect to the staple
        lambda1[mu][X] = -c1*lambda1[mu][X]*base_field.gauge[mu][X];
      }
    }

    // level1 staple
    foralldir(mu) foralldir(nu) if(mu!=nu){
      staple_dir_derivative(level2[nu][mu], level2[mu][nu], lambda1[mu], result1[nu][mu], result1[mu][nu], mu, nu);
    }



    // level2 exponential
    foralldir(mu) foralldir(nu) if(mu!=nu) {
      onsites(ALL){
        element<sun> m0, m1, qn, eQ, Q;
        Q = -c2*base_field.gauge[mu][X]*staples2[nu][mu][X];
        project_antihermitean(Q);

        m0 = result1[nu][mu][X]*base_field.gauge[mu][X];
        exp_and_derivative(Q, m0, lambda2[nu][mu][X], eQ, exp_steps);
        project_antihermitean(lambda2[nu][mu][X]);

        // First derivative term, R in R*exp(Q)*L
        result[mu][X] += eQ*result1[nu][mu][X];

        // second derivative term, the first link in the plaquette
        result[mu][X] -= c2*staples2[nu][mu][X]*lambda2[nu][mu][X];
        
        // Now update Lambda to the derivative with respect to the staple
        lambda2[nu][mu][X] = -c2*lambda2[nu][mu][X]*base_field.gauge[mu][X];
      }
    }

    // level2 staple
    foralldir(mu) foralldir(nu) if(mu!=nu){
      foralldir(rho) if(rho!=mu) if(rho!=nu){
        direction eta;
        foralldir(e) if(e!=mu && e!=nu && e!=rho) eta = e;
        staple_dir_derivative(level3[eta][mu], level3[eta][rho], lambda2[nu][mu], result2[eta][mu], result2[eta][rho], mu, rho);
      }
    }



    // level3 exponential
    foralldir(mu) foralldir(nu) if(mu!=nu) {
      onsites(ALL){
        element<sun> m0, m1, qn, eQ, Q;
        Q = -c3*base_field.gauge[mu][X]*staples3[nu][mu][X];
        project_antihermitean(Q);

        m0 = result2[nu][mu][X]*base_field.gauge[mu][X];
        exp_and_derivative(Q, m0, lambda3[nu][mu][X], eQ, exp_steps);
        project_antihermitean(lambda3[nu][mu][X]);

        // First derivative term, R in R*exp(Q)*L
        result[mu][X] += eQ*result2[nu][mu][X];

        // second derivative term, the first link in the plaquette
        result[mu][X] -= c3*staples3[nu][mu][X]*lambda3[nu][mu][X];
        
        // Now update Lambda to the derivative with respect to the staple
        lambda3[nu][mu][X] = -c3*lambda3[nu][mu][X]*base_field.gauge[mu][X];
      }
    }

    // level3 staple
    foralldir(mu) foralldir(nu) if(mu!=nu){
      staple_dir_derivative(base_field.gauge[mu], base_field.gauge[nu], lambda3[nu][mu], result[mu], result[nu], mu, nu);
    }


    // Add to the base gauge momentum
    base_field.add_momentum(result);
  }

  void draw_momentum(){
    base_field.draw_momentum();
  }
  void zero_momentum(){
    base_field.zero_momentum();
  }

  void backup(){
    base_field.backup();
  }

  // Restore the previous backup
  void restore_backup(){
    base_field.restore_backup();
  }


  field<sun> & get_momentum(int dir){
    return base_field.get_momentum(dir);
  }
  field<sun> & get_gauge(int dir){
    return base_field.get_gauge(dir);
  }
};


#endif





#endif