#ifndef SMEARING_H
#define SMEARING_H


#include "datatypes/sun.h"
#include "hmc/gauge_field.h"


template<typename sun>
struct stout_smeared_field {
  using gauge_type = sun;
  using fund_type = sun;
  using basetype = typename sun::base_type;
  static constexpr int Nf = sun::size;
  static constexpr int N = sun::size;

  double c;
  int exp_steps = 5;

  gauge_field<sun> &fundamental;
  field<sun> gauge[NDIM];

  stout_smeared_field(gauge_field<fund_type>  &f, double _c) 
    : fundamental(f), c(_c){}
  stout_smeared_field(gauge_field<fund_type>  &f, double _c, int nstep) 
    : fundamental(f), c(_c), exp_steps(nstep) {}
  stout_smeared_field(stout_smeared_field &r)
    : fundamental(r.fundamental), c(r.c), exp_steps(r.exp_steps){}


  // Represent the fields
  void refresh(){
    foralldir(dir){
      gauge[dir].check_alloc();
    }
    foralldir(dir){
      field<sun> staple[NDIM];
      staple[dir] = calc_staples(fundamental.gauge, dir);
      onsites(ALL){
        staple[dir][X] = -c*fundamental.gauge[dir][X]*staple[dir][X];
        project_antihermitean(staple[dir][X]);
        staple[dir][X].exp(exp_steps);
        gauge[dir][X] = fundamental.gauge[dir][X]*staple[dir][X];
      }
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


  void add_momentum(field<squarematrix<N,cmplx<basetype>>> (&force)[NDIM]){
    field<squarematrix<N,cmplx<basetype>>> result[NDIM];
    field<sun> staple[NDIM];
    field<sun> Qfield[NDIM];
    field<sun> Lambda[NDIM];
    foralldir(dir){
      result[dir][ALL] = 0;
      staple[dir] = calc_staples(fundamental.gauge, dir);
      Qfield[dir] = calc_staples(fundamental.gauge, dir);
      onsites(ALL){
        Qfield[dir][X] = -c*fundamental.gauge[dir][X]*staple[dir][X];
        project_antihermitean(Qfield[dir][X]);
      }
    }
    foralldir(dir){
      onsites(ALL){
        element<sun> m0, m1, qn, iq, ex1, Q, R, L;
        Q = Qfield[dir][X];
        R = fundamental.gauge[dir][X];
        L = force[dir][X];
        Lambda[dir][X] = 0;
        m1=0; ex1=0; qn=1;
        iq = cmplx(0,1)*Q;
        double n=1.0;
        m0 = cmplx(0,1)*R*L;
        // gamma matrix (morningstar paper, eq 74)
        for(int k=0; k<exp_steps; k++){
          m1 = m0*qn+iq*m1;
          ex1 = ex1+n*qn; //ex1 = e^(iQ) 
          qn = qn*iq;
          // lambda = 1/(n+1)!*m1
          n=n*1.0/((double)(k+1));
          //output0 << " * " << 1.0/n << "\n":
          Lambda[dir][X] = c*Lambda[dir][X] + m1*n;
          project_antihermitean(Lambda[dir][X]);
        }
        result[dir][X] = R*L*ex1;
      }

      onsites(ALL){
        element<sun> lambda, stp, R, tmp;
        stp = staple[dir][X];
        R = fundamental.gauge[dir][X];
        lambda = Lambda[dir][X];

        result[dir][X] -= R*stp*lambda;
        Lambda[dir][X] = lambda*R;
      }
    }
    foralldir(dir1) foralldir(dir2){
      field<sun> stapleder2, stapleder3;
      onsites(ALL){
        element<sun> U1, U2, U3, U4, L, L2, tmp1, tmp2, tmp3, tmp4;
        U1 = fundamental.gauge[dir1][X];
        U2 = fundamental.gauge[dir2][X+dir1];
        U3 = fundamental.gauge[dir1][X+dir2];
        U4 = fundamental.gauge[dir2][X];
        L = Lambda[dir1][X];
        L2 = Lambda[dir1][X+dir2];

        stapleder2[X] = U2*(L*U4*U3).conjugate();
        stapleder3[X] = U3*U2.conjugate()*L*U4;
        result[dir2][X] -= U4*U3*U2.conjugate()*L;
        stapleder2[X] = U2*L*U4.conjugate()*U1;
        result[dir1][X] -= U1*U2*L.conjugate()*U4.conjugate();
        result[dir2][X] -= U4*L.conjugate()*U2.conjugate()*U1.conjugate();
      }
    }
    fundamental.add_momentum(result);
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