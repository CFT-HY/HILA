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
  int smear_steps = 1;
  int exp_steps = 10;

  gauge_field<sun> &fundamental;
  field<sun> gauge[NDIM];

  stout_smeared_field(gauge_field<fund_type>  &f, double coeff) 
    : fundamental(f), c(coeff){}
  stout_smeared_field(gauge_field<fund_type>  &f, double coeff, int nsteps)
    : fundamental(f), c(coeff), smear_steps(nsteps) {}
  stout_smeared_field(gauge_field<fund_type>  &f, double coeff, int nsteps, int expsteps)
    : fundamental(f), c(coeff), exp_steps(expsteps) {}
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
      onsites(ALL){
        Qfield[dir][X] = -c*fundamental.gauge[dir][X]*staple[dir][X];
        project_antihermitean(Qfield[dir][X]);
      }
    }
    foralldir(dir){
      onsites(ALL){
        element<sun> m0, m1, qn, ex1, Q, R, L;
        Q = Qfield[dir][X];
        R = fundamental.gauge[dir][X];
        L = force[dir][X];
        //L = gauge[dir][X].conjugate()*force[dir][X];
        ex1=1; qn=Q;
        double n=1.0;
        m0 = L*R;
        m1 = m0;
        if(exp_steps==0){
          Lambda[dir][X] = 0;
        } else {
          Lambda[dir][X] = m1;
          ex1=ex1+Q;
        }
        // gamma matrix (morningstar paper, eq 74)
        for(int k=2; k<=exp_steps; k++){
          n=n*1.0/((double)k);
          m1 = m0*qn+Q*m1;
          qn = qn*Q;
          ex1 = ex1+n*qn; //ex1 = e^(iQ)
          // lambda = 1/(n+1)!*m1
          //output0 << " * " << 1.0/n << "\n";
          
          Lambda[dir][X] = Lambda[dir][X] + m1*n;
        }
        project_antihermitean(Lambda[dir][X]);
        result[dir][X] = ex1*L;
      }

      onsites(ALL){
        element<sun> lambda, stp, R, tmp;
        stp = staple[dir][X];
        R = fundamental.gauge[dir][X];
        lambda = Lambda[dir][X];

        result[dir][X] -= c*stp*lambda;
        Lambda[dir][X] = -c*lambda*R;
      }
    }
    foralldir(dir1) foralldir(dir2) if(dir1!=dir2){
      field<sun> stapleder2, stapleder3;
      onsites(ALL){
        element<sun> U1, U2, U3, U4, L, L2;
        U1 = fundamental.gauge[dir1][X];
        U2 = fundamental.gauge[dir2][X+dir1];
        U3 = fundamental.gauge[dir1][X+dir2];
        U4 = fundamental.gauge[dir2][X];
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

      onsites(ALL){
        result[dir2][X] += stapleder2[X-dir1];
        result[dir1][X] += stapleder3[X-dir2];
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