#ifndef CG_ALG
#define CG_ALG

///////////////////////////////////////////////////////
/// Generalized Conjugate gradient algorithm for fields
///
/// Idea is to solve for field2 in the equation 
/// field1 = operator * field2 where operator can be a 
/// staggered dirac operator or something else. 
///////////////////////////////////////////////////////

#include "../dirac.h"
#include<iostream>

#define MAXITERS 10000


template<typename vector, typename op>
void GG_invert(vector &in, vector &out, op &M){
  vector r, p, Dp, DDp;
  double pDp = 0, rr = 0, rrnew = 0, rr_start=0;
  double alpha, beta;
  double target_rr, source_norm=0;
  double accuracy = 1e-8;
  
  onsites(ALL){
    source_norm += norm_squared(in[X]);
  }
        
  target_rr = accuracy*accuracy * source_norm;

  M.apply(out, Dp);
  M.dagger(Dp, DDp);
  onsites(ALL){
      r[X] = in[X] - DDp[X];
      p[X] = r[X];
  }

  onsites(ALL){
      rr += (r[X]*r[X]).re;
  }
  rr_start = rr;

  for (int i = 0; i < MAXITERS; i++){
    pDp=rrnew=0;
    M.apply(p, Dp);
    M.dagger(Dp, DDp);
    onsites(ALL){
        pDp += (p[X]*DDp[X]).re;
    }

    alpha=rr/pDp;

    onsites(ALL){
        out[X] = out[X] + alpha*p[X];
        r[X] = r[X] - alpha*DDp[X];
    }
    onsites(ALL){ 
        rrnew += (r[X]*r[X]).re;
    }
    #ifdef DEBUG
    printf("CG iter %d, node %d, %g %g %g %g\n", i, mynode(), rrnew, rr_start, target_rr, pDp);
    #endif
    if( rrnew < target_rr )
        return;
    beta = rrnew/rr;
    p[ALL] = beta*p[X] + r[X];
    rr = rrnew;
  }
}


template<typename vector, typename Op>
class CG{
  private:
    Op & M; // The operator to invert
  public:

    // Constructor: initialize mass, gauge and eta
    CG(Op & op) : M(op) {};

    void apply(vector &in, vector &out)
    {
      GG_invert(in, out, M);
    }
};

#endif
