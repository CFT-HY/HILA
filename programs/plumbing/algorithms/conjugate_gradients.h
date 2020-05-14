#ifndef CG_ALG
#define CG_ALG

///////////////////////////////////////////////////////
/// Generalized Conjugate gradient algorithm for fields
///
/// Idea is to solve for field2 in the equation 
/// field1 = operator * field2 where operator can be a 
/// staggered dirac operator or something else. 
///////////////////////////////////////////////////////

#include "../fermion/staggered.h"
#include<iostream>

#define MAXITERS 10000


// Implement the apply function of the conjugate gradient operator
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
      rr += r[X].rdot(r[X]);
  }
  rr_start = rr;

  for (int i = 0; i < MAXITERS; i++){
    pDp=rrnew=0;
    M.apply(p, Dp);
    M.dagger(Dp, DDp);
    onsites(ALL){
        pDp += p[X].rdot(DDp[X]);
    }

    alpha=rr/pDp;

    onsites(ALL){
        out[X] = out[X] + alpha*p[X];
        r[X] = r[X] - alpha*DDp[X];
    }
    onsites(ALL){ 
        rrnew += r[X].rdot(r[X]);
    }
    #ifdef DEBUG_CG
    printf("CG iter %d, node %d, %g %g %g %g\n", i, mynode(), rrnew, rr_start, target_rr, pDp);
    #endif
    if( rrnew < target_rr )
        return;
    beta = rrnew/rr;
    p[ALL] = beta*p[X] + r[X];
    rr = rrnew;
  }
}


/// The conjugate gradient operator. Applies square the inverse of an operator on a vector
template<typename vector, typename Op>
class CG{
  private:
    Op & M; // The operator to invert
  public:

    // Constructor: iniialize the operator
    CG(Op & op) : M(op) {};

    void apply(vector &in, vector &out)
    {
      GG_invert(in, out, M);
    }
};

#endif
