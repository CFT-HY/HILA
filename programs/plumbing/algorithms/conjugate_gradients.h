#ifndef CG_ALG
#define CG_ALG

///////////////////////////////////////////////////////
/// Generalized Conjugate gradient algorithm for fields
///
/// Idea is to solve for field2 in the equation 
/// field1 = operator * field2 where operator can be a 
/// staggered dirac operator or something else. 
///////////////////////////////////////////////////////

#include<iostream>

constexpr int CG_DEFAULT_MAXITERS = 10000;
constexpr double CG_DEFAULT_ACCURACY = 1e-8;


// Implement the apply function of the conjugate gradient operator
template<typename vector, typename op>
void GG_invert(vector &in, vector &out, op &M,
  double accuracy, int max_iters)
{
  vector r, p, Dp, DDp;
  double pDp = 0, rr = 0, rrnew = 0, rr_start=0;
  double alpha, beta;
  double target_rr, source_norm=0;
  
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

  for (int i = 0; i < max_iters; i++){
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
template<typename Op>
class CG{
  private:
    // The operator to invert
    Op & M; 
    // desired relative accuracy
    double accuracy = CG_DEFAULT_ACCURACY; 
    // maximum number of iterations
    double maxiters = CG_DEFAULT_MAXITERS;
  public:

    using vector_type = typename Op::vector_type;

    // Constructor: iniialize the operator
    CG(Op & op) : M(op) {};
    CG(Op & op, double _accuracy) : M(op) {accuracy = _accuracy;};
    CG(Op & op, double _accuracy, int _maxiters) : M(op) 
    {accuracy = _accuracy; maxiters=_maxiters;};

    void apply(field<vector_type> &in, field<vector_type> &out)
    {
      GG_invert(in, out, M, accuracy, maxiters);
    }
};

#endif
