#ifndef CG_ALG
#define CG_ALG

///////////////////////////////////////////////////////
/// Generalized Conjugate gradient algorithm for fields
///
/// Idea is to solve for field2 in the equation 
/// field1 = operator * field2 where operator can be a 
/// staggered dirac operator or something else. 
///////////////////////////////////////////////////////


#include <sstream>
#include<iostream>

constexpr int CG_DEFAULT_MAXITERS = 10000;
constexpr double CG_DEFAULT_ACCURACY = 1e-10;


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
      field<vector_type> r, p, Dp, DDp;
      r.copy_boundary_condition(in);
      p.copy_boundary_condition(in);
      Dp.copy_boundary_condition(in);
      DDp.copy_boundary_condition(in);
      double pDp = 0, rr = 0, rrnew = 0, rr_start=0;
      double alpha, beta;
      double target_rr, source_norm=0;

      onsites(M.par){
        source_norm += norm_squared(in[X]);
      }

      target_rr = accuracy*accuracy * source_norm;

      M.apply(out, Dp);
      M.dagger(Dp, DDp);
      onsites(M.par){
        r[X] = in[X] - DDp[X];
        p[X] = r[X];
      }

      onsites(M.par){
        rr += norm_squared(r[X]);
      }
      rr_start = rr;

      for (int i = 0; i < maxiters; i++){
        pDp=rrnew=0;
        M.apply(p, Dp);
        M.dagger(Dp, DDp);
        onsites(M.par){
          pDp += norm_squared(Dp[X]);
        }

        alpha=rr/pDp;

        onsites(M.par){
          out[X] = out[X] + alpha*p[X];
          r[X] = r[X] - alpha*DDp[X];
        }
        onsites(M.par){ 
          rrnew += norm_squared(r[X]);
        }
        #ifdef DEBUG_CG
        printf("CG iter %d, node %d, %g %g %g %g\n", i, mynode(), rrnew, rr_start, target_rr, pDp);
        #endif
        if( rrnew < target_rr )
            return;
        beta = rrnew/rr;
        p[M.par] = beta*p[X] + r[X];
        rr = rrnew;
      }

    }
};

#endif
