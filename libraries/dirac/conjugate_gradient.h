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
#include <iostream>

constexpr int CG_DEFAULT_MAXITERS = 10000;
constexpr double CG_DEFAULT_ACCURACY = 1e-12;

/// The conjugate gradient operator. Applies the inverse square of an operator on a vector
template <typename Op> class CG {
  private:
    // The operator to invert
    Op &M;
    // desired relative accuracy
    double accuracy = CG_DEFAULT_ACCURACY;
    // maximum number of iterations
    double maxiters = CG_DEFAULT_MAXITERS;

  public:
    /// Get the type the operator applies to
    using vector_type = typename Op::vector_type;

    /// Constructor: initialize the operator
    CG(Op &op) : M(op){};
    /// Constructor: operator and accuracy
    CG(Op &op, double _accuracy) : M(op) { accuracy = _accuracy; };
    /// Constructor: operator, accuracy and maximum number of iterations
    CG(Op &op, double _accuracy, int _maxiters) : M(op) {
        accuracy = _accuracy;
        maxiters = _maxiters;
    };

    /// The apply() -member runs the full conjugate gradient
    /// The operators themselves have the same structure.
    /// The conjugate gradient operator is Hermitean, so there is
    /// no dagger().
    void apply(Field<vector_type> &in, Field<vector_type> &out) {
        int i;
        struct timeval start, end;
        Field<vector_type> r, p, Dp, DDp;
        r.copy_boundary_condition(in);
        p.copy_boundary_condition(in);
        Dp.copy_boundary_condition(in);
        DDp.copy_boundary_condition(in);
        out.copy_boundary_condition(in);
        double pDp = 0, rr = 0, rrnew = 0, rr_start = 0;
        double alpha, beta;
        double target_rr, source_norm = 0;

        gettimeofday(&start, NULL);

        onsites(M.par) { source_norm += norm_squared(in[X]); }

        target_rr = accuracy * accuracy * source_norm;

        M.apply(out, Dp);
        M.dagger(Dp, DDp);
        onsites(M.par) {
            r[X] = in[X] - DDp[X];
            p[X] = r[X];
        }

        onsites(M.par) { rr += norm_squared(r[X]); }
        rr_start = rr;

        for (i = 0; i < maxiters; i++) {
            pDp = rrnew = 0;
            M.apply(p, Dp);
            M.dagger(Dp, DDp);
            onsites(M.par) { pDp += norm_squared(Dp[X]); }

            alpha = rr / pDp;

            onsites(M.par) {
                out[X] = out[X] + alpha * p[X];
                r[X] = r[X] - alpha * DDp[X];
            }
            onsites(M.par) { rrnew += norm_squared(r[X]); }
#ifdef DEBUG_CG
            output0 << "CG step " i << ", residue " << sqrt(rrnew / target_rr) << "\n";
#endif
            if (rrnew < target_rr)
                break;
            beta = rrnew / rr;
            p[M.par] = beta * p[X] + r[X];
            rr = rrnew;
        }

        gettimeofday(&end, NULL);
        double timing =
            1e-3 * (end.tv_usec - start.tv_usec) + 1e3 * (end.tv_sec - start.tv_sec);

        output0 << "Conjugate Gradient: " << i << " steps in " << timing << "ms, ";
        output0 << "relative residue:" << rrnew / source_norm << "\n";
    }
};

#endif
