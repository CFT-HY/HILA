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

#define MAXITERS 1
//valid operators for CG engine 
struct staggered_dirac;

template<typename T>
class CG_engine{
    public:
    template<typename mtype, typename vtype> 
    void solve(
        const mtype gauge[NDIM],
        const double mass,
        const vtype &v_in,
        vtype &v_out)
    {
        std::cout << "Incorrect args in call to CG_engine.solve(..args..) for the selected operator";
    }
    template<typename mtype, typename vtype>
    void act(
        const mtype gauge[NDIM],
        const double mass,
        const vtype &v_in,
        vtype &v_out)
    {}
};

template<>
class CG_engine<staggered_dirac>{
    public:
    template<typename mtype, typename vtype> 
    void solve(
        const mtype gauge[NDIM],
        const double mass,
        const vtype &b,
        vtype &sol)
    {
        vtype r, p, Dp;
        double pDDp = 0, rr = 0, rrnew = 0;
        double alpha, beta;

        onsites(ALL){
            r[X] = b[X];
            p[X] = b[X];
        }

        for (int i = 0; i < MAXITERS; i++){
            dirac_stagggered(gauge, mass, p, Dp); //give current p, get transformed version in Dp
            rr=pDDp=0;

            onsites(ALL){
                rr += norm_squared(r[X]);
                pDDp += norm_squared(Dp[X]);
            }

            alpha = rr / pDDp;
            rrnew = 0;

            onsites(ALL){
                sol[X] = sol[X] + alpha*p[X]; //update solution vector
                r[X] = r[X] - alpha*Dp[X]; //update residuals
                rrnew += norm_squared(r[X]);
            }

            if (rrnew < 10) return;

            beta = rrnew/rr;
            p[ALL] = r[X] + beta*p[X];
        }
    }

    template<typename mtype, typename vtype>
    void act(
        const mtype gauge[NDIM],
        const double mass,
        const vtype &v_in,
        vtype &v_out)
    {
        dirac_stagggered_alldim(gauge,mass,v_in,v_out);
    }
};

#endif