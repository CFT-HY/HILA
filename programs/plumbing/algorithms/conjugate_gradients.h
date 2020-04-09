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

#define MAXITERS 10
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

template<> //not optimized yet
class CG_engine<staggered_dirac>{
    public:
    template<typename mtype, typename vtype> 
    void solve(
        const mtype gauge[NDIM],
        const double mass,
        const vtype &b,
        vtype &x_0)
    {
        vtype r, p, Dp;
        double pDDp = 0, rr = 0, rrnew = 0;
        double alpha, beta;
        dirac_stagggered(gauge, mass, x_0, Dp);
        onsites(ALL){
            r[X] = b[X] - Dp[X];
            p[X] = b[X];
        }

        for (int i = 0; i < MAXITERS; i++){
            dirac_stagggered(gauge, mass, p, Dp);
            rr=pDDp=rrnew=0;
            //note: unsure about whether it should be pDDp or pDp  
            onsites(ALL){
                rr += norm_squared(r[X]);
                pDDp += norm_squared(p[X].conjugate()*Dp[X]);
            }

            alpha=rr/pDDp;

            onsites(ALL){
                x_0[X] = x_0[X] + alpha*p[X];
                r[X] = r[X] - alpha*Dp[X]; 
            }
            onsites(ALL){ 
                rrnew += norm_squared(r[X]); 
            }
            beta = rrnew/rr;
            p[ALL] = beta*p[X] + r[X];
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
