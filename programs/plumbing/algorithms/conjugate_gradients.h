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
        vtype &b,
        vtype &x_0)
    {
        vtype r, p, Dp;
        double pDDp = 0, rr = 0, rrnew = 0;
        double alpha, beta;
        double target_rr, source_norm=0;
        double accuracy = 1e-8;
        
        onsites(ALL){
            source_norm += norm_squared(b[X]);
        }
        
        target_rr = accuracy*accuracy * source_norm;

        dirac_stagggered(gauge, mass, x_0, Dp);
        onsites(ALL){
            r[X] = b[X] - Dp[X];
            p[X] = r[X];
        }

        onsites(ALL){
            rr += norm_squared(r[X]);
        }

        for (int i = 0; i < MAXITERS; i++){
            pDDp=rrnew=0;
            dirac_stagggered(gauge, mass, p, Dp);
            //note: unsure about whether it should be pDDp or pDp  
            onsites(ALL){
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
            #ifdef DEBUG
            printf("CG iter %d, node %d, %g %g %g\n", i, mynode(), alpha, rrnew, target_rr);
            #endif
            if( rrnew < target_rr )
                return;
            beta = rrnew/rr;
            p[ALL] = beta*p[X] + r[X];
            rr = rrnew;
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
