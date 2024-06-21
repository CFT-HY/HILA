/** @file gradient_flow.h */

#ifndef GRADIENT_FLOW_H_
#define GRADIENT_FLOW_H_

#include "hila.h"
#include "energy_and_topo_charge_clover.h"
#include "energy_and_topo_charge_log.h"
#include "wilson_plaquette_action.h"
#include "bulk_prevention_action.h"
#include "wilson_line_and_force.h"
#include "improved_action.h"
#include "clover_action.h"
#include "string_format.h"



template <typename group,typename atype=hila::arithmetic_type<group>>
void get_gf_force(const GaugeField<group>& U,VectorField<Algebra<group>>& E) {
    // wrapper for force computation routine to be used for gradient flow

    atype eps=2.0; // factor 2.0 to switch from unoriented to oriented plaquettes
                   // (factor is usually absorbed in \beta, but gradient flow force is
                   // computed from action term with \beta-factor stripped off) 


    //get_force_wplaq(U,E,eps); // force for Wilson plaquette action

    //get_force_bp(U,E,eps);    // force for bulk-preventing (BP) action (should only be used
                                // for gradient flow if BP is also used for HMC evolution)


    // improved action:
    // DBW2 action parameters:
    //atype c12=-1.4088;     // rectangle weight
    //atype c11=1.0-8.0*c12; // plaquette weight

    // Iwasaki action parameters:
    atype c12=-0.331;     // rectangle weight
    atype c11=1.0-8.0*c12; // plaquette weight

    // LW action parameters:
    //atype c12=-1.0/12.0;     // rectangle weight
    //atype c11=1.0-8.0*c12; // plaquette weight

    // Wilson plaquette action parameters:
    //atype c12=0;
    //atype c11=1.0;
    get_force_impr(U,E,eps*c11,eps*c12); // force for improved action:  
                                            // 2.0 * ( c11*\sum_{unoriented plaquettes P} ReTr(P) 
                                            //        + c12*\sum_{unoriented 1x2-rectangles R}  ReTr(R) )

    //get_force_clover(U,E,eps); // force for clover action (not useful by itself)

    // 2/9*BP+7/9*Wilson (cancels O(a^4) terms)
    //get_force_bp(U,E,eps*2.0/9.0);
    //get_force_wplaq_add(U,E,eps*7.0/9.0);
};

template <typename group,typename atype=hila::arithmetic_type<group>>
atype measure_gf_s(const GaugeField<group>& U) {

    //atype res=measure_s_wplaq(U);

    //atype res=measure_s_bp(U);
    
    // DBW2 action parameters:
    //atype c12=-1.4088;     // rectangle weight
    //atype c11=1.0-8.0*c12; // plaquette weight

    // Iwasaki action parameters:
    atype c12=-0.331;     // rectangle weight
    atype c11=1.0-8.0*c12; // plaquette weight

    // LW action parameters:
    //atype c12=-1.0/12.0;     // rectangle weight
    //atype c11=1.0-8.0*c12; // plaquette weight

    // Wilson plaquette action parameters:
    //atype c12=0;
    //atype c11=1.0;

    atype res=measure_s_impr(U,c11,c12);

    //atype res=measure_s_clover(U);

    // 2/9*BP+7/9*Wilson (cancels O(a^4) terms)
    //atype res=measure_s_bp(U,E,eps*2.0/9.0)+measure_s_wplaq(U,E,eps*7.0/9.0);

    return res;
}

template <typename group,typename atype=hila::arithmetic_type<group>>
atype measure_dE_wplaq_dt(const GaugeField<group>& U) {
    Reduction<atype> declov=0;
    declov.allreduce(false).delayed(true);
    VectorField<Algebra<group>> K,Kc;
    get_gf_force(U,K);
    get_force_wplaq(U,Kc,2.0);
    foralldir(d) {
        onsites(ALL) {
            declov+=Kc[d][X].dot(K[d][X]);
        }
    }
    return declov.value();
}

template <typename group,typename atype=hila::arithmetic_type<group>>
atype measure_dE_clov_dt(const GaugeField<group>& U) {
    Reduction<atype> declov=0;
    declov.allreduce(false).delayed(true);
    VectorField<Algebra<group>> K,Kc;
    get_gf_force(U,K);
    get_force_clover(U,Kc,2.0);
    foralldir(d) {
        onsites(ALL) {
            declov+=Kc[d][X].dot(K[d][X]);
        }
    }
    return declov.value();
}

template <typename group,typename atype=hila::arithmetic_type<group>>
void measure_gradient_flow_stuff(const GaugeField<group>& V,atype flow_l,atype t_step) {
    // perform measurements on flowed gauge configuration V at flow scale flow_l
    // [t_step is the flow time integration step size used in last gradient flow step]
    // and print results in formatted form to standard output
    static bool first=true;
    if(first) {
        // print legend for flow measurement output
        hila::out0<<"LWFLMEAS  l(ambda)        S-flow        S-plaq        E_plaq    dE_plaq/dl         E_clv     dE_clv/dl     Qtopo_clv         E_log     Qtopo_log   [t step size]\n";
        first=false;
    }
    atype slocal=measure_gf_s(V)/(lattice.volume()*NDIM*(NDIM-1)/2); // average action per plaquette 
    atype plaq=measure_s_wplaq(V)/(lattice.volume()*NDIM*(NDIM-1)/2); // average wilson plaquette action
    atype eplaq=plaq*NDIM*(NDIM-1)*group::size(); // naive energy density (based on wilson plaquette action)

    // average energy density and toplogical charge from 
    // symmetric log definition of field strength tensor :
    atype qtopolog,elog;
    measure_topo_charge_and_energy_log(V,qtopolog,elog);
    elog/=lattice.volume();

    // average energy density and toplogical charge from 
    // clover definition of field strength tensor :
    atype qtopocl,ecl;
    measure_topo_charge_and_energy_clover(V,qtopocl,ecl);
    ecl/=lattice.volume();

    // derivative of plaquette energy density w.r.t. to flow time :
    atype deplaqdt=measure_dE_wplaq_dt(V)/lattice.volume();

    // derivative of clover energy density w.r.t. to flow time :
    atype declovdt=measure_dE_clov_dt(V)/lattice.volume();

    // print formatted results to standard output :
    hila::out0<<string_format("WFLMEAS  % 9.3f % 0.6e % 0.6e % 0.6e % 0.6e % 0.6e % 0.6e % 0.6e % 0.6e % 0.6e       [%0.5f]",flow_l,slocal,plaq,eplaq,0.25*flow_l*deplaqdt,ecl,0.25*flow_l*declovdt,qtopocl,elog,qtopolog,t_step)<<'\n';
}

template <typename group,typename atype=hila::arithmetic_type<group>>
atype do_gradient_flow_adapt(GaugeField<group>& V,atype l_start,atype l_end,atype atol=1.0e-7,atype rtol=1.0e-7,atype tstep=0.001) {
    // wilson flow integration from flow scale l_start to l_end using 3rd order
    // 3-step Runge-Kutta (RK3) from arXiv:1006.4518 (cf. appendix C of
    // arXiv:2101.05320 for derivation of this Runge-Kutta method)
    // and embedded RK2 for adaptive step size

    // translate flow scale interval [l_start,l_end] to corresponding
    // flow time interval [t,tmax] :
    atype t=l_start*l_start/8.0;
    atype tmax=l_end*l_end/8.0;
    atype lstab=0.095; // stability limit
    atype step=min(min(tstep,0.51*(tmax-t)),lstab);  //initial step size

    VectorField<Algebra<group>> k1,k2;
    GaugeField<group> V2,V0;
    Field<atype> reldiff;
    atype maxreldiff,maxstep;

    // RK3 coefficients from arXiv:1006.4518 :
    // correspond to standard RK3 with Butcher-tableau 
    // (cf. arXiv:2101.05320, Appendix C)
    //  0  |   0     0     0
    //  #  |  1/4    0     0
    //  #  | -2/9   8/9    0
    // -------------------------
    //     |  1/4    0    3/4
    //
    atype a11=0.25;
    atype a21=-17.0/36.0,a22=8.0/9.0;
    atype a33=0.75;

    // RK2 coefficients :
    // cf. Alg. 6 and Eqn. (13)-(14) in arXiv:2101.05320 to see
    // how these are obtained from standard RK2 with Butcher-tableau
    //  0  |   0     0
    //  #  |  1/4    0
    // -----------------
    //     |  -1     2
    //
    atype b21=-1.25,b22=2.0;

    V0=V;
    bool stop=false;

    while(t<tmax&&!stop) {

        if(t+step>=tmax) {
            step=tmax-t;
            stop=true;
        } else {
            if(t+2.0*step>=tmax) {
                step=0.51*(tmax-t);
            }
        }

        get_gf_force(V,k1);
        foralldir(d) onsites(ALL) {
            // first steps of RK3 and RK2 are the same :
            V[d][X]=chexp(k1[d][X]*(step*a11))*V[d][X];
        }

        get_gf_force(V,k2);
        foralldir(d) onsites(ALL) {
            // second step of RK2 :
            V2[d][X]=chexp(k2[d][X]*(step*b22)+k1[d][X]*(step*b21))*V[d][X];

            // second step of RK3 :
            k2[d][X]=k2[d][X]*(step*a22)+k1[d][X]*(step*a21);
            V[d][X]=chexp(k2[d][X])*V[d][X];
        }

        get_gf_force(V,k1);
        foralldir(d) onsites(ALL) {
            // third step of RK3 :
            V[d][X]=chexp(k1[d][X]*(step*a33)-k2[d][X])*V[d][X];
        }

        // determine maximum difference between RK3 and RK2, 
        // relative to desired accuracy :
        maxreldiff=0;
        foralldir(d) {
            reldiff[ALL]=(V[d][X]-V2[d][X]).norm()/(atol+rtol*V[d][X].norm());
            maxreldiff=max(maxreldiff,reldiff.max());
        }
        maxreldiff/=(atype)(group::size()*group::size());

        // max. allowed step size to achieve desired accuracy :
        maxstep=min(step/pow(maxreldiff,1.0/3.0),1.0);
        if(step>maxstep) {
            // repeat current iteration if step size was larger than maxstep
            V=V0;
        } else {
            // proceed to next iteration
            t+=step;
            V0=V;
        }
        // adjust step size for next iteration to better match accuracy goal :
        step=min(0.9*maxstep,lstab);
    }

    V.reunitarize_gauge();
    return step;
}

#endif