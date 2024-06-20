/** @file wilson_flow.h */

#ifndef WILSON_FLOW_H_
#define WILSON_FLOW_H_

#include "hila.h"

/*
* definition of get_wf_force() will be needed, e.g.: 
template <typename group>
void get_wf_force(const GaugeField<group>& U,VectorField<Algebra<group>>& E) {
    // force computation routine to be used for wilson flow
    //get_force(U,E); // force for Wilson plaquette action
    //get_force_bp(U,E); // force for bulk-preventing action

    // DBW2 action parameters:
    //ftype c12=-1.4088;     // rectangle weight
    //ftype c11=1.0-8.0*c12; // plaquette weight

    // Iwasaki action parameters:
    ftype c12=-0.331;     // rectangle weight
    ftype c11=1.0-8.0*c12; // plaquette weight
    
    // LW action parameters:
    //ftype c12=-1.0/12.0;     // rectangle weight
    //ftype c11=1.0-8.0*c12; // plaquette weight

    // Wilson plaquette action parameters:
    //ftype c12=0;
    //ftype c11=1.0;

    get_force_impr_f2(U,E,c11,c12); // force for improved action  c11*\sum_{unoriented plaquettes P} ReTr(P) + c12*\sum_{unoriented 1x2-rectangles R}  ReTr(R) 
};
*/

template <typename group,typename atype=hila::arithmetic_type<group>>
atype do_wilson_flow_adapt(GaugeField<group>& V,atype l_start,atype l_end,void (*get_wf_force)(const GaugeField<group>&,VectorField<Algebra<group>>&),atype atol=1.0e-7,atype rtol=1.0e-7,atype tstep=0.001) {
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

        get_wf_force(V,k1);
        foralldir(d) onsites(ALL) {
            // first steps of RK3 and RK2 are the same :
            V[d][X]=chexp(k1[d][X]*(step*a11))*V[d][X];
        }

        get_wf_force(V,k2);
        foralldir(d) onsites(ALL) {
            // second step of RK2 :
            V2[d][X]=chexp(k2[d][X]*(step*b22)+k1[d][X]*(step*b21))*V[d][X];

            // second step of RK3 :
            k2[d][X]=k2[d][X]*(step*a22)+k1[d][X]*(step*a21);
            V[d][X]=chexp(k2[d][X])*V[d][X];
        }

        get_wf_force(V,k1);
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