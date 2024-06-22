/** @file improved_action.h */

#ifndef IMPROVED_ACTION_H_
#define IMPROVED_ACTION_H_

#include "hila.h"
#include "gauge/wilson_line_and_force.h"

/* 
* parameters for for different improved actions:
// DBW2 action parameters:
//ftype c12=-1.4088;     // rectangle weight
//ftype c11=1.0-8.0*c12; // plaquette weight

// Iwasaki action parameters:
//ftype c12=-0.331;     // rectangle weight
//ftype c11=1.0-8.0*c12; // plaquette weight

// LW action parameters:
//ftype c12=-1.0/12.0;     // rectangle weight
//ftype c11=1.0-8.0*c12; // plaquette weight

// Wilson plaquette action parameters:
//ftype c12=0;
//ftype c11=1.0;
*/


// functions for improved action -S_{impr} = \p.beta/N*( c11*\sum_{plaquettes P} ReTr(P) + c12*\sum_{1x2-rectangles R}  ReTr(R) )
/*
template <typename group,typename atype=hila::arithmetic_type<group>>
void get_force_impr_add(const GaugeField<group>& U,VectorField<Algebra<group>>& K,atype c11,atype c12) {
    // compute the force for the improved action -S_{impr}=\beta/N*(c11*ReTr(plaq)+c12*ReTr(rect))
    // and write it to K

    // plaquette part:
    foralldir(dir1) foralldir(dir2) if(dir1<dir2) {
        Direction path[4]={dir1,dir2,-dir1,-dir2};
        get_wloop_force_add(U,path,c11,K);
    }

    if(c12!=0) {
        // rectangle part:
        foralldir(dir1) foralldir(dir2) if(dir1!=dir2) {
            Direction path[6]={dir1,dir2,dir2,-dir1,-dir2,-dir2};
            get_wloop_force_add(U,path,c12,K);
        }
    }
}
*/

template <typename group,typename atype=hila::arithmetic_type<group>>
atype measure_s_impr(const GaugeField<group>& U, atype c11, atype c12) {
    // measure the improved action for dir1<dir2
    Reduction<hila::arithmetic_type<group>> stot=0;
    stot.allreduce(false).delayed(true);
    Field<group> tP;
    foralldir(dir1) foralldir(dir2) if(dir1<dir2) {
        U[dir2].start_gather(dir1,ALL);
        U[dir1].start_gather(dir2,ALL);
        onsites(ALL) {
            //plaquettes:
            tP[X]=U[dir1][X]*U[dir2][X+dir1]*(U[dir2][X]*U[dir1][X+dir2]).dagger();

            //plaquette part:
            stot+=c11*(1.0-real(trace(tP[X]))/group::size());
        }

        tP.start_gather(dir1,ALL);
        tP.start_gather(dir2,ALL);

        onsites(ALL) {
            //2x1-rectangle part:
            stot+=c12*(1.0-real(trace(U[dir1][X]*tP[X+dir1]*U[dir1][X].dagger()*tP[X]))/group::size());
            //1x2-rectangle part:
            stot+=c12*(1.0-real(trace(tP[X]*U[dir2][X]*tP[X+dir2]*U[dir2][X].dagger()))/group::size());
        }
    }
    return stot.value();
}

template <typename group,typename atype=hila::arithmetic_type<group>>
void get_force_impr_add(const GaugeField<group>& U,VectorField<Algebra<group>>& K,atype c11,atype c12) {
    // compute the force for the improved action -S_{impr}=\beta/N*(c11*ReTr(plaq)+c12*ReTr(rect))
    // in an even faster way and add the results to K
    Field<group> ustap;
    Field<group> lstap;
    Field<group> tstap;
    bool first=true;
    foralldir(dir1) {
        first=true;
        foralldir(dir2) if(dir1!=dir2) {
            U[dir2].start_gather(dir1,ALL);
            U[dir1].start_gather(dir2,ALL);

            // get upper (dir1,dir2) and lower (dir1,-dir2) staples
            onsites(ALL) {
                lstap[X]=(U[dir1][X]*U[dir2][X+dir1]).dagger()*U[dir2][X];
                ustap[X]=U[dir2][X+dir1]*(U[dir2][X]*U[dir1][X+dir2]).dagger();
            }

            lstap.start_gather(-dir2,ALL);

            // sum the staples for dir1-link
            if(first) {
                onsites(ALL) {
                    tstap[X]=ustap[X];
                    tstap[X]+=lstap[X-dir2];
                }
                first=false;
            } else {
                onsites(ALL) {
                    tstap[X]+=ustap[X];
                    tstap[X]+=lstap[X-dir2];
                }
            }

            // compute rectangle contribution to force (if any) and add it to K
            if(c12!=0) {
                // compose rectangle and store it in ustap
                onsites(ALL) ustap[X]=lstap[X-dir2].dagger()*ustap[X];

                // corresponding path
                Direction path[6]={-dir2,dir1,dir2,dir2,-dir1,-dir2};

                // compute rectangle force and add it to K
                get_wloop_force_from_wl_add(U,path,6,ustap,c12,K);
            }
        }
        // compute plaquette contribution to force and add it to K
        onsites(ALL) {
            K[dir1][X]-=(U[dir1][X]*tstap[X]).project_to_algebra_scaled(c11);
        }
    }
}

template <typename group,typename atype=hila::arithmetic_type<group>>
void get_force_impr(const GaugeField<group>& U,VectorField<Algebra<group>>& K,atype c11,atype c12) {
    // compute the force for the improved action -S_{impr}=\beta/N*(c11*ReTr(plaq)+c12*ReTr(rect))
    // in an even faster way and write the results to K
    foralldir(d1) {
        K[d1][ALL]=0;
    }
    get_force_impr_add(U,K,c11,c12);
}


#endif