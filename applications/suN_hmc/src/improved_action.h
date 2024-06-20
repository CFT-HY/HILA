/** @file improved_action.h */

#ifndef IMPROVED_ACTION_H_
#define IMPROVED_ACTION_H_

#include "hila.h"
#include "wilson_line_and_force.h"

// functions for improved action -S_{impr} = \p.beta/N*( c11*\sum_{plaquettes P} ReTr(P) + c12*\sum_{1x2-rectangles R}  ReTr(R) )

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

template <typename group,typename atype=hila::arithmetic_type<group>>
void get_force_impr(const GaugeField<group>& U,VectorField<Algebra<group>>& K,atype c11,atype c12) {
    // compute the force for the improved action -S_{impr}=\beta/N*(c11*ReTr(plaq)+c12*ReTr(rect))
    // and write it to K
    foralldir(d1) {
        K[d1][ALL]=0;
    }
    get_force_impr_add(U,K,c11,c12);
}

template <typename group,typename atype=hila::arithmetic_type<group>>
void update_E_impr(const GaugeField<group>& U,VectorField<Algebra<group>>& E,atype c11,atype c12,atype delta) {
    // compute the force for the improved action -S_{impr}=\beta/N*(c11*ReTr(plaq)+c12*ReTr(rect))
    // and use it to evolve the momentum field E
    auto eps=delta/group::size();
    get_force_impr_add(U,E,eps*c11,eps*c12);
}


template <typename group,typename atype=hila::arithmetic_type<group>>
void get_force_impr_f_add(const GaugeField<group>& U,VectorField<Algebra<group>>& K,atype c11,atype c12) {
    // compute the force for the improved action -S_{impr}=\beta/N*(c11*ReTr(plaq)+c12*ReTr(rect))
    // in a faster way and add the results to K
    Field<group> ustap;
    Field<group> lstap;
    foralldir(dir1) foralldir(dir2) if(dir1!=dir2) {
        U[dir2].start_gather(dir1,ALL);
        U[dir1].start_gather(dir2,ALL);
        // get upper (dir1,dir2) and lower (dir1,-dir2) staples
        onsites(ALL) {
            lstap[X]=(U[dir1][X]*U[dir2][X+dir1]).dagger()*U[dir2][X];
            ustap[X]=U[dir2][X+dir1]*(U[dir2][X]*U[dir1][X+dir2]).dagger();
        }

        lstap.start_gather(-dir2,ALL);

        // compute plaquette contribution to force and add it to K
        onsites(ALL) {
            K[dir1][X]-=(U[dir1][X]*(ustap[X]+lstap[X-dir2])).project_to_algebra_scaled(c11);
        }

        if(c12!=0) {
            // rectangle contribution to force
            Direction path[6]={-dir2,dir1,dir2,dir2,-dir1,-dir2};

            // compose rectangle and store it in ustap
            onsites(ALL) ustap[X]=lstap[X-dir2].dagger()*ustap[X];

            // compute rectangle force and add it to K
            get_wloop_force_add(U,path,ustap,c12,K);
        }
    }
}

template <typename group,typename atype=hila::arithmetic_type<group>>
void get_force_impr_f(const GaugeField<group>& U,VectorField<Algebra<group>>& K,atype c11,atype c12) {
    // compute the force for the improved action -S_{impr}=\beta/N*(c11*ReTr(plaq)+c12*ReTr(rect))
    // in a faster way and write the results to K
    foralldir(d1) {
        K[d1][ALL]=0;
    }
    get_force_impr_f_add(U,K,c11,c12);
}

template <typename group,typename atype=hila::arithmetic_type<group>>
void update_E_impr_f(const GaugeField<group>& U,VectorField<Algebra<group>>& E,atype c11,atype c12,atype delta) {
    // compute the force for the improved action -S_{impr}=\beta/N*(c11*ReTr(plaq)+c12*ReTr(rect))
    // in a faster way and write the results to E

    auto eps=delta/group::size();
    get_force_impr_f_add(U,E,eps*c11,eps*c12);
}

template <typename group,typename atype=hila::arithmetic_type<group>>
atype measure_s_impr_f2(const GaugeField<group>& U, atype c11, atype c12) {
    // measure the improved action for dir1<dir2
    // (just to have same normalization as with plaquette action)
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
void get_force_impr_f2_add(const GaugeField<group>& U,VectorField<Algebra<group>>& K,atype c11,atype c12) {
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

            // sum the dir1-staples
            if(first) {
                onsites(ALL) {
                    tstap[X]=ustap[X];
                    tstap[X]+=lstap[X-dir2];
                }
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
                get_wloop_force_add(U,path,ustap,c12,K);
            }
            first=false;
        }
        // compute plaquette contribution to force and add it to K
        onsites(ALL) {
            K[dir1][X]-=(U[dir1][X]*tstap[X]).project_to_algebra_scaled(c11);
        }
    }
}

template <typename group,typename atype=hila::arithmetic_type<group>>
void get_force_impr_f2(const GaugeField<group>& U,VectorField<Algebra<group>>& K,atype c11,atype c12) {
    // compute the force for the improved action -S_{impr}=\beta/N*(c11*ReTr(plaq)+c12*ReTr(rect))
    // in an even faster way and write the results to K
    foralldir(d1) {
        K[d1][ALL]=0;
    }
    get_force_impr_f2_add(U,K,c11,c12);
}


template <typename group,typename atype=hila::arithmetic_type<group>>
void update_E_impr_f2(const GaugeField<group>& U,VectorField<Algebra<group>>& E,atype c11,atype c12,atype delta) {
    // compute the force for the improved action -S_{impr}=\beta/N*(c11*ReTr(plaq)+c12*ReTr(rect))
    // in an even faster way and write the results to K

    auto eps=delta/group::size();
    get_force_impr_f2_add(U,E,eps*c11,eps*c12);
}


#endif