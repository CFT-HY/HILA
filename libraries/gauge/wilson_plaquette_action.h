/** @file wilson_plaquette_action.h */

#ifndef WILSON_PLAQUETTE_ACTION_H_
#define WILSON_PLAQUETTE_ACTION_H_

#include "hila.h"

// functions for Wilson's plaquette action -S_{plaq}=\beta/N * \sum_{plaq} ReTr(plaq)

template <typename T>
void rstaplesum(const GaugeField<T> &U, out_only Field<T> &staples, Direction d1) {

    Field<T> lower;
    bool first = true;
    foralldir(d2) if (d2 != d1) {

        // anticipate that these are needed
        // not really necessary, but may be faster
        U[d2].start_gather(d1, ALL);
        U[d1].start_gather(d2, ALL);

        // calculate first lower 'u' of the staple sum
        onsites(ALL) {
            lower[X] = (U[d1][X] * U[d2][X + d1]).dagger() * U[d2][X];
        }
        lower.start_gather(-d2, ALL);

        // calculate then the upper 'n'
        if (first) {
            onsites(ALL) {
                staples[X] = U[d2][X + d1] * (U[d2][X] * U[d1][X + d2]).dagger();
            }
            first = false;
        } else {
            onsites(ALL) {
                staples[X] += U[d2][X + d1] * (U[d2][X] * U[d1][X + d2]).dagger();
            }
        }

        // add lower
        onsites(ALL) {
            staples[X] += lower[X - d2];
        }
        
    }
}

template <typename group>
double measure_s_wplaq(const GaugeField<group> &U) {
    // measure the total Wilson plaquette action for the gauge field U
    Reduction<double> plaq = 0;
    plaq.allreduce(false).delayed(true);
    foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {
        U[dir2].start_gather(dir1, ALL);
        U[dir1].start_gather(dir2, ALL);
        onsites(ALL) {
            plaq += 1.0 - real(trace(U[dir1][X] * U[dir2][X + dir1] *
                                     (U[dir2][X] * U[dir1][X + dir2]).dagger())) /
                              group::size();
        }
    }
    return plaq.value();
}

template <typename group, typename atype = hila::arithmetic_type<group>>
double measure_s_wplaq(const GaugeField<group> &U, out_only atype &max_plaq) {
    // measure the total and maximal Wilson plaquette action for the gauge field U
    Reduction<double> plaq = 0;
    max_plaq = -1.0;
    Field<atype> P;
    plaq.allreduce(false).delayed(true);
    foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {
        U[dir2].start_gather(dir1, ALL);
        U[dir1].start_gather(dir2, ALL);
        onsites(ALL) {
            P[X] = 1.0 - real(trace(U[dir1][X] * U[dir2][X + dir1] *
                                    (U[dir2][X] * U[dir1][X + dir2]).dagger())) /
                             group::size();
            plaq += (double)P[X];
        }
        atype tmax_plaq = P.max();
        if(tmax_plaq>max_plaq) {
            max_plaq = tmax_plaq;
        }
    }
    return plaq.value();
}

template <typename group, typename atype = hila::arithmetic_type<group>>
void get_force_wplaq_add(const GaugeField<group> &U, VectorField<Algebra<group>> &K,
                         atype eps = 1.0) {
    // compute the force for the plaquette action and write it to K

    Field<group> staple;
    foralldir(d) {
        rstaplesum(U, staple, d);
        onsites(ALL) {
            K[d][X] -= (U[d][X] * staple[X]).project_to_algebra_scaled(eps);
        }
    }

    /*
    Field<group> tP,tP1,tP2;
    foralldir(d1) {
        foralldir(d2) if(d1<d2) {
            U[d2].start_gather(d1,ALL);
            U[d1].start_gather(d2,ALL);
            onsites(ALL) {
                tP[X]=U[d1][X]*U[d2][X+d1]*(U[d2][X]*U[d1][X+d2]).dagger();
                tP1[X]=(tP[X]*U[d2][X]).dagger()*U[d2][X]; // parallel transport tP[X].dagger() to X+d2 
                tP2[X]=U[d1][X].dagger()*tP[X]*U[d1][X]; // parallel transport tP[X] to X+d1
            }
            tP1.start_gather(-d2,ALL);
            tP2.start_gather(-d1,ALL);
            onsites(ALL) {
                K[d1][X]-=(tP1[X-d2]+tP[X]).project_to_algebra_scaled(eps);
                K[d2][X]-=(tP2[X-d1]-tP[X]).project_to_algebra_scaled(eps);
            }
        }
    }
    */
}

template <typename group, typename atype = hila::arithmetic_type<group>>
void get_force_wplaq(const GaugeField<group> &U, out_only VectorField<Algebra<group>> &K,
                     atype eps = 1.0) {
    // compute the force for the plaquette action and write it to K
    Field<group> staple;
    foralldir(d) {
        rstaplesum(U, staple, d);
        onsites(ALL) {
            K[d][X] = (U[d][X] * staple[X]).project_to_algebra_scaled(-eps);
        }
    }
}


#endif