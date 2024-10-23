/** @file log_plaquette_action.h */

#ifndef LOG_PLAQUETTE_ACTION_H_
#define LOG_PLAQUETTE_ACTION_H_

#include "hila.h"
#include "gauge/wilson_line_and_force.h"

// functions for log clover gauge action

template <typename group, typename atype = hila::arithmetic_type<group>>
void get_log_plaq_mat(const GaugeField<group> &U, Direction dir1, Direction dir2,
                      out_only Field<group> &C) {
    // sets C=log(U) where U is the plaquette variable panned by dir1 and dir2
    U[dir2].start_gather(dir1, ALL);
    U[dir1].start_gather(dir2, ALL);

    onsites(ALL) {
        // log of dir1-dir2-plaquette that starts and ends at X; corresponds to F[dir1][dir2]
        // at center location X+dir1/2+dir2/2 of plaquette:
        C[X] = log((U[dir1][X] * U[dir2][X + dir1] * (U[dir2][X] * U[dir1][X + dir2]).dagger()))
                   .expand();
    }
}

template <typename group>
double measure_s_log_plaq(const GaugeField<group> &U) {
    // measure the log action for dir1<dir2
    // (just to have same normalization as with plaquette action)
    Reduction<double> stot = 0;
    stot.allreduce(false).delayed(true);
    Field<group> C;
    foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {
        get_log_plaq_mat(U, dir1, dir2, C);
        onsites(ALL) {
            stot += 0.5 * C[X].squarenorm();
        }
    }
    return stot.value() / group::size();
}

template <typename group, typename atype = hila::arithmetic_type<group>>
void get_force_log_plaq_add(const GaugeField<group> &U, VectorField<Algebra<group>> &K,
                       atype eps = 1.0) {
    // compute gauge force for clover action and add result to K
    Field<group> C;

    foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {
        // get matrix log of the plaquette matrix in (dir1,dir2)-plane
        get_log_plaq_mat(U, dir1, dir2, C);

        // define the path around the plaquette
        std::vector<Direction> path = {dir1, dir2, -dir1, -dir2};

        // the force on the link (Y,dir) from the log of plaquette C with foot point X is:
        // F[dir][Y] = -d_{Y,dir}(-0.5*Tr(C*C))
        //           = 0.5*Tr((d_{Y,dir}C)*C) + 0.5*Tr(C*(d_{Y,dir}C))
        //           = Tr((d_{Y,dir}C)*C)
        // with: d_{Y,dir}log(P[X]) = (d_{Y,dir}P[X])*P[X].dagger()
        // furthermore: (d^a_{X,dir[0]}P[X])*P[X].dagger() = T^a
        //              (d^a_{X+dir[0],dir[1]}P[X])*P[X].dagger() = 
        //                                             U[dir[0]][X].dagger() * T^a * U[dir[0]][X]
        //              ...
        //              ...
        // We therefore need only C to compute the log-action force:
        get_wloop_force_from_wl_add(U, path, C, eps, K);
    }
}

template <typename group, typename atype = hila::arithmetic_type<group>>
void get_force_log_plaq(const GaugeField<group> &U, out_only VectorField<Algebra<group>> &K,
                        atype eps = 1.0) {
    // determine gauge force for clover action and store result in K
    foralldir(d1) {
        K[d1][ALL] = 0;
    }
    get_force_log_plaq_add(U, K, eps);
}


#endif