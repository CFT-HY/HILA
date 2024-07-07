/** @file log_clover_action.h */

#ifndef LOG_CLOVER_ACTION_H_
#define LOG_CLOVER_ACTION_H_

#include "hila.h"
#include "gauge/wilson_line_and_force.h"

// functions for log clover gauge action

template <typename group, typename atype = hila::arithmetic_type<group>>
void get_log_clover_mat(const GaugeField<group> &U, Direction dir1, Direction dir2,
                        Field<group> &Csum) {
    // sets Csum=(C[0]+C[1]+C[2]+C[3])/4 where C[0],C[1],C[2],C[3] are the logs of the clover-leaf
    // matrices in counter-clockwise order

    U[dir2].start_gather(dir1, ALL);
    U[dir1].start_gather(dir2, ALL);

    Field<group> tF;

    onsites(ALL) {
        // log of dir1-dir2-plaquette that starts and ends at X; corresponds to F[dir1][dir2]
        // at center location X+dir1/2+dir2/2 of plaquette:
        Csum[X] = log((U[dir1][X] * U[dir2][X + dir1] * (U[dir2][X] * U[dir1][X + dir2]).dagger()))
                     .expand();
        // parallel transport to X+dir1
        tF[X] = U[dir1][X].dagger() * Csum[X] * U[dir1][X];
    }
    tF.start_gather(-dir1, ALL);
    // add parallel transported leaf-matrix from X-dir1 to the leaf-matrix at X
    onsites(ALL) {
        Csum[X] += tF[X - dir1];
    }

    onsites(ALL) {
        // parallel transport two-leaves sum from X to X+dir2
        tF[X] = U[dir2][X].dagger() * Csum[X] * U[dir2][X];
    }
    tF.start_gather(-dir2, ALL);
    // add parallel transported two-leaves sum from X-dir2 to two-leaves sum at X and divide by 4
    onsites(ALL) {
        Csum[X] += tF[X - dir2];
        Csum[X] *= 0.25;
    }
}

template <typename group>
double measure_s_log(const GaugeField<group> &U) {
    // measure the log action for dir1<dir2
    // (just to have same normalization as with plaquette action)
    Reduction<double> stot = 0;
    stot.allreduce(false).delayed(true);
    Field<group> Csum;
    foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {
        get_log_clover_mat(U, dir1, dir2, Csum);
        onsites(ALL) {
            stot += 0.5 * Csum[X].squarenorm();
        }
    }
    return stot.value() / group::size();
}

template <typename group, typename atype = hila::arithmetic_type<group>>
void get_force_log_add(const GaugeField<group> &U, VectorField<Algebra<group>> &K,
                       atype eps = 1.0) {
    // compute gauge force for clover action and add result to K
    Field<group> Csum;

    foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {
        // get the 4 clover leave matrices C[][X] (ordered counter-clockwise in (dir1,dir2)-plane)
        // and their anti-hermitian traceless sum Csum[X]
        get_log_clover_mat(U, dir1, dir2, Csum);

        // define for each clover leaf matrix the correspondin path of gauge links
        std::vector<Direction> paths[4] = {{dir1, dir2, -dir1, -dir2},
                                           {dir2, -dir1, -dir2, dir1},
                                           {-dir1, -dir2, dir1, dir2},
                                           {-dir2, dir1, dir2, -dir1}};

        // the force on the link (Y,dir) from the clover term with foot point X is:
        // F[dir][Y] = -d_{Y,dir}(-0.5*Tr(Csum*Csum))
        //           = 0.5*Tr((d_{Y,dir}Csum)*Csum) + 0.5*Tr(Csum*(d_{Y,dir}Csum))
        //           = Tr((d_{Y,dir}Csum)*Csum)
        //           = 0.25*\sum_k Tr((d_{Y,dir}log(C[k][X]))*Csum[X])
        // with: d_{Y,dir}log(C[k][X]) = (d_{Y,dir}C[k][X])*C[k][X].dagger()
        //       where C[k][X] is the k-th clover leaf matrix,
        // furthermore: (d^a_{X,dir[0]}C[0][X])*C[0][X].dagger() = T^a
        //              (d^a_{X+dir[0],dir[1]}C[0][X])*C[0][X].dagger() = 
        //                                             U[dir[0]][X].dagger() * T^a * U[dir[0]][X]
        //              ...
        //              ...
        // We therefore need only Csum/4 to compute the log-action force:
        for (int k = 0; k < 4; ++k) {
            // compute the gauge force from the first 4 links of all the Wilson loops that are
            // summed in tC[k] (is valid to do so since first 4 links are the same for all Wilson
            // loops and projection on Lie-algebra is linear)
            get_wloop_force_from_wl_add(U, paths[k], Csum, 0.25*eps, K);
        }
    }
}

template <typename group, typename atype = hila::arithmetic_type<group>>
void get_force_log(const GaugeField<group> &U, VectorField<Algebra<group>> &K, atype eps = 1.0) {
    // determine gauge force for clover action and store result in K
    foralldir(d1) {
        K[d1][ALL] = 0;
    }
    get_force_log_add(U, K, eps);
}


#endif