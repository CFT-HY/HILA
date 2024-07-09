/** @file clover_action.h */

#ifndef CLOVER_ACTION_H_
#define CLOVER_ACTION_H_

#include "hila.h"
#include "gauge/wilson_line_and_force.h"

// functions for clover gauge action

template <typename group, typename atype = hila::arithmetic_type<group>>
void get_clover_leaves(const GaugeField<group> &U, Direction d1, Direction d2,
                       Field<group>(out_only &C)[4], out_only Field<group> &Csum) {
    // assignes the clover leaf matrices in counter-clockwise order to C[0],C[1],C[2],C[3] and sets
    // Csum to the full Clover matrix, i.e. to the anti-hermitian tracless part of
    // (C[0]+C[1]+C[2]+C[3])/4

    U[d2].start_gather(d1, ALL);
    U[d1].start_gather(d2, ALL);
    Field<group> tF;

    onsites(ALL) {
        // d1-d2-plaquette that starts and ends at X
        C[0][X] = U[d1][X] * U[d2][X + d1] * (U[d2][X] * U[d1][X + d2]).dagger();
        C[0][X] *= 0.25;

        // parallel transport C[0][X] to X+d1 (use C[2][X] as temp. storage)
        C[2][X] = U[d1][X].dagger() * C[0][X] * U[d1][X];
    }

    C[2].start_gather(-d1, ALL);

    onsites(ALL) {
        // parallel transport C[0][X] to X+d2 (use Csum as temp. storage)
        Csum[X] = U[d2][X].dagger() * C[0][X] * U[d2][X];
    }

    Csum.start_gather(-d2, ALL);

    onsites(ALL) {
        C[1][X] = C[2][X - d1];

        // parallel transport C[1][X] to X+d2 (use tF as temp. storage)
        tF[X] = U[d2][X].dagger() * C[1][X] * U[d2][X];
    }

    tF.start_gather(-d2, ALL);

    onsites(ALL) {
        C[3][X] = Csum[X - d2];
    }

    onsites(ALL) {
        C[2][X] = tF[X - d2];

        Csum[X] = C[0][X];
        Csum[X] += C[1][X];
        Csum[X] += C[2][X];
        Csum[X] += C[3][X];

        // Lie-algebra project Csum[X]
        Csum[X] -= Csum[X].dagger();
        Csum[X] *= 0.5;
        Csum[X] -= trace(Csum[X]) / group::size();
    }
}

template <typename group>
double measure_s_clover(const GaugeField<group> &U) {
    // measure the clover action for dir1<dir2
    // (just to have same normalization as with plaquette action)
    Reduction<double> stot = 0;
    stot.allreduce(false).delayed(true);
    Field<group> C[4];
    Field<group> Csum;
    foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {
        get_clover_leaves(U, dir1, dir2, C, Csum);
        onsites(ALL) {
            stot += 0.5 * Csum[X].squarenorm();
        }
    }
    return stot.value() / group::size();
}

template <typename group, typename atype = hila::arithmetic_type<group>>
void get_force_clover_add(const GaugeField<group> &U, VectorField<Algebra<group>> &K,
                          atype eps = 1.0) {
    // compute gauge force for clover action and add result to K
    Field<group> C[4];
    Field<group> Csum;

    Field<group> tC;
    foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {
        // get the 4 clover leaf matrices C[][X] (ordered counter-clockwise in (dir1,dir2)-plane)
        // and their anti-hermitian traceless sum Csum[X]
        get_clover_leaves(U, dir1, dir2, C, Csum);

        // define for each clover leaf matrix the correspondin path of gauge links
        std::vector<Direction> paths[4] = {{dir1, dir2, -dir1, -dir2},
                                           {dir2, -dir1, -dir2, dir1},
                                           {-dir1, -dir2, dir1, dir2},
                                           {-dir2, dir1, dir2, -dir1}};

        // the force on the link (Y,dir) from the clover term with foot point X is:
        // F[dir][Y] = d_{Y,dir}(0.5*Tr(Csum*Csum))
        //           = 0.5*Tr((d_{Y,dir}Csum)*Csum) + 0.5*Tr(Csum*(d_{Y,dir}Csum))
        //           = Tr((d_{Y,dir}Csum)*Csum)
        //           = 0.5*\sum_k (Tr((d_{Y,dir}C[k][X])*Csum[X]) -
        //           Tr((d_{Y,dir}C[k][X].dagger())*Csum[X])) = \sum_k
        //           Tr((d_{Y,dir}C[k][X])*Csum[X])
        for (int k = 0; k < 4; ++k) {
            // multiply the clover leave matrix C[k][X] by Csum[X] --> yields (weighted) sum
            // of Wilson loops with foot point X and whose first 4 links are always those of C[k][X]
            onsites(ALL) {
                tC[X] = C[k][X] * Csum[X];
            }
            // compute the gauge force from the first 4 links of all the Wilson loops that are
            // summed in tC[k] (is valid to do so since first 4 links are the same for all Wilson
            // loops and projection on Lie-algebra is linear)
            get_wloop_force_from_wl_add(U, paths[k], tC, eps, K);
        }
    }
}

template <typename group, typename atype = hila::arithmetic_type<group>>
void get_force_clover(const GaugeField<group> &U, out_only VectorField<Algebra<group>> &K,
                      atype eps = 1.0) {
    // determine gauge force for clover action and store result in K
    foralldir(d1) {
        K[d1][ALL] = 0;
    }
    get_force_clover_add(U, K, eps);
}


#endif