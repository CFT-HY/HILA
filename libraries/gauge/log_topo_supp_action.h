/** @file log_topo_supp_action.h */

#ifndef LOG_TOPO_SUPP_ACTION_H_
#define LOG_TOPO_SUPP_ACTION_H_

#include "hila.h"
// #include "gauge/log_clover_action.h"

// functions for topo.-sup. gauge action

template <typename group, typename atype = hila::arithmetic_type<group>>
void get_log_p_and_c_mat(const GaugeField<group> &U, Direction dir1, Direction dir2,
                         out_only Field<group> &C, out_only Field<group> &Csum) {
    // sets Csum=(C[0]+C[1]+C[2]+C[3])/4 where C[0],C[1],C[2],C[3] are the logs of the clover-leaf
    // matrices in counter-clockwise order

    U[dir2].start_gather(dir1, ALL);
    U[dir1].start_gather(dir2, ALL);

    Field<group> tF;

    onsites(ALL) {
        // log of dir1-dir2-plaquette that starts and ends at X; corresponds to F[dir1][dir2]
        // at center location X+dir1/2+dir2/2 of plaquette:
        C[X] = log((U[dir1][X] * U[dir2][X + dir1] * (U[dir2][X] * U[dir1][X + dir2]).dagger()))
                      .expand();
        // parallel transport to X+dir1
        tF[X] = U[dir1][X].dagger() * C[X] * U[dir1][X];
    }
    tF.start_gather(-dir1, ALL);
    // set Csum to leaf-matrix at X
    onsites(ALL) {
        Csum[X] = C[X];
    }
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
double measure_s_topo_supp(const GaugeField<group> &U) {
    // double teps = 1.0 / (4.0 * M_PI * M_PI);
    
    Field<group> F[6];
    //Field<group> Fsum[6];

    int k = 0;
    foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {
        // get the 4 clover leave matrices C[][X] (ordered counter-clockwise in (dir1,dir2)-plane)
        // and their anti-hermitian traceless sum Csum[X]
        get_log_plaq_mat(U, dir1, dir2, F[k]);

        ++k;
    }

    int kmap[6] = {5, 4, 3, 2, 1, 0};

    //int kmap[6] = {0, 1, 2, 3, 4, 5};

    //int psign[6] = {1, -1, 1, 1, -1, 1};

    Reduction<double> stot = 0;
    stot.allreduce(false).delayed(true);
    for (k = 0; k < 3; ++k) {
        onsites(ALL) {
            double ttstot = -real(mul_trace(F[k][X], F[kmap[k]][X]));
            stot += ttstot * ttstot;
        }
    }

    return stot.value() / group::size();
}

template <typename group, typename atype = hila::arithmetic_type<group>>
void get_force_topo_supp_add(const GaugeField<group> &U, VectorField<Algebra<group>> &K,
                             atype eps = 1.0) {
    // compute gauge force for clover action and add result to K
    // atype qtoponf = 1.0 / (4.0 * M_PI * M_PI);
    atype teps = 2.0 * eps;

    Field<group> F[6];
    //Field<group> Fsum[6];
    
    int kmap[6] = {5, 4, 3, 2, 1, 0};

    //int kmap[6] = {0, 1, 2, 3, 4, 5};

    //int psign[6] = {1, -1, 1, 1, -1, 1};
    Field<atype> Qdens[3];


    int k = 0;
    foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {
        // get the 4 clover leave matrices C[][X] (ordered counter-clockwise in (dir1,dir2)-plane)
        // and their anti-hermitian traceless sum Csum[X]
        get_log_plaq_mat(U, dir1, dir2, F[k]);

        ++k;
    }

    
    for (k = 0; k < 3; ++k) {
        onsites(ALL) {
            Qdens[k][X] = -real(mul_trace(F[k][X], F[kmap[k]][X]));
        }
    }

    for (k = 0; k < 3; ++k) {
        onsites(ALL) F[k][X] *= Qdens[k][X];
        onsites(ALL) F[kmap[k]][X] *= Qdens[k][X];
    }
    


    k = 0;
    foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {

        std::vector<Direction> path = {dir1, dir2, -dir1, -dir2};

        get_wloop_force_from_wl_add(U, path, F[kmap[k]], teps, K);

        /*
        // define for each clover leave matrix the correspondin path of gauge links
        std::vector<Direction>
            paths[4] = {{dir1, dir2, -dir1, -dir2},
                        {dir2, -dir1, -dir2, dir1},
                        {-dir1, -dir2, dir1, dir2},
                        {-dir2, dir1, dir2, -dir1}};

        for (int m = 0; m < 4; ++m) {
            get_wloop_force_from_wl_add(U, paths[m], F[kmap[k]], 0.25 * teps, K);
        }
        */
        ++k;
    }
}

template <typename group, typename atype = hila::arithmetic_type<group>>
void get_force_topo_supp(const GaugeField<group> &U, out_only VectorField<Algebra<group>> &K,
                         atype eps = 1.0) {
    // determine gauge force for clover action and store result in K
    foralldir(d1) {
        K[d1][ALL] = 0;
    }
    get_force_topo_supp_add(U, K, eps);
}


#endif