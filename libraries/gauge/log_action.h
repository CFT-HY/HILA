/** @file log_action.h */

#ifndef LOG_ACTION_H_
#define LOG_ACTION_H_

#include "hila.h"
#include "gauge/wilson_line_and_force.h"

// functions for clover gauge action

template <typename group,typename atype=hila::arithmetic_type<group>>
void get_log_clover_mat(const GaugeField<group>& U,Direction dir1,Direction dir2,Field<group>& Csum) {
    // sets Csum=(C[0]+C[1]+C[2]+C[3])/4 where C[0],C[1],C[2],C[3] are the logs of the clover-leave
    // matrices in counter-clockwise order 

    U[dir2].start_gather(dir1,ALL);
    U[dir1].start_gather(dir2,ALL);

    Field<group> tF0,tF1;

    onsites(ALL) {
        // log of dir1-dir2-plaquette that starts and ends at X; corresponds to F[dir1][dir2]
        // at center location X+dir1/2+dir2/2 of plaquette:
        tF0[X]=log((U[dir1][X]*U[dir2][X+dir1]*(U[dir2][X]*U[dir1][X+dir2]).dagger())).expand();
        // parallel transport to X+dir1
        tF1[X]=U[dir1][X].dagger()*tF0[X]*U[dir1][X];
    }

    tF1.start_gather(-dir1,ALL);
    onsites(ALL) {
        tF0[X]+=tF1[X-dir1];
    }

    U[dir2].start_gather(-dir2,ALL);
    tF0.start_gather(-dir2,ALL);
    onsites(ALL) {
        // get F[dir1][dir2] at X from average of the (parallel transported) F[dir1][dir2] from 
        // the centers of all dir1-dir2-plaquettes that touch X :
        Csum[X]=(tF0[X]+U[dir2][X-dir2].dagger()*tF0[X-dir2]*U[dir2][X-dir2])*0.25;
    }

}

template <typename group,typename atype=hila::arithmetic_type<group>>
atype measure_s_log(const GaugeField<group>& U) {
    // measure the log action for dir1<dir2
    // (just to have same normalization as with plaquette action)
    Reduction<double> stot=0;
    stot.allreduce(false).delayed(true);
    Field<group> Csum;
    foralldir(dir1) foralldir(dir2) if(dir1<dir2) {
        get_log_clover_mat(U,dir1,dir2,Csum);
        onsites(ALL) {
            stot+=-0.5*real(mul_trace(Csum[X],Csum[X]));
        }
    }
    return (atype)stot.value()/group::size();
}

template <typename group,typename atype=hila::arithmetic_type<group>>
void get_force_log_add(const GaugeField<group>& U,VectorField<Algebra<group>>& K, atype eps=1.0) {
    // compute gauge force for clover action and add result to K
    Field<group> Csum;

    foralldir(dir1) foralldir(dir2) if(dir1<dir2) {
        // get the 4 clover leave matrices C[][X] (ordered counter-clockwise in (dir1,dir2)-plane)
        // and their anti-hermitian traceless sum Csum[X]
        get_log_clover_mat(U,dir1,dir2,Csum);

        // define for each clover leave matrix the correspondin path of gauge links
        std::vector<Direction> paths[4]={{dir1,dir2,-dir1,-dir2},{dir2,-dir1,-dir2,dir1},{-dir1,-dir2,dir1,dir2},{-dir2,dir1,dir2,-dir1}};

        // the force on the link (Y,dir) from the clover term with foot point X is:
        // F[dir][Y] = d_{Y,dir}(0.5*Tr(Csum*Csum)) 
        //           = 0.5*Tr((d_{Y,dir}Csum)*Csum) + 0.5*Tr(Csum*(d_{Y,dir}Csum))
        //           = Tr((d_{Y,dir}Csum)*Csum)
        //           = 0.25*\sum_k Tr((d_{Y,dir}log(C[k][X]))*Csum[X])
        // with: d_{Y,dir}log(C[k][X]) = (d_{Y,dir}C[k][X])*C[k][X].dagger()
        //       where C[k][X] is the k-th clover leave matrix, 
        // Furthermore: (d^a_{Y,dir[0]}C[0][X])*C[0][X].dagger() = T^a
        //              (d^a_{Y,dir[1]}C[0][X])*C[0][X].dagger() = U[dir[0]][X].dagger() * T^a * U[dir[0]][X]
        //              ...
        //              ...
        // We therefore need only Csum/4 to compute the log-action force:
        Csum[ALL]*=0.25;
        for(int k=0; k<4; ++k) {
            // compute the gauge force from the first 4 links of all the Wilson loops that are summed in tC[k]
            // (is valid to do so since first 4 links are the same for all Wilson loops and projection on 
            // Lie-algebra is linear)
            get_wloop_force_from_wl_add(U,paths[k],Csum,eps,K);
        }
    }
}

template <typename group,typename atype=hila::arithmetic_type<group>>
void get_force_log(const GaugeField<group>& U,VectorField<Algebra<group>>& K, atype eps=1.0) {
    // determine gauge force for clover action and store result in K
    foralldir(d1) {
        K[d1][ALL]=0;
    }
    get_force_log_add(U,K,eps);
}


#endif