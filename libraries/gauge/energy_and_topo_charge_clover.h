/** @file energy_and_topo_charge_clover.h */

#ifndef ENERGY_AND_TOPO_CHARGE_CLOVER_H_
#define ENERGY_AND_TOPO_CHARGE_CLOVER_H_

#include "hila.h"


template <typename group, typename atype = hila::arithmetic_type<group>>
void measure_topo_charge_and_energy_clover(const GaugeField<group> &U, atype &qtopo_out,
                                           atype &energy_out) {
    // measure topological charge and field strength energy of the gauge field, using the
    // clover matrices as components of the field strength tensor

    Reduction<double> qtopo = 0;
    Reduction<double> energy = 0;
    qtopo.allreduce(false).delayed(true);
    energy.allreduce(false).delayed(true);

#if NDIM == 4
    Field<group> F[6];
    // F[0]: F[0][1], F[1]: F[0][2], F[2]: F[0][3],
    // F[3]: F[1][2], F[4]: F[1][3], F[5]: F[2][3]
    /*
    int k=0;
    foralldir(dir1) foralldir(dir2) if(dir1<dir2) {
        onsites(ALL) {
            // clover operator as in eq. (2.9) of Nuclear Physics B259 (1985) 572-596
            F[k][X]=0.25*(U[dir1][X]*U[dir2][X+dir1]*(U[dir2][X]*U[dir1][X+dir2]).dagger()
                +U[dir2][X]*(U[dir2][X-dir1]*U[dir1][X-dir1+dir2]).dagger()*U[dir1][X-dir1]
                +(U[dir2][X-dir1-dir2]*U[dir1][X-dir1]).dagger()*U[dir1][X-dir1-dir2]*U[dir2][X-dir2]
                +U[dir2][X-dir2].dagger()*U[dir1][X-dir2]*U[dir2][X+dir1-dir2]*U[dir1][X].dagger());
            // Lie-algebra projection:
            //F[k][X]=0.5*(F[k][X]-F[k][X].dagger());
            F[k][X]-=F[k][X].dagger();
            F[k][X]*=0.5;
            F[k][X]-=trace(F[k][X])/group::size();
        }
        ++k;
    }
    */

    Field<group> tF0, tF1;

    int k = 0;
    foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {
        U[dir2].start_gather(dir1, ALL);
        U[dir1].start_gather(dir2, ALL);

        onsites(ALL) {
            // dir1-dir2-plaquette that starts and ends at X; corresponds to F[dir1][dir2]
            // at center location X+dir1/2+dir2/2 of plaquette:
            tF0[X] = U[dir1][X] * U[dir2][X + dir1] * (U[dir2][X] * U[dir1][X + dir2]).dagger();
            // project to Lie-algebra (anti-hermitian trace-free)
            tF0[X] -= tF0[X].dagger();
            tF0[X] *= 0.5;
            tF0[X] -= trace(tF0[X]) / group::size();
            // parallel transport to X+dir1
            tF1[X] = U[dir1][X].dagger() * tF0[X] * U[dir1][X];
        }

        tF1.start_gather(-dir1, ALL);
        onsites(ALL) {
            tF0[X] += tF1[X - dir1];
        }

        U[dir2].start_gather(-dir2, ALL);
        tF0.start_gather(-dir2, ALL);
        onsites(ALL) {
            // get F[dir1][dir2] at X from average of the (parallel transported) F[dir1][dir2] from
            // the centers of all dir1-dir2-plaquettes that touch X :
            F[k][X] =
                (tF0[X] + U[dir2][X - dir2].dagger() * tF0[X - dir2] * U[dir2][X - dir2]) * 0.25;
        }
        ++k;
    }

    onsites(ALL) {
        qtopo += real(mul_trace(F[0][X], F[5][X]));
        qtopo += -real(mul_trace(F[1][X], F[4][X]));
        qtopo += real(mul_trace(F[2][X], F[3][X]));

        energy += F[0][X].squarenorm();
        energy += F[1][X].squarenorm();
        energy += F[2][X].squarenorm();
        energy += F[3][X].squarenorm();
        energy += F[4][X].squarenorm();
        energy += F[5][X].squarenorm();
    }
#endif
    qtopo_out = (atype)qtopo.value() / (4.0 * M_PI * M_PI);
    energy_out = (atype)energy.value();
}

#endif