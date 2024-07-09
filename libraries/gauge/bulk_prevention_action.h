/** @file bulk_prevention_action.h */

#ifndef BULK_PREVENTION_ACTION_H_
#define BULK_PREVENTION_ACTION_H_

#include "hila.h"

// functions for bulk-preventing action, cf. arXiv:2306.14319 (with n=2)
template <typename T>
T ch_inv(const T &U) {
    // compute inverse of the square matrix U, using the
    // Cayley-Hamilton theorem (Faddeev-LeVerrier algorithm)
    T tB[2];
    int ip = 0;
    tB[ip] = 1.;
    auto tc = trace(U);
    tB[1 - ip] = U;
    for (int k = 2; k <= T::size(); ++k) {
        tB[1 - ip] -= tc;
        mult(U, tB[1 - ip], tB[ip]);
        tc = trace(tB[ip]) / k;
        ip = 1 - ip;
    }
    return tB[ip] / tc;
}

template <typename T>
T bp_UAmat(const T &U) {
    // compute U*A(U) with the A-matrix from Eq. (B3) of arXiv:2306.14319 for n=2
    T tA1 = U;
    tA1 += 1.;
    tA1 *= 0.5;
    T tA2 = ch_inv(tA1);
    mult(tA2, tA2.dagger(), tA1);
    return U * tA1 * tA1 * tA2;
}

template <typename T>
T bp_iOsqmat(const T &U) {
    // compute matrix inside the trace on r.h.s. of Eq. (B1) of arXiv:2306.14319 for n=2
    // (without identiy matrix subtraction)
    T tA1 = U;
    tA1 += 1.;
    tA1 *= 0.5;
    T tA2 = ch_inv(tA1);
    mult(tA2, tA2.dagger(), tA1);
    return tA1 * tA1;
}

template <typename group>
double measure_s_bp(const GaugeField<group> &U) {
    // measure the BP plaquette action
    Reduction<double> plaq = 0;
    plaq.allreduce(false).delayed(true);
    foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {
        U[dir2].start_gather(dir1, ALL);
        U[dir1].start_gather(dir2, ALL);
        onsites(ALL) {
            plaq += real(trace(bp_iOsqmat(U[dir1][X] * U[dir2][X + dir1] *
                                          (U[dir2][X] * U[dir1][X + dir2]).dagger()))) /
                        group::size() -
                    1.0;
        }
    }
    return plaq.value();
}

template <typename group, typename atype = hila::arithmetic_type<group>>
void get_force_bp_add(const GaugeField<group> &U, VectorField<Algebra<group>> &K, atype eps = 1.0) {
    // compute force for BP action for n=2 according to Eq. (B5) of arXiv:2306.14319 and add it to K
    Field<group> fmatp;
    Field<group> fmatmd1;
    Field<group> fmatmd2;
    atype teps = 2.0 * eps;
    foralldir(d1) {
        foralldir(d2) if (d1 < d2) {
            U[d2].start_gather(d1, ALL);
            U[d1].start_gather(d2, ALL);
            onsites(ALL) {
                fmatp[X] = bp_UAmat(U[d1][X] * U[d2][X + d1] * (U[d2][X] * U[d1][X + d2]).dagger());
                fmatmd1[X] = (fmatp[X] * U[d2][X]).dagger() *
                             U[d2][X]; // parallel transport fmatp[X].dagger() to X+d2
                fmatmd2[X] =
                    U[d1][X].dagger() * fmatp[X] * U[d1][X]; // parallel transport fmatp[X] to X+d1
            }
            fmatmd1.start_gather(-d2, ALL);
            fmatmd2.start_gather(-d1, ALL);
            onsites(ALL) {
                K[d1][X] -= (fmatmd1[X - d2] + fmatp[X]).project_to_algebra_scaled(teps);
                K[d2][X] -= (fmatmd2[X - d1] - fmatp[X]).project_to_algebra_scaled(teps);
            }
        }
    }
}

template <typename group, typename atype = hila::arithmetic_type<group>>
void get_force_bp(const GaugeField<group> &U, out_only VectorField<Algebra<group>> &K,
                  atype eps = 1.0) {
    // compute force for BP action for n=2 according to Eq. (B5) of arXiv:2306.14319 and write it to
    // K
    foralldir(d1) {
        K[d1][ALL] = 0;
    }
    get_force_bp_add(U, K, eps);
}


#endif