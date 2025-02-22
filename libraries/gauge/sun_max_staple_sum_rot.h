/** @file sun_max_staple_sum_rot.h */

#ifndef SUN_MAX_STAPLE_SUM_ROT_H
#define SUN_MAX_STAPLE_SUM_ROT_H

#include <numeric>
#include "sun_matrix.h"
#include "plumbing/shuffle.h"
#include "tools/floating_point_epsilon.h"


/**
 * @brief Find \f$ SU(N) \f$ matrix U that maximizes ReTr(U*staple) for given staple
 *
 * @tparam T float/double/long double
 * @tparam N Number of colors
 * @param U \f$ SU(N) \f$ matrix that maximizes ReTr(U*staple)
 * @param staple Staple 
 */
template <typename T, int N>
SU<N, T> suN_max_staple_sum_rot(const SU<N, T> &staple) {

    SU<N, T> tl0(staple);
    SU<N, T> U(1.0);

    // get a random permutation of the list of possible matrix indices
    std::array<int, N> tperm;
    std::iota(tperm.begin(), tperm.end(), 0);
    hila::shuffle(tperm);

    // run thourgh the different SU(2) subgroups of SU(N) (in the order given by the permuatation
    // iperm), to build up the SU(N) matrix that that has maximum overlap with l0:
    int niter = 15 * (N - 1) * N;
    if (N == 2)
        niter = 1;

    int iter, tt1, tt2, t1, t2,ic1;
    Complex<T> tvnn11, tvnn12, tvnn21, tvnn22, ttl01, ttl02;
    T tu1, tu2, tu3, tu4, k, ki, lnorm, dlnorm, dilnorm;
    for (iter = 0; iter < niter; ++iter) {
        for (tt1 = 0; tt1 < N - 1; ++tt1) {
            t1 = tperm[tt1];
            for (tt2 = tt1 + 1; tt2 < N; ++tt2) {
                t2 = tperm[tt2];

                tvnn11 = tl0.e(t1, t1);
                tvnn12 = tl0.e(t1, t2);
                tvnn21 = tl0.e(t2, t1);
                tvnn22 = tl0.e(t2, t2);

                // get inverse of SU(2) projection of the 2x2 matrix tvnnXX:

                // determine k=sqrt(det(tvnnXX)):
                tu1 = 0.5 * (tvnn11.real() + tvnn22.real());
                tu2 = 0.5 * (tvnn12.imag() + tvnn21.imag());
                tu3 = 0.5 * (tvnn12.real() - tvnn21.real());
                tu4 = 0.5 * (tvnn11.imag() - tvnn22.imag());
                k = sqrt(tu1 * tu1 + tu2 * tu2 + tu3 * tu3 + tu4 * tu4);

                // determine tui=inverse(tvnnXX/k):
                //  normalize tu:
                ki = 1. / k;
                tu1 *= ki;
                tu2 *= ki;
                tu3 *= ki;
                tu4 *= ki;
                // invert tu:
                // (re-use variables tvnnXX to hold tui)
                tvnn11.re = tu1;
                tvnn11.im = -tu4;
                tvnn12.re = -tu3;
                tvnn12.im = -tu2;
                tvnn21.re = tu3;
                tvnn21.im = -tu2;
                tvnn22.re = tu1;
                tvnn22.im = tu4;
                for (ic1 = 0; ic1 < N; ++ic1) {
                    // apply SU(2) subgroup rotation "tui" to tl0:
                    ttl01 = tl0.e(t1, ic1);
                    ttl02 = tl0.e(t2, ic1);
                    tl0.e(t1, ic1) = tvnn11 * ttl01 + tvnn12 * ttl02;
                    tl0.e(t2, ic1) = tvnn21 * ttl01 + tvnn22 * ttl02;
                }

                for (ic1 = 0; ic1 < N; ++ic1) {
                    // apply SU(2) subgroup rotation "tui" to U:
                    ttl01 = U.e(t1, ic1);
                    ttl02 = U.e(t2, ic1);
                    U.e(t1, ic1) = tvnn11 * ttl01 + tvnn12 * ttl02;
                    U.e(t2, ic1) = tvnn21 * ttl01 + tvnn22 * ttl02;
                }
            }
        }
        // check whether tl0 has reached the form H + c*id:
        ttl01 = tl0.e(N - 1, N - 1);
        lnorm = ttl01.squarenorm();
        dlnorm = 0;
        dilnorm = ttl01.imag();
        for (t1 = 0; t1 < N - 1; ++t1) {
            ttl01 = tl0.e(t1, t1);
            lnorm += ttl01.squarenorm();
            dlnorm += ::squarenorm(ttl01.imag() - dilnorm);
            for (t2 = t1 + 1; t2 < N; ++t2) {
                dlnorm += 0.5 * ::squarenorm(tl0.e(t1, t2) - tl0.e(t2, t1).conj());
            }
        }
        lnorm = ::sqrt(lnorm);
        dlnorm = ::sqrt(dlnorm);
        if (dlnorm < fp<T>::epsilon * (lnorm + (T)N) * (T)N) {
            break;
        }
    }

    //U.reunitarize();

    return U;
}

/**
 * @brief \f$ SU(N) \f$ overrelaxation using \f$ SU(N) \f$ matrix obtained 
 * from max_staple_sum_rot() routine
 *
 * @tparam T float/double/long double
 * @tparam N Number of colors
 * @param U \f$ SU(N) \f$ link to perform overrelaxation on
 * @param staple Staple to compute overrelaxation with
 */
template <typename T, int N, typename Btype>
int suN_overrelax(SU<N, T> &U, const SU<N, T> &staple, Btype beta) {
    SU<N, T> tstaple = staple.dagger();
    SU<N, T> tU = suN_max_staple_sum_rot(tstaple);
    if(N==2) {
        U = tU * U.dagger() * tU;
        return 1;
    } else {
        tU = tU * U.dagger() * tU;
        T ds = real(mul_trace((tU - U), tstaple));
        if(ds > 0 || hila::random() < exp(beta/N * ds)) {
            U = tU;
            return 1;
        }
    }
    return 0;
}

#endif