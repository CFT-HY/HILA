/** @file sun_overrelax.h */

#ifndef SUN_OVERRELAX_H
#define SUN_OVERRELAX_H

#include "sun_matrix.h"
#include "su2.h"

#ifdef 

/**
 * @brief \f$ SU(N) \f$ Overrelaxation using SU(2) subgroups
 *
 * @tparam T float/double/long double
 * @tparam N Number of colors
 * @param U \f$ SU(N) \f$ link to perform overrelaxation on
 * @param staple Staple to compute overrelaxation with
 */
template <typename T, int N>
void suN_overrelax(SU<N, T> &U, const SU<N, T> &staple) {
    /* Do overrelaxation by SU(2) subgroups */
    SU2<T> a;
    SU<2, T> u;
    SU<N, T> action = U * staple.dagger();

    for (int ina = 0; ina < N - 1; ina++)
        for (int inb = ina + 1; inb < N; inb++) {

            /* decompose the action into SU(2) subgroups using Pauli matrix
             * expansion
             * The SU(2) hit matrix is represented as
             * a0 + i * Sum j (sigma j * aj)
             */
            a = project_from_matrix(action, ina, inb);

            a.normalize();

            /* Now we need to multiply twice with the su2 element, thus;
             * do now a^2; for example
             *    a^2_1 = a0 a1 + a1 a0 - i a2 a3 s2 s3 - i a2 a3 s3 s2
             *          = 2 a0 a1
             * do also complex conjugate here
             */

            auto r = a.d * a.d - a.a * a.a - a.b * a.b - a.c * a.c;
            a.a *= -2 * a.d;
            a.b *= -2 * a.d;
            a.c *= -2 * a.d;
            a.d = r;

            /* Elements of SU(2) matrix */

            u = a.convert_to_2x2_matrix();

            /* Do SU(2) hit on U and action (to overrelax)  */

            U.mult_by_2x2_left(ina, inb, u);
            action.mult_by_2x2_left(ina, inb, u);

        } /* hits */
          /* check_unitarity( U );
           */
}


/**
 * @brief \f$ SU(N) \f$ full overrelaxation using SVD
 *
 * Following de Forcrand and Jahn, hep-lat/0503041
 *
 * Algorithm:
 * - let U be link matrix, S sum of staples so that action = -beta/N Re Tr U S^*
 * - Calculate svd  S = u D v.dagger(), where D is diagonal matrix of singular values, u,v \in U(N)
 * - Now M = u v.dagger() maximizes Tr(M S.dagger()) (in U(N)), but not SU(N) in general
 * - compute det S = rho exp(i\phi)
 * - Now e^(-i phi/N) u v.dagger() is SU(N)
 * - Optimize: find diagonal P = {e^{i\theta_j}} in SU(N), i.e. sum_i theta_i = 0 (mod 2pi), which
 *             maximizes ReTr( S.dagger() e^{-i phi/N} u P v.dagger()).
 * - Now matrix Z = e^(-i phi/N) u P v.dagger() in SU(N), and overrelax update
 *   U_new = Z U_old.dagger() Z
 * - accept/reject with change in action
 *
 *
 * Maximize ReTr( S.dagger() e^{-i phi/N} u P v.dagger()) = ReTr( e^{-i phi/N} D P ) =
 *          sum_j Re e^{i phi/N} e^{i\theta_j} D_j      j = 0 .. (N-1)
 * assuming D_min is smallest, def \theta_min = -\bar\theta = -sum_j \theta_j   j != i_min
 *
 * Taking derivatives d/d\theta_j and setting these == 0, we obtain N-1 eqns which can be linearized
 *     Re e^{-i \phi/N} (i e^{i\theta_j}) D_j - Re e^{-i \phi/N} (i e^{-i \bar\theta}) D_min = 0
 *                      (i - \theta_j)                           (i + \bar\theta)
 *
 * Let c == cos(phi/N), s == sin(phi/N)
 *    =>    (s - c \theta_j) D_j + (-s - c\bar\theta) D_min = 0
 *    =>    sum_k (D_j delta_jk + D_min ) \theta_k = (s/c)(D_j - D_min)
 *
 * This is a matrix equation of form M = (A + b F) theta = (A + w w^T) theta = rhs
 * where A is diagonal and F is filled with 1, and w = sqrt(b)[1 1 1 ..]^T.
 *
 * This can be solved with the Sherman-Morrison formula
 */

#define deForcrandJahn

template <typename T, int N, typename Btype>
int suN_full_overrelax(SU<N, T> &U, const SU<N, T> &S, Btype beta) {

    auto svdS = S.svd();
    auto phiN = S.det().arg() / N;

#ifdef deForcrandJahn
    // find smallest singular value
    int imin = 0;
    double smin = svdS.singularvalues.e(0);
    for (int i = 1; i < N; i++) {
        if (svdS.singularvalues.e(i) < smin) {
            smin = svdS.singularvalues.e(i);
            imin = i;
        }
    }

    DiagonalMatrix<N - 1, double> diag;
    int j = 0;
    for (int i = 0; i < N; i++) {
        if (i != imin)
            diag.e(j++) = svdS.singularvalues.e(i);
    }

    // Do the Sherman-Morrison inverse - diag will contain the phase angles
    diag.asVector() = hila::linalg::invert_diagonal_plus_constant_matrix(diag, smin) *
           (tan(phiN) * (diag - smin).asVector());

    // and construct diagonal P - absorb e^{-i phi/N} to P
    DiagonalMatrix<N, Complex<T>> P;
    j = 0;
    for (int i = 0; i < N; i++) {
        if (i != imin)
            P.e(i) = expi(diag.e(j++) - phiN);
    }
    P.e(imin) = expi(-trace(diag) - phiN);
#else
    // This is Narayanan-Neuberger method.  It has slightly worse acceptance
    // to Forcrand + Jahn
    DiagonalMatrix<N, Complex<T>> P;
    P = expi(-phiN);
#endif

    // make unitary matrix Z
    auto Z = svdS.U * P * svdS.V.dagger();

    // candidate for new U
    Z = Z * U.dagger() * Z;

    // exp(old-new) > random()
    if (exp(beta / N * real(mul_trace(Z - U, S.dagger()))) > hila::random()) {
        // accepted
        U = Z;
        return 1;
    } else
        return 0;
}


#endif