#ifndef SUN_OVERRELAX_H
#define SUN_OVERRELAX_H

#include "sun_matrix.h"
#include "su2.h"


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

#endif