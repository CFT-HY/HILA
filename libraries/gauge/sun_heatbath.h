#ifndef HILA_SUN_HEATBATH_H
#define HILA_SUN_HEATBATH_H

/** @file sun_heatbath.h */

/* Kennedy-Pendleton quasi heat bath on SU(2) subgroups */

/* This does the heatbath with
 *       beta ( 1 - 1/Ncol Tr U S^+ )  S is the sum of staples.
 *  ->  -beta/Ncol Tr u (US+) = -beta/Ncol u_ab (US+)_ba
 *  here u is a SU(2) embedded in SU(N); u = su2 * 1.
 */

#include "sun_matrix.h"
#include "su2.h"

/**
 * @brief \f$ SU(N) \f$ heatbath
 * @details Kennedy-Pendleton quasi heat bath on \f$ SU(2)\f$ subgroups
 * @tparam T Group element type such as Real or Complex
 * @tparam N Number of colors
 * @param U \f$ SU(N) \f$ link to perform heatbath on
 * @param staple Staple to compute heatbath with
 * @param beta
 * @return double
 */
template <typename T, int N>
double suN_heatbath(SU<N, T> &U, const SU<N, T> &staple, double beta) {
    // K-P quasi-heat bath by SU(2) subgroups

    double utry, uhit;
    double xr1, xr2, xr3, xr4;

    double r, r2, rho, z;
    double al, d, xl, xd;
    int k, test;
    double b3;
    SU<N, T> action;
    SU<2, T> h2x2;
    SU2<T> v, a, h;
    double pi2 = M_PI * 2.0;

    b3 = beta / N;

    utry = uhit = 0.0;

    /* now for the qhb updating */

    action = U * staple.dagger();

    /*  SU(2) subgroups */
    for (int ina = 0; ina < N - 1; ina++)
        for (int inb = ina + 1; inb < N; inb++) {

            /* decompose the action into SU(2) subgroups using
             * Pauli matrix expansion
             * The SU(2) hit matrix is represented as
             * a0 + i * Sum j (sigma j * aj)
             * Let us actually store dagger, v = action(a,b).dagger
             */
            v.d = action.e(ina, ina).re + action.e(inb, inb).re;
            v.c = -(action.e(ina, ina).im - action.e(inb, inb).im);
            v.a = -(action.e(ina, inb).im + action.e(inb, ina).im);
            v.b = -(action.e(ina, inb).re - action.e(inb, ina).re);

            z = sqrt(v.det());

            v /= z;

            /* end norm check--trial SU(2) matrix is a0 + i a(j)sigma(j)*/
            /* test            */
            /* now begin qhb */
            /* get four random numbers */

            xr1 = log(1.0 - hila::random());
            xr2 = log(1.0 - hila::random());
            xr3 = hila::random();
            xr4 = hila::random();

            xr3 = cos(pi2 * xr3);


            //   generate a0 component of suN matrix
            //
            //   first consider generating an su(2) matrix h
            //   according to exp(bg/Ncol * re tr(h*s))
            //   rewrite re tr(h*s) as re tr(h*v)z where v is
            //   an su(2) matrix and z is a real normalization constant
            //   let v = z*v. (z is 2*xi in k-p notation)
            //   v is represented in the form v(0) + i*sig*v (sig are pauli)
            //   v(0) and vector v are real
            //
            //   let a = h*v and now generate a
            //   rewrite beta/3 * re tr(h*v) * z as al*a0
            //   a0 has prob(a0) = n0 * sqrt(1 - a0**2) * exp(al * a0)


            al = b3 * z;

            // let a0 = 1 - del**2
            // get d = del**2
            // such that prob2(del) = n1 * del**2 * exp(-al*del**2)


            d = -(xr2 + xr1 * xr3 * xr3) / al;

            /* monte carlo prob1(del) = n2 * sqrt(1 - 0.5*del**2)
             * then prob(a0) = n3 * prob1(a0)*prob2(a0)
             */

            if ((1.00 - 0.5 * d) <= xr4 * xr4) {
                if (al > 2.0) { /* k-p algorithm */
                    test = 0;
                    for (k = 0; k < 40 && !test; k++) {
                        /*  get four random numbers */
                        xr1 = log(1.0 - hila::random());
                        xr2 = log(1.0 - hila::random());
                        xr3 = hila::random();
                        xr4 = hila::random();

                        xr3 = cos(pi2 * xr3);

                        d = -(xr2 + xr1 * xr3 * xr3) / al;
                        if ((1.00 - 0.5 * d) > xr4 * xr4)
                            test = 1;
                    }
                    utry += k;
                    uhit++;
#ifndef __CUDA_ARCH__
                    // TODO: Make parallel-warning out of this
                    // hila::out0 << "site took 20 kp hits\n";
#endif
                } else { /* now al <= 2.0 */

                    /* creutz algorithm */
                    xl = exp((double)(-2.0 * al));
                    xd = 1.0 - xl;
                    test = 0;
                    double a0;
                    for (k = 0; k < 40 && test == 0; k++) {
                        /*        get two random numbers */
                        xr1 = hila::random();
                        xr2 = hila::random();

                        r = xl + xd * xr1;
                        a0 = 1.00 + log((double)r) / al;
                        if ((1.0 - a0 * a0) > xr2 * xr2)
                            test = 1;
                    }
                    d = 1.0 - a0;
                    utry += k;
                    uhit++;
#ifndef __CUDA_ARCH__
                    // hila::out0 << "site  took 20 creutz hits\n";
#endif
                } /* endif al */

            } else {
                /* direct hit */
                utry += 1.0;
                uhit += 1.0;
            }

            /*  generate full su(2) matrix and update link matrix*/

            /* find a0  = 1 - d*/
            a.d = 1.0 - d;
            /* compute r */
            r2 = 1.0 - a.d * a.d;
            r2 = fabs(r2);
            r = sqrt(r2);

            /* compute a3 */
            a.c = (2.0 * hila::random() - 1.0) * r;

            /* compute a1 and a2 */
            rho = r2 - a.c * a.c;
            rho = sqrt(fabs(rho));

            /*xr2 is a random number between 0 and 2*pi */
            xr2 = pi2 * hila::random();
            a.a = rho * cos((double)xr2);
            a.b = rho * sin((double)xr2);

            /* now do the updating.  h = a*v^dagger, new u = h*u */
            h = a * v; // v was daggerized above

            /* Elements of SU(2) matrix */
            h2x2 = h.convert_to_2x2_matrix();

            /* update the link and 'action' */

            U.mult_by_2x2_left(ina, inb, h2x2);
            action.mult_by_2x2_left(ina, inb, h2x2);

        } /*  hits */

    /* check_unitarity( U );
     */
    return (uhit / utry);

} /* site */

#endif
