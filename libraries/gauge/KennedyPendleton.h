#ifndef KENNEDY_PENDLETON_H
#define KENNEDY_PENDLETON_H

#include "plumbing/hila.h"
#include "datatypes/su2.h"

/*************************************************************
 *                                                           *
 *     Kennedy-Pendleton su2 heatbath code                   *
 *     p contains the 'staple'                               *
 *     act contains the resulting action                     *
 *                                                           *
 *     act = u*p', weight exp[-u*p']                         *
 *                                                           *
 ************************************************************/

#define looplim 100

// this is might be too much to have in a single kernel

template <typename T>
SU2<T> KennedyPendleton(SU2<T> &staple) {
    int nloops;
    SU2<T> a, ua;
    double e, f, r1, r2, r3, rd, theta;
    const double pi2 = M_PI * 2.0;

    /* first, normalize and move to u [-] because want exp[up] */

    ua = staple / (-sqrt(det(staple)));

    /*     generate random su2 variable according to Pendleton
     *     and Kennedy, PL 156B(1985)393
     */

    nloops = 0;
    do {
        r1 = log(1.0 - hila::random());
        r2 = log(1.0 - hila::random());
        r3 = cos(pi2 * hila::random());
        r3 *= r3;
        e = 1.0 + (r1 + r2 * r3) * b;
        r3 = hila::random();
        f = e - 2.0 * sqr(r3) + 1.0;
    } while (nloops++ < looplim && (f <= 0.0));

    if (nloops >= looplim) {
        // printf(" ---> K-P loop limit %d exceeded!\n",looplim);
        // printf(" ---> staple magnitude %g\n",1/b);
        /* exit(1); */
        e = 1e-9;
    }

    /* generate random vector on s2, radius sqrt(1.0-a.d**2) */
    a.d = e;
    rd = 1.0 - sqr(e);
    if (rd <= 0.0) {
        // printf(" ---> neg. K-P sqrt, val %g\n",rd);
        rd = 0.0;
        // fflush(stdout);
    }
    r3 = 2.0 * hila::random() - 1.0;
    a.c = sqrt(rd) * r3;
    theta = pi2 * hila::random();
    rd = sqrt(rd * (1.0 - sqr(r3)));
    a.b = rd * cos(theta);
    a.a = rd * sin(theta);


    /* a = -u*p' => u = -a*p = a*ua */
    return a * ua;
}


#endif