#ifndef KENNEDY_PENDLETON_H
#define KENNEDY_PENDLETON_H

#include "plumbing/hila.h"
#include "datatypes/su2.h"

// this is might be too much to have in a single kernel

/**
 * @brief \f$ SU(2) \f$ heatbath
 * @details Generate an \f$ SU(2)\f$ matrix according to Kennedy-Pendleton heatbath. 
 * This assumes that the gauge link U to be generated contributes to the action as Re Tr U.S,
 * where S is the staple. Note sign and normalization of the staple (beta is absorbed in S).
 * @tparam T float-type
 * @param staple Staple to compute heatbath with. NB: In general the staple is NOT an SU(2) matrix, 
 * however we can describe it using the SU2 class because (i) The sum of SU(2) matrices is proportional to SU(2), 
 * (ii) SU2 class does not require that the determinant is always normalized to 1.
 * 
 * @return SU2<T> 
 */
template <typename T>
SU2<T> KennedyPendleton(const SU2<T> &staple) {

    // action Tr U.S => weight exp(-Tr U.S)  (for SU2 the trace is real).
    // Following K-P we need to bring this to form exp(xi Tr h.u), 
    // where h = U (our link) and u is in eq. (10) of K-P paper.
    // But since the staple is proportional to SU(2) matrix, 
    // this just means we normalize the staple to SU(2) and flip its sign.
    // The K-P algorithm generates a new SU2 matrix a = h.u 
    // from which the new link h is u^+.a

    const T xi = sqrt(det(staple));
    const SU2<T> u = -staple / xi; // this is now in SU(2)

    // 'alpha' that appears in eq (15)
    T alpha = 2.0*xi;

    // Now generate new a = h.u, section 4 of K-P

    int nloops;
    const int looplim = 100;
    const double pi2 = M_PI * 2.0;
    double r1, r2, r3, delta;
    nloops = 0;
    do {
        r1 = -log(1.0 - hila::random()) / alpha;
        r2 = -log(1.0 - hila::random()) / alpha;
        r3 = cos(pi2 * hila::random());
        r3 *= r3;
        delta = r1*r3 + r2;

        r3 = hila::random(); // R'''
        nloops++;
    } while (nloops < looplim && r3*r3 > 1.0 - 0.5*delta);

    T a0 = 1.0 - delta;
    if (nloops >= looplim) {
        // Failed to find new a0...
        // TODO some error or warning
        a0 = 1e-9;
    }

    // Generate random vector on S2, radius sqrt(1.0-a0**2)
    T rd = 1.0 - a0*a0;

    if (rd <= 0.0) {
        // negative radius^2, something went wrong...
        // TODO some error or warning
        rd = 0.0;
    }

    // ... and put them in the 'a' matrix. Notation gets messy, reminder:
    // SU2 = d + a i\sigma_1 + b i\sigma_2 + c i\sigma_3$
    SU2<T> a;
    a.d = a0;

    r3 = 2.0 * hila::random() - 1.0;
    a.c = sqrt(rd) * r3;
    T theta = pi2 * hila::random();
    rd = sqrt(rd * (1.0 - sqr(r3)));
    a.b = rd * cos(theta);
    a.a = rd * sin(theta);
    // Could this be optimized by removing sin, cos? http://corysimon.github.io/articles/uniformdistn-on-sphere/

    // Now just get the new link variable
    return u.dagger() * a;
}


#endif