#include <cmath>
#include "plumbing/random.h"
#include "plumbing/defs.h"

///////////////////////////////////////////////////////////////////
/// gaussian rng generation routines
/// By default these give random numbers with variance 1.0, i.e.
/// probability distribution is
///           exp( -x*x/2 ), so < x^2 > = 1
/// If you want random numbers with variance sigma, multiply the
/// result by  sqrt(sigma):
///       sqrt(sigma) * gaussrand();
////////////////////////////////////////////////////////////////////

#define VARIANCE 1.0

static constexpr double pi = 3.14159265358979;

double hila::gaussrand2(double &out2) {

    double phi, urnd, r;
    phi = 2.0 * pi * hila::random();

    // this should not really trigger
    do {
        urnd = 1.0 - hila::random();
    } while (urnd == 0.0);

    r = sqrt(-::log(urnd) * (2.0 * VARIANCE));
    out2 = r * cos(phi);
    return r * sin(phi);
}

#if !defined(CUDA) && !defined(HIP)

double hila::gaussrand() {
    static double second;
    static bool draw_new = true;
    if (draw_new) {
        draw_new = false;
        return hila::gaussrand2(second);
    }
    draw_new = true;
    return second;
}

#else

// Cuda and other stuff which does not accept
// static variables - just throw away another gaussian number.

//#pragma hila loop function contains rng
double hila::gaussrand() {
    double second;
    return hila::gaussrand2(second);
}

#endif
