#ifndef RANDOM
#define RANDOM

#include "plumbing/defs.h"
#include <cmath>

// gaussian rng generation routines ----------------------


#pragma hila loop_function
static double gaussian_ran(double variance=0.5)
{
  static double second;
  static bool draw_new = true;
  if( draw_new ) {
    double first, phi, urnd, r;
    phi = 2.0 * 3.141592654 * (double) hila_random();
    urnd = (double)(1.0-hila_random());
    r  = sqrt( -log(urnd) * (2.0 * variance) );
    first = (double)(r*sin(phi));
    second = (r*cos(phi));

    draw_new = false;
    return first;
  } else {
    draw_new = true;
    return second;
  }
}

#define VARIANCE 0.5
#define MYSINF(X) sin(X)
#define MYCOSF(X) cos(X)

template<typename radix>
#pragma hila loop_function
radix gaussian_ran2(radix* out2) 
{
  double phi, urnd, r;
  phi = 2.0 * 3.141592654 * (double) hila_random();
  urnd = (double)(1.0-hila_random());
  r  = sqrt( -log(urnd) * (2.0 * VARIANCE) );
  *out2 = (r*MYCOSF(phi));
  return (radix)(r*MYSINF(phi));
}


#endif