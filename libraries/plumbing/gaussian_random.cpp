#ifndef RANDOM_H_
#define RANDOM_H_

#include <cmath>
#include "plumbing/random.h"
#include "plumbing/defs.h"

// gaussian rng generation routines ----------------------


static constexpr double pi = 3.14159265358979;

double gaussian_ran(double variance=0.5) {
  static double second;
  static bool draw_new = true;
  if( draw_new ) {
    double first, phi, urnd, r;
    phi = 2.0 * pi * hila_random();
    urnd = (1.0-hila_random());
    r  = sqrt( -log(urnd) * (2.0 * variance) );
    first =  r*sin(phi);
    second =  r*cos(phi);

    draw_new = false;
    return first;
  } else {
    draw_new = true;
    return second;
  }
}

#define VARIANCE 0.5

double gaussian_ran2(double & out2) {

  double phi, urnd, r;
  phi = 2.0 * pi * hila_random();
  urnd =  1.0 - hila_random();
  r  = sqrt( -log(urnd) * (2.0 * VARIANCE) );
  out2 = r*cos(phi);
  return r*sin(phi);
}

float gaussian_ran2(float & out2) {
  double phi, urnd, r;
  phi = 2.0 * pi * hila_random();
  urnd =  1.0 - hila_random();
  r  = sqrt( -log(urnd) * (2.0 * VARIANCE) );
  out2 = (float) (r*cos(phi));
  return (float) (r*sin(phi));
}

#endif