#include "sun.h"
#include "../plumbing/defs.h"
#include <cmath>

#define VARIANCE 0.5
#define MYSINF(X) sin(X)
#define MYCOSF(X) cos(X)

template<typename radix>
inline radix gaussian_ran2 (radix* out2) 
{
  double phi, urnd, r;
  phi = 2.0 * 3.141592654 * (double) hila_random();
  urnd = (double)(1.0-hila_random());
  r  = sqrt( -log(urnd) * (2.0 * VARIANCE) );
  *out2 = (r*MYCOSF(phi));
  return (radix)(r*MYSINF(phi));
}

template<typename radix>
radix SU<2,radix>::sqr(){
    return a*a + b*b + c*c + d*d;
}

template<typename radix>
SU<2,radix> & SU<2,radix>::normalize(){
    radix sq = this.sqr();
    a /= sq;
    b /= sq;
    c /= sq;
    d /= sq;   
    return *this;
}

template<typename radix>
SU<2,radix> & SU<2,radix>::reunitarize(){
    return this.normalize();
}

template<typename radix>
SU<2,radix> & SU<2,radix>::random(){
    radix one, two;
    one = gaussian_ran2(&two);
    a = one;
    b = two;
    one = gaussian_ran2(&two);
    c = one;
    d = two;
    return this.normalize();
}

template<typename radix>
SU<2,radix> & SU<2,radix>::inv(){
    return *this;
}




