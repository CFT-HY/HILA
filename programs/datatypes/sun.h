#include "cmplx.h"
#include "general_matrix.h"
#include "../plumbing/defs.h"
#include "../plumbing/mersenne.h" //has to be included
#include <cmath>

//////////////////
///
/// Template class for SU matrices
/// 
/// Uses specialized format for cases n = 2 and n = 3 
/// (parametrized by pauli matrices) for efficiency
/// otherwise behaves like a regular complex matrix
///
/// Implementation is found in header file so that linker 
/// finds specializations
//////////////////

//TODO: SU(2) vector class + conversion to and from matrix rep.

template<int n, typename radix>
class SU : public matrix<n,n,cmplx<radix>> {
    public: 
        void reunitarize();
};

template<typename radix>
class SU<2, radix> { 
    public: 

        SU(){};
        ~SU(){};
        SU<2,radix> & normalize(); //normalize elements
        SU<2,radix> & reunitarize(); 
        SU<2,radix> & random(); //generate random SU2 element 
        SU<2,radix> & inv(); //invert matrix and return reference to self
        radix sqr(); //calculate square of all elements
        radix det(); //determinant
        radix tr(); //trace 

        SU<2,radix> operator + (const SU<2,radix> &); //return copies of result 
        SU<2,radix> operator - (const SU<2,radix> &);
        SU<2,radix> operator * (const SU<2,radix> &);
        SU<2,radix> operator + (const radix &);
        SU<2,radix> operator - (const radix &);
        SU<2,radix> operator * (const radix &);
        SU<2,radix> & operator += (const SU<2,radix> &); //same ops as above, except store result in this matrix
        SU<2,radix> & operator -= (const SU<2,radix> &);
        SU<2,radix> & operator *= (const SU<2,radix> &);
        SU<2,radix> & operator += (const radix &);
        SU<2,radix> & operator -= (const radix &);
        SU<2,radix> & operator *= (const radix &);
        SU<2,radix> & operator = (const SU<2,radix> &);

    private:

        radix a, b, c, d;      
};

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
    radix sq = sqrt(this->sqr());
    a /= sq;
    b /= sq;
    c /= sq;
    d /= sq;   
    return *this;
}

template<typename radix>
SU<2,radix> & SU<2,radix>::reunitarize(){
    return this->normalize();
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
    return this->normalize();
}

template<typename radix>
SU<2,radix> & SU<2,radix>::inv(){
    return *this;
}
