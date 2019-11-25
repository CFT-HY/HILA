#include "cmplx.h"
#include "general_matrix.h"
#include "../plumbing/defs.h"
#include "../plumbing/mersenne.h" //has to be included
#include <cmath>

// matrix multiplication routines ------------------------

#define nn_a(x,y) (x.d*y.a + x.a*y.d - x.b*y.c + x.c*y.b)
#define nn_b(x,y) (x.d*y.b + x.b*y.d - x.c*y.a + x.a*y.c)
#define nn_c(x,y) (x.d*y.c + x.c*y.d - x.a*y.b + x.b*y.a)
#define nn_d(x,y) (x.d*y.d - x.a*y.a - x.b*y.b - x.c*y.c)

// gaussian rng generation routines ----------------------

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

//--------------------------------------------------------

//////////////////
///
/// SU(N) matrix class
/// 
///
/// Implementations are found in this header file, since non-specialized
/// templates have to be visible to the tranlation units that use them
///
//////////////////

//TODO: SU(2) vector class + conversion to and from matrix rep.
//TODO: write template for lhs scalar times matrix 

template<int n, typename radix>
class SU : public matrix<n,n,cmplx<radix>> {
    public: 
        void reunitarize();
};

template<typename radix>
class SU2<radix> { 
    public: 

        SU2() : a(0), b(0), c(0), d(1) {}
        SU2(radix * vals) : a(vals[0]), b(vals[1]), c(vals[2]), d(vals[3]) {}
        ~SU2(){}

        SU2<radix> & normalize(); //normalize elements
        SU2<radix> & reunitarize(); 
        SU2<radix> & random(); //generate random SU2 element 
        SU2<radix> & inv(); //invert matrix and return reference to self
        SU2<radix> & adj(); 
        
        inline radix sqr(); //calculate square of all elements
        inline radix tr(); //trace 
        radix det(); //determinant

        SU2<radix> operator + (const SU2<radix> &); //basic operations. These versions return separate copy of result
        SU2<radix> operator - (const SU2<radix> &);
        inline SU2<radix> operator * (const SU2<radix> &);
        SU2<radix> operator * (const radix &);
        SU2<radix> operator / (const radix &);
        SU2<radix> & operator += (const SU2<radix> &); //same ops as above, except store result in lhs matrix
        SU2<radix> & operator -= (const SU2<radix> &);
        SU2<radix> & operator *= (const SU2<radix> &);
        SU2<radix> & operator *= (const radix &);
        SU2<radix> & operator = (const SU2<radix> &);

    private:

        radix a, b, c, d;  

};

template<typename radix>
radix SU2<radix>::sqr(){
    return a*a + b*b + c*c + d*d;
}

template<typename radix>
SU2<radix> & SU2<radix>::normalize(){
    radix sq = sqrt(this->sqr());
    a /= sq;
    b /= sq;
    c /= sq;
    d /= sq;   
    return *this;
}

template<typename radix>
SU2<radix> & SU2<radix>::reunitarize(){
    return this->normalize();
}

template<typename radix>
SU2<radix> & SU2<radix>::random(){
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
SU2<radix> & SU2<radix>::inv(){
    a *= static_cast<regex>(-1);
    b *= static_cast<regex>(-1);
    c *= static_cast<regex>(-1);
    d *= static_cast<regex>(-1);
    return *this;
}

inline SU2<radix> SU2<radix>::operator * (const SU2<radix> & y){
    SU2<radix> r;
    r.a = nn_a(*this,y); 
    r.b = nn_b(*this,y); 
    r.c = nn_c(*this,y); 
    r.d = nn_d(*this,y); 
    return r;
};

inline SU2<radix> SU2<radix>::operator *= (const SU2<radix> & y){
    a = nn_a(*this,y); 
    b = nn_b(*this,y); 
    c = nn_c(*this,y); 
    d = nn_d(*this,y); 
    return *this;
};

SU2<radix> & SU2<radix>::operator *= (const radix & rhs){ 
    a *= static_cast<regex>(rhs);
    b *= static_cast<regex>(rhs);
    c *= static_cast<regex>(rhs);
    d *= static_cast<regex>(rhs);
    return *this;
};

