#ifndef SUNM
#define SUNM

#include "cmplx.h"
#include "general_matrix.h"
#include "sun_vector.h"
#include "../plumbing/defs.h"
#include "../plumbing/mersenne.h" //has to be included
#include <cmath>

// SU3 crossprod calculation routines 

#define CMULJJ(a,b,c) do { (c).re =  (a).re*(b).re - (a).im*(b).im; \
   		        (c).im = -(a).re*(b).im - (a).im*(b).re; } while(0)
#define CMULJJ_ADD(a,b,c) do { (c).re +=  (a).re*(b).re - (a).im*(b).im; \
   		            (c).im += -(a).re*(b).im - (a).im*(b).re; } while(0)
#define CMULJJ_SUB(a,b,c) do { (c).re -=  (a).re*(b).re - (a).im*(b).im; \
   	                (c).im += (a).re*(b).im + (a).im*(b).re; } while(0)

#define su3_vector_crossprod_conj(av, bv, res)  \
do {                                            \
                                                \
    CMULJJ((av).c[1], (bv).c[2], (res).c[0]);         \
    CMULJJ_SUB((av).c[2], (bv).c[1], (res).c[0]);     \
                                                \
    CMULJJ((av).c[2], (bv).c[0], (res).c[1]);         \
    CMULJJ_SUB((av).c[0], (bv).c[2], (res).c[1]);     \
                                                \
    CMULJJ((av).c[0], (bv).c[1], (res).c[2]);         \
    CMULJJ_SUB((av).c[1], (bv).c[0], (res).c[2]);     \
} while (0)

// SU2 matrix multiplication routines ------------------------

#define nn_a(x,y) (x.d*y.a + x.a*y.d - x.b*y.c + x.c*y.b)
#define nn_b(x,y) (x.d*y.b + x.b*y.d - x.c*y.a + x.a*y.c)
#define nn_c(x,y) (x.d*y.c + x.c*y.d - x.a*y.b + x.b*y.a)
#define nn_d(x,y) (x.d*y.d - x.a*y.a - x.b*y.b - x.c*y.c)

#define na_a(x,y) (-x.d*y.a + x.a*y.d + x.b*y.c - x.c*y.b)
#define na_b(x,y) (-x.d*y.b + x.b*y.d + x.c*y.a - x.a*y.c)
#define na_c(x,y) (-x.d*y.c + x.c*y.d + x.a*y.b - x.b*y.a)
#define na_d(x,y) ( x.d*y.d + x.a*y.a + x.b*y.b + x.c*y.c)

#define an_a(x,y) ( x.d*y.a - x.a*y.d + x.b*y.c - x.c*y.b)
#define an_b(x,y) ( x.d*y.b - x.b*y.d + x.c*y.a - x.a*y.c)
#define an_c(x,y) ( x.d*y.c - x.c*y.d + x.a*y.b - x.b*y.a)
#define an_d(x,y) ( x.d*y.d + x.a*y.a + x.b*y.b + x.c*y.c)

#define aa_a(x,y) (-x.d*y.a - x.a*y.d - x.b*y.c + x.c*y.b)
#define aa_b(x,y) (-x.d*y.b - x.b*y.d - x.c*y.a + x.a*y.c)
#define aa_c(x,y) (-x.d*y.c - x.c*y.d - x.a*y.b + x.b*y.a)
#define aa_d(x,y) ( x.d*y.d - x.a*y.a - x.b*y.b - x.c*y.c)

// gaussian rng generation routines ----------------------

#define VARIANCE 0.5
#define MYSINF(X) sin(X)
#define MYCOSF(X) cos(X)

template<typename radix>
radix gaussian_ran2 (radix* out2) 
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
/// SU(N) matrix class
/// Implementations are found in this header file, since non-specialized
/// templates have to be visible to the tranlation units that use them
//////////////////

template<int n, typename radix>
using SU = matrix<n,n,cmplx<radix> >;

template<typename radix> 
class SU2; 

template<typename radix> 
class SU3; 

template<typename radix> 
class SU2adjoint; 

template<typename radix> 
class SU3 : public matrix<3,3,cmplx<radix> > {
    public:
    //constructor from two Su vectors
    SU3 (const SU_vector<3, radix> & a, const SU_vector<3, radix> & b){
        SU_vector<3, radix> c; //last column of matrix to be made from cross product
        const SU_vector<3, radix> ai[3] = { a, b, c };
        int i,j;
        su3_vector_crossprod_conj(a,b,c);
        #ifdef CUDA
        #pragma unroll
        #endif
        for (i = 0; i < 3; i++){
        #ifdef CUDA
        #pragma unroll
        #endif
            for (j = 0; j < 3; j++){
                this->c[i][j] = ai[i].c[j];
            }
        }
    }
};

template<typename radix>
class SU2 { 
    public: 

        SU2() : a(0), b(0), c(0), d(1) {}

        SU2(radix * vals) : a(vals[0]), b(vals[1]), c(vals[2]), d(vals[3]) { normalize(); }

        SU2(const SU2vector<radix> & rhs){
            b = rhs.c[0].re;				
            a = rhs.c[0].im;				
            d = rhs.c[1].re;				
            c =-rhs.c[1].im;           
        };

        SU2(const SU2<radix> & rhs){
            b = rhs.b;				
            a = rhs.a;				
            d = rhs.d;				
            c = rhs.c;           
        };

        friend SU2<radix> operator * (const SU2<radix> & x, const SU2<radix> & y){
            SU2<radix> r;
            r.a = nn_a(x,y); r.b = nn_b(x,y); 
            r.c = nn_c(x,y); r.d = nn_d(x,y);
            return r;
        }
 
        friend SU2<radix> operator * (const SU2<radix> & x, const SU2adjoint<radix> & y){
            SU2<radix> r;
            r.a = na_a(x,y.ref); r.b = na_b(x,y.ref); \
            r.c = na_c(x,y.ref); r.d = na_d(x,y.ref);
            return r;
        }
 
        friend SU2<radix> operator * (const SU2adjoint<radix> & x, const SU2<radix> & y){
            SU2<radix> r;
            r.a = an_a(x.ref,y); r.b = an_b(x.ref,y); \
            r.c = an_c(x.ref,y); r.d = an_d(x.ref,y);
            return r;
        }

        friend SU2<radix> operator * (const SU2adjoint<radix> & x, const SU2adjoint<radix> & y){
            SU2<radix> r;
            r.a = aa_a(x.ref, y.ref); r.b = aa_b(x.ref, y.ref); \
            r.c = aa_c(x.ref, y.ref); r.d = aa_d(x.ref, y.ref);
            return r;
        }

        SU2<radix> & operator = (const SU2<radix> &);
        SU2<radix> & operator = (const SU2adjoint<radix> &);
        SU2<radix> & operator = (const SU2vector<radix> &);
        SU2<radix> & normalize();
        SU2<radix> & reunitarize();  
        SU2<radix> & random(); 
        SU2<radix> & inv(); 
        SU2adjoint<radix> & adj(); 
        
        radix sqr() const;
        radix tr() const;
        radix det() const; 

        SU2<radix> operator + (const SU2<radix> &); //basic operations. These versions return new matrix as result
        SU2<radix> operator - (const SU2<radix> &);
        SU2<radix> & operator += (const SU2<radix> &); //same ops as above, except store result in lhs matrix
        SU2<radix> & operator -= (const SU2<radix> &);
        SU2<radix> & operator *= (const SU2<radix> &);

        SU2<radix> & operator += (const SU2adjoint<radix> &); //same ops as above, except store result in lhs matrix
        SU2<radix> & operator -= (const SU2adjoint<radix> &);
        SU2<radix> & operator *= (const SU2adjoint<radix> &);

        SU2<radix> & operator *= (const radix &);
        SU2<radix> operator * (const radix &);
    private:
        radix a, b, c, d;  
};


template<typename radix>
class SU2adjoint {
    public:
    SU2adjoint (const SU2<radix> & rhs) : ref(rhs) {} ;
    const SU2<radix> & ref;
};

template<typename radix>
SU2adjoint<radix> & SU2<radix>::adj(){
    return SU2adjoint<radix>(*this);
}; 

template<typename radix>
radix SU2<radix>::sqr() const {
    return a*a + b*b + c*c + d*d;
}

template<typename radix>
radix SU2<radix>::det() const {
    return a*a + b*b + c*c + d*d;
}

template<typename radix>
radix SU2<radix>::tr() const {
    return 2*d;
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
    a *= static_cast<radix>(-1);
    b *= static_cast<radix>(-1);
    c *= static_cast<radix>(-1);
    d *= static_cast<radix>(-1);
    return *this;
}

template<typename radix>
SU2<radix> & SU2<radix>::operator = (const SU2<radix> & rhs){
    a = rhs.a;
    b = rhs.b;
    c = rhs.c;
    d = rhs.d;
    return *this;
};

template<typename radix>
SU2<radix> & SU2<radix>::operator = (const SU2adjoint<radix> & rhs){
    a = -rhs.a;
    b = -rhs.b;
    c = -rhs.c;
    d = rhs.d;
    return *this;
};

template<typename radix>
SU2<radix> & SU2<radix>::operator *= (const SU2<radix> & y){
    a = nn_a((*this),y); 
    b = nn_b((*this),y); 
    c = nn_c((*this),y); 
    d = nn_d((*this),y); 
    return *this;
}

template<typename radix>
SU2<radix> & SU2<radix>::operator *= (const SU2adjoint<radix> & y){
    a = na_a((*this),y); 
    b = na_b((*this),y); 
    c = na_c((*this),y); 
    d = na_d((*this),y); 
    return *this;
}

template<typename radix>
SU2<radix> & SU2<radix>::operator *= (const radix & rhs){ 
    a *= rhs;
    b *= rhs;
    c *= rhs;
    d *= rhs;
    return *this;
};

template<typename radix>
SU2<radix> SU2<radix>::operator * (const radix & rhs){ 
    SU2<radix> r;
    r.a = a*rhs;
    r.b = b*rhs;
    r.c = c*rhs;
    r.d = d*rhs;
    return r;
};

template<typename radix>
SU2<radix> SU2<radix>::operator + (const SU2<radix> & y){
    SU2<radix> r;
    r.a = a + y.a; r.b = b + y.b; \
    r.c = c + y.c; r.d = d + y.d;
    return r;
};

template<typename radix>
SU2<radix> SU2<radix>::operator - (const SU2<radix> & y){
    SU2<radix> r;
    r.a = a - y.a; r.b = b - y.b; \
    r.c = c - y.c; r.d = d - y.d;
    return r;
};

#endif 