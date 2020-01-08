#ifndef SU_VEC
#define SU_VEC

#include "cmplx.h"

template<int N, typename radix>
class SU_vector {
    public:
    cmplx<radix> c[N];
};

template<typename radix> 
class SU2; 

template<typename radix> 
class SU2vector {
    public:
    cmplx<radix> c[2];
    SU2vector(const SU2<radix> & m){
        c[0].re = m.b;				    
        c[0].im = m.a;			       	
        c[1].re = m.d;			       	
        c[1].im =-m.c;
    }
};

#endif