#ifndef SU_VEC
#define SU_VEC

#include "cmplx.h"
#include "array.h"

template<int N, typename radix> 
class SU_vector : public array<N,cmplx<radix>> {};

#endif