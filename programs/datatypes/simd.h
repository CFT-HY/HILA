#ifndef SIMD_H
#define SIMD_H

#include <immintrin.h>

/// Here simd operations: AVX 2 and AVX512

const int SIMDV_SIZE = 256;   // AVX2, do this in makefile

template <typename T>
using 

template <typename T>
class simdv {
public:
  enum { num_floats = SIMDV_SIZE/32, num_doubles = SIMDV_SIZE/64 };
  
  __m256 

  simdv() 
} 




#endif // SIMD_H