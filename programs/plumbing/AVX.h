#ifndef AVX_H
#define AVX_H

#include <immintrin.h>


/// A new vector class is necessary, the intrinsic types
/// are not always proper types. This encapsulates them
/// and defines basic arithmetic
struct avxdvector {
  __m256d c;

  avxdvector() = default;
  avxdvector(const avxdvector & a) =default;

  constexpr avxdvector(__m256d x):c(x) {}
  avxdvector(double x):c(_mm256_broadcast_sd(&x)) {}

  // Cast to base type interpred as a sum, implements
  // the sum reduction
  #pragma transformer loop_function
  operator double() { 
    __m256d s = _mm256_hadd_pd(c,c);
    return ((double*)&s)[0] + ((double*)&s)[2];
  }

  #pragma transformer loop_function
  avxdvector & operator+= (const avxdvector & lhs) {
    c += lhs.c;
    return *this;
  }
};


/* Define operations for the vector type */

#pragma transformer loop_function
avxdvector operator+(const avxdvector & a, const avxdvector & b) {
  return avxdvector(_mm256_add_pd(a.c, b.c));
}


#pragma transformer loop_function
double operator+(double & a, const avxdvector & b) {
  __m256d s = _mm256_hadd_pd(b.c,b.c);
  return a + ((double*)&s)[0] + ((double*)&s)[2];
}


#endif