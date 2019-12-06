#ifndef AVX_H
#define AVX_H

#include <immintrin.h>


/// A new vector class is necessary, the intrinsic types
/// are not always proper types. This encapsulates them
/// and defines basic arithmetic
struct avxdvector {
  __m256d c;

  constexpr avxdvector(__m256d x):c(x) {}
  
};


/* Define operations for the vector type */

#pragma transformer loop_function
avxdvector operator+(const avxdvector & a, const avxdvector & b) {
  return avxdvector(_mm256_add_pd(a.c, b.c));
}




#endif