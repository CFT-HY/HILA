#ifndef AVX_H
#define AVX_H

#include <immintrin.h>
#include "../vectorclass/vectorclass.h"


#define VECTORIZED
constexpr static int max_vector_size = 8;


#pragma transformer loop_function
inline double reduce_sum(Vec4d v){
  float sum = 0;
  for(int i=0; i<4; i++)
    sum += v[i];
  return sum;
}

#pragma transformer loop_function
inline double reduce_sum(Vec8f v){
  float sum = 0;
  for(int i=0; i<8; i++)
    sum += v[i];
  return sum;
}



namespace vectorized
{
  template<int vector_len>
  inline void permute(int *perm, void * element, int n_elements){}

  template<>
  inline void permute<4>(int *perm, void * element, int n_elements){
    __m256d * e = (__m256d *) element;
    for( int v=0; v<n_elements; v++ ){
      __m256d t = e[v];
      for( int i=0; i<4; i++ )
        e[v][i] =  t[perm[i]];
    }
  }

  // Vector length in bytes
  constexpr int sizeofvector = 32;
}


#endif