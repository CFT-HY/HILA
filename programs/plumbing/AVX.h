#ifndef AVX_H
#define AVX_H

#include "../plumbing/defs.h"
#include <immintrin.h>
#include "../vectorclass/vectorclass.h"
#include "../vectorclass/vectormath_exp.h"
#include "../vectorclass/vectormath_trig.h"
#include "../vectorclass/vectormath_hyp.h"


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

#pragma transformer loop_function
inline double reduce_sum(Vec8i v){
  float sum = 0;
  for(int i=0; i<8; i++)
    sum += v[i];
  return sum;
}

#pragma transformer loop_function
inline double reduce_prod(Vec4d v){
  float sum = 1;
  for(int i=0; i<4; i++)
    sum *= v[i];
  return sum;
}

#pragma transformer loop_function
inline double reduce_prod(Vec8f v){
  float sum = 1;
  for(int i=0; i<8; i++)
    sum *= v[i];
  return sum;
}

#pragma transformer loop_function
inline double reduce_prod(Vec8i v){
  float sum = 1;
  for(int i=0; i<8; i++)
    sum *= v[i];
  return sum;
}


/// Define modulo operator for integer vector
#pragma transformer loop_function
inline Vec8i operator%( const Vec8i &lhs, const int &rhs)
{
  Vec8i r;
  for(int i=0; i<8; i++)
    r.insert(i, lhs[i] % rhs);
  return r;
}


#pragma transformer loop_function
inline Vec4d hila_random_Vec4d(){
  Vec4d r;
  for(int i=0; i<4; i++){
    r.insert(i,mersenne());
  }
  return r;
}

#pragma transformer loop_function
inline Vec8f hila_random_Vec8f(){
  Vec8f r;
  for(int i=0; i<8; i++){
    r.insert(i,mersenne());
  }
  return r;
}



namespace vectorized
{
  template<int vector_len>
  inline void permute(int *perm, void * element, int n_elements){
    assert(false && "permute not implemented for vector size");
  }

  template<>
  inline void permute<4>(int *perm, void * element, int n_elements){
    __m256d * e = (__m256d *) element;
    for( int v=0; v<n_elements; v++ ){
      __m256d t = e[v];
      for( int i=0; i<4; i++ )
        e[v][i] =  t[perm[i]];
    }
  }

  template<>
  inline void permute<8>(int *perm, void * element, int n_elements){
    __m256 * e = (__m256 *) element;
    for( int v=0; v<n_elements; v++ ){
      __m256 t = e[v];
      for( int i=0; i<8; i++ )
        e[v][i] =  t[perm[i]];
    }
  }

  // Vector length in bytes
  constexpr int sizeofvector = 32;
}


#endif