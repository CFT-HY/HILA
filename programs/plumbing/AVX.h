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


inline double reduce_sum(Vec4d v){
  double sum = 0;
  double store[4];
  v.store(&(store[0]));
  for(int i=0; i<4; i++)
    sum += store[i];
  return sum;
}

inline double reduce_sum(Vec8f v){
  double sum = 0;
  float store[8];
  v.store(&(store[0]));
  for(int i=0; i<8; i++)
    sum += store[i];
  return sum;
}

inline double reduce_sum(Vec8i v){
  double sum = 0;
  int store[8];
  v.store(&(store[0]));
  for(int i=0; i<8; i++)
    sum += store[i];
  return sum;
}

inline double reduce_prod(Vec4d v){
  double sum = 1;
  double store[4];
  v.store(&(store[0]));
  for(int i=0; i<4; i++)
    sum *= store[i];
  return sum;
}

inline double reduce_prod(Vec8f v){
  double sum = 0;
  float store[8];
  v.store(&(store[0]));
  for(int i=0; i<8; i++)
    sum *= store[i];
  return sum;
}

inline double reduce_prod(Vec8i v){
  double sum = 0;
  int store[8];
  v.store(&(store[0]));
  for(int i=0; i<8; i++)
    sum *= store[i];
  return sum;
}


/// Define modulo operator for integer vector
inline Vec8i operator%( const Vec8i &lhs, const int &rhs)
{
  Vec8i r;
  int tvec1[8], tvec2[8];
  lhs.store(&(tvec1[0]));
  for(int i=0; i<8; i++)
    tvec2[i] = tvec1[i] % rhs;
  r.load(&(tvec2[0]));
  return r;
}

inline Vec4i operator%( const Vec4i &lhs, const int &rhs)
{
  Vec4i r;
  int tvec1[4], tvec2[4];
  lhs.store(&(tvec1[0]));
  for(int i=0; i<4; i++)
    tvec2[i] = tvec1[i] % rhs;
  r.load(&(tvec2[0]));
  return r;
}


inline Vec4d hila_random_Vec4d(){
  Vec4d r;
  double tvec[4];
  for(int i=0; i<4; i++){
    tvec[i] = mersenne();
  }
  r.load(&(tvec[0]));
  return r;
}

inline Vec8f hila_random_Vec8f(){
  Vec8f r;
  float tvec[8];
  for(int i=0; i<8; i++){
    tvec[i] = mersenne();
  }
  r.load(&(tvec[0]));
  return r;
}



namespace vectorized
{
  template<int vector_len>
  inline void permute(int *perm, void * element, const int n_elements){
    assert(false && "permute not implemented for vector size");
  }

  template<>
  inline void permute<4>(int *perm, void * element, const int n_elements){
    Vec4d * e = (Vec4d *) element;
    for( int v=0; v<n_elements; v++ ){
      double t1[4], t2[4];
      e[v].store(&(t1[0]));
      for( int i=0; i<4; i++ )
        t2[i] =  t1[perm[i]];
      e[v].load(&(t2[0]));
    }
  }

  template<>
  inline void permute<8>(int *perm, void * element, const  int n_elements){
    Vec8f * e = (Vec8f *) element;
    for( int v=0; v<n_elements; v++ ){
      float t1[8], t2[8];
      e[v].store(&(t1[0]));
      for( int i=0; i<8; i++ )
        t2[i] =  t1[perm[i]];
      e[v].load(&(t2[0]));
    }
  }

}


#endif