#ifndef AVX_H
#define AVX_H

#include <immintrin.h>

#define VECTORIZED
constexpr int max_vector_size = 8;


/// A new vector class is necessary, the intrinsic types
/// are not always proper types. This encapsulates them
/// and defines basic arithmetic
struct avxdvector {
  __m256d c;

  avxdvector() = default;
  avxdvector(const avxdvector & a) =default;

  constexpr avxdvector(__m256d x):c(x) {}
  avxdvector(const double & x):c(_mm256_broadcast_sd(&x)) {}

  // Cast to base type interpred as a sum, implements
  // the sum reduction
  #pragma transformer loop_function
  avxdvector & operator+= (const avxdvector & lhs) {
    c += lhs.c;
    return *this;
  }

  #pragma transformer loop_function
  double reduce_sum(){
    __m256d s = _mm256_hadd_pd(c,c);
    return ((double*)&s)[0] + ((double*)&s)[2];
  }

  #pragma transformer loop_function
  double reduce_prod(){
    double m = 1;
    for(int i=0; i<4; i++){
      m*=((double*)&c)[i];
    }
    return m;
  }

  avxdvector operator-() const {return _mm256_xor_pd(c, _mm256_set1_pd(-0.0)); }

};



/* Define operations for the vector type */

#pragma transformer loop_function
inline avxdvector operator+(const avxdvector & a, const avxdvector & b) {
  return avxdvector(_mm256_add_pd(a.c, b.c));
}

#pragma transformer loop_function
inline avxdvector operator-(const avxdvector & a, const avxdvector & b) {
  return avxdvector(_mm256_sub_pd(a.c, b.c));
}

#pragma transformer loop_function
inline avxdvector operator*(const avxdvector & a, const avxdvector & b) {
  return avxdvector(_mm256_mul_pd(a.c, b.c));
}

#pragma transformer loop_function
inline avxdvector operator/(const avxdvector & a, const avxdvector & b) {
  return avxdvector(_mm256_div_pd(a.c, b.c));
}



/// A new vector class is necessary, the intrinsic types
/// are not always proper types. This encapsulates them
/// and defines basic arithmetic
struct avxfvector {
  __m256 c;

  avxfvector() = default;
  avxfvector(const avxfvector & a) =default;

  constexpr avxfvector(__m256 x):c(x) {}
  avxfvector(const float & x):c(_mm256_broadcast_ss(&x)) {}

  // Cast to base type interpred as a sum, implements
  // the sum reduction
  #pragma transformer loop_function
  avxfvector & operator+= (const avxfvector & lhs) {
    c += lhs.c;
    return *this;
  }

  #pragma transformer loop_function
  float reduce_sum(){
    float sum = 0;
    for(int i=0; i<8; i++)
      sum += c[i];
    return sum;
  }

  #pragma transformer loop_function
  float reduce_prod(){
    float m = 1;
    for(int i=0; i<4; i++){
      m*=((float*)&c)[i];
    }
    return m;
  }

  avxfvector operator-() const {return _mm256_xor_ps(c, _mm256_set1_ps(-0.0)); }

};


/* Define operations for the vector type */

#pragma transformer loop_function
inline avxfvector operator+(const avxfvector & a, const avxfvector & b) {
  return avxfvector(_mm256_add_ps(a.c, b.c));
}

#pragma transformer loop_function
inline avxfvector operator-(const avxfvector & a, const avxfvector & b) {
  return avxfvector(_mm256_sub_ps(a.c, b.c));
}

#pragma transformer loop_function
inline avxfvector operator*(const avxfvector & a, const avxfvector & b) {
  return avxfvector(_mm256_mul_ps(a.c, b.c));
}

#pragma transformer loop_function
inline avxfvector operator/(const avxfvector & a, const avxfvector & b) {
  return avxfvector(_mm256_div_ps(a.c, b.c));
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