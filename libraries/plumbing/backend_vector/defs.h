#ifndef VECTOR_DEFS_H
#define VECTOR_DEFS_H

#include "plumbing/defs.h"
#include <immintrin.h>
#include "../../vectorclass/vectorclass.h"
#include "../../vectorclass/vectormath_exp.h"
#include "../../vectorclass/vectormath_trig.h"
#include "../../vectorclass/vectormath_hyp.h"

#define VECTORIZED

#ifndef VECTOR_SIZE
#define VECTOR_SIZE 32
#endif

// Define random number generator
inline double hila_random(){ return mersenne(); }
inline void seed_random(int seed) {

  /* First seed the generator */
  seed_mersenne(seed);

  /* "Warm up" to create a full state */
	for(int i=0; i<543210; i++)
    mersenne();
}


// Trivial synchronization
inline void synchronize_threads(){}



/// Implements test for basic in types, similar to 
/// std::is_arithmetic, but allows the backend to add
/// it's own basic types (such as AVX vectors)
template< class T >
struct is_avx_vector : std::integral_constant<
  bool,
  std::is_same<T,Vec4d>::value ||
  std::is_same<T,Vec8f>::value ||
  std::is_same<T,Vec8i>::value ||
  std::is_same<T,Vec8d>::value ||
  std::is_same<T,Vec16f>::value ||
  std::is_same<T,Vec16i>::value 
> {};

template< class T >
struct is_arithmetic : std::integral_constant<
  bool,
  std::is_arithmetic<T>::value ||
  is_avx_vector<T>::value
> {};






/*** The next section contains basic operations for vectors ***/

// Norm squared
inline Vec4d norm_squared(Vec4d val){
  return val*val;
}

inline Vec8f norm_squared(Vec8f val){
  return val*val;
}

inline Vec8i norm_squared(Vec8i val){
  return val*val;
}

inline Vec8d norm_squared(Vec8d val){
  return val*val;
}

inline Vec16f norm_squared(Vec16f val){
  return val*val;
}

inline Vec16i norm_squared(Vec16i val){
  return val*val;
}

// Reductions
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

inline int64_t reduce_sum(Vec8i v){
  int64_t sum = 0;
  int store[8];
  v.store(&(store[0]));
  for(int i=0; i<8; i++)
    sum += store[i];
  return sum;
}

inline double reduce_sum(Vec8d v){
  double sum = 0;
  double store[8];
  v.store(&(store[0]));
  for(int i=0; i<8; i++)
    sum += store[i];
  return sum;
}

inline double reduce_sum(Vec16f v){
  double sum = 0;
  float store[16];
  v.store(&(store[0]));
  for(int i=0; i<16; i++)
    sum += store[i];
  return sum;
}

inline int64_t reduce_sum(Vec16i v){
  int64_t sum = 0;
  int store[16];
  v.store(&(store[0]));
  for(int i=0; i<16; i++)
    sum += store[i];
  return sum;
}

inline int64_t reduce_sum(Vec4q v){
  int64_t sum = 0;
  int64_t store[4];
  v.store(&(store[0]));
  for(int i=0; i<4; i++)
    sum += store[i];
  return sum;
}

inline int64_t reduce_sum(Vec8q v){
  int64_t sum = 0;
  int64_t store[4];
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

inline double reduce_prod(Vec8d v){
  double sum = 1;
  double store[8];
  v.store(&(store[0]));
  for(int i=0; i<8; i++)
    sum *= store[i];
  return sum;
}

inline double reduce_prod(Vec16f v){
  double sum = 0;
  float store[16];
  v.store(&(store[0]));
  for(int i=0; i<16; i++)
    sum *= store[i];
  return sum;
}

inline double reduce_prod(Vec16i v){
  double sum = 0;
  int store[16];
  v.store(&(store[0]));
  for(int i=0; i<16; i++)
    sum *= store[i];
  return sum;
}

inline double reduce_prod(Vec4q v){
  double sum = 0;
  int64_t store[4];
  v.store(&(store[0]));
  for(int i=0; i<4; i++)
    sum *= store[i];
  return sum;
}

inline double reduce_prod(Vec8q v){
  double sum = 0;
  int64_t store[8];
  v.store(&(store[0]));
  for(int i=0; i<8; i++)
    sum *= store[i];
  return sum;
}



// Return the 
template <typename base_t, typename vector_t, typename T, typename vecT>
T reduce_sum_in_vector(const vecT & vt) {
  constexpr int nvec = sizeof(vecT)/sizeof(vector_t);
  static_assert( nvec == sizeof(T)/sizeof(base_t), "Mismatch in vectorized type sizes");
  T res;
  auto * vptr = (const vector_t *)(&vt);
  base_t * bptr = (base_t *)(&res);
  for (int i=0; i<nvec; i++) {
    bptr[i] = reduce_sum(vptr[i]);
  }
  return res;
}



/// If vector elements are implemented in the c++ code,
/// reductions to base variables need to be supported.
/// Will this lead to problematic behavior?
template<typename Vec>
inline double& operator+=(double &lhs, const Vec rhs)
{
  lhs += reduce_prod(rhs);
  return lhs;
}

template<typename Vec>
inline float& operator+=(float &lhs, const Vec rhs)
{
  lhs += reduce_prod(rhs);
  return lhs;
}


// Define modulo operator for integer vector
inline Vec16i operator%( const Vec16i &lhs, const int &rhs)
{
  Vec16i r;
  int tvec1[16], tvec2[16];
  lhs.store(&(tvec1[0]));
  for(int i=0; i<16; i++)
    tvec2[i] = tvec1[i] % rhs;
  r.load(&(tvec2[0]));
  return r;
}


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


// Random numbers
// Since you cannot specialize by return type,
// it needs to be a struct...
#if VECTOR_SIZE == 32
template<typename T>
inline auto hila_random_vector(){
  Vec8f r;
  float tvec[8];
  for(int i=0; i<8; i++){
    tvec[i] = mersenne();
  }
  r.load(&(tvec[0]));
  return r;
};

template<>
inline auto hila_random_vector<double>(){
  Vec4d r;
  double tvec[4];
  for(int i=0; i<4; i++){
    tvec[i] = mersenne();
  }
  r.load(&(tvec[0]));
  return r;
};


#elif VECTOR_SIZE == 64

template<typename T>
inline auto hila_random_vector(){
  Vec16f r;
  float tvec[16];
  for(int i=0; i<16; i++){
    tvec[i] = mersenne();
  }
  r.load(&(tvec[0]));
  return r;
};

template<>
inline auto hila_random_vector<double>(){
  Vec8d r;
  double tvec[8];
  for(int i=0; i<8; i++){
    tvec[i] = mersenne();
  }
  r.load(&(tvec[0]));
  return r;
};

#endif





#endif
