#ifndef AVX_H
#define AVX_H

#include "../plumbing/defs.h"
#include <immintrin.h>
#include "../vectorclass/vectorclass.h"
#include "../vectorclass/vectormath_exp.h"
#include "../vectorclass/vectormath_trig.h"
#include "../vectorclass/vectormath_hyp.h"

#include "../plumbing/memory.h"


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

inline double reduce_sum(Vec16i v){
  double sum = 0;
  int store[16];
  v.store(&(store[0]));
  for(int i=0; i<16; i++)
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


/// Define modulo operator for integer vector
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


inline Vec8d hila_random_Vec8d(){
  Vec8d r;
  double tvec[8];
  for(int i=0; i<8; i++){
    tvec[i] = mersenne();
  }
  r.load(&(tvec[0]));
  return r;
}

inline Vec16f hila_random_Vec16f(){
  Vec16f r;
  float tvec[16];
  for(int i=0; i<16; i++){
    tvec[i] = mersenne();
  }
  r.load(&(tvec[0]));
  return r;
}








/// Utility for returning mapping a field element type into 
/// a corresponding vector. This is not used directly as a type
template <typename T>
struct field_info{
  constexpr static int vector_size = 1;
  constexpr static int base_type_size = 1;
  constexpr static int elements = 1;

  using base_type = double;
  using vector_type = Vec4d;
};



template <typename T>
class field_storage {
  public:

    // Array of structures implementation
    T * fieldbuf = NULL;

    void allocate_field( const int field_alloc_size ) {
      fieldbuf = (T*) allocate_field_mem( sizeof(T) * field_alloc_size);
      #pragma acc enter data create(fieldbuf)
    }

    void free_field() {
      #pragma acc exit data delete(fieldbuf)
      free_field_mem((void *)fieldbuf);
      fieldbuf = nullptr;
    }

    #pragma transformer loop_function
    inline T get(const int i, const int field_alloc_size) const
    {
      // There is some problem with directly assigning intrinsic vectors, at least.
      // This is a universal workaround, but could be fixed by assigning element
      // by element
      using vectortype = typename field_info<T>::vector_type;
      using basetype = typename field_info<T>::base_type;
      T value;
      basetype *vp = (basetype *) (fieldbuf + i);
      vectortype *valuep = (vectortype *)(&value);
      for( int e=0; e<field_info<T>::elements; e++ ){
        valuep[e].load(vp+e*field_info<T>::vector_size);
      }
      //std::memcpy( &value, fieldbuf+i, sizeof(T) );
      return value;
    }

    #pragma transformer loop_function
    inline void set(const T value, const int i, const int field_alloc_size) 
    {
      using vectortype = typename field_info<T>::vector_type;
      using basetype = typename field_info<T>::base_type;
      basetype *vp = (basetype *) (fieldbuf + i);
      vectortype *valuep = (vectortype *)(&value);
      for( int e=0; e<field_info<T>::elements; e++ ){
        valuep[e].store((vp + e*field_info<T>::vector_size));
      }
      //std::memcpy( fieldbuf+i, &value, sizeof(T) );
    }
};



#endif
