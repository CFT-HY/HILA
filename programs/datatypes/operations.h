
#ifndef OPS_H
#define OPS_H

///////////////////////////////////////////
// General operation definitions for dtypes
////////////////////////////////////////////

//complex conjugate
template<typename T>
#pragma transformer loop_function
inline T conj(T rhs){
  return rhs;
}

//transpose 
template<typename T>
#pragma transformer loop_function
inline T trans(T rhs){
  return rhs;
}

template<typename T>
#pragma transformer loop_function
inline auto type_norm_sq(T val){
  return val*val;
}

template<typename T> 
#pragma transformer loop_function
inline auto type_norm_sq(cmplx<T> val){
  return val.norm_sq();
}

#endif 