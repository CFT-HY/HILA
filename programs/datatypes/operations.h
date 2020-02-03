
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

template<typename T, std::enable_if_t<is_arithmetic<T>::value, int> = 0>
#pragma transformer loop_function
inline auto norm_sq(T val){
  return val*val;
}


#endif 