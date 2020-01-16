
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
inline double norm_sq(T val){
  return static_cast<double>(val*val);
}

#endif 