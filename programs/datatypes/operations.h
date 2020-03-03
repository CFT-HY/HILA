
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
inline T norm_squared(T val){
  return val.norm_sq();
}

inline double norm_squared(double val){
  return val*val;
}

inline float norm_squared(float val){
  return val*val;
}



#endif 