
#ifndef OPS_H
#define OPS_H

///////////////////////////////////////////
// The routines are used as operation modifiers
// i.e. they change how +,-,* are performed
// between members of the same data type.   
// By default they do not change any behaviour 
// unless they are specialized for a dtype. 
// For example A*trans(B) multiplies A by the
// transpose of the elements of B when they are 
// matrices. 
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
  return val*val;
}

#endif 