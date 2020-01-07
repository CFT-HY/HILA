
#ifndef OPS_H
#define OPS_H

///////////////////////////////////////////
// These routines define the default
// behaviour of trans, conj, etc.
// routines. These routines alter
// how operations are performed between
// members of the same dtype. If they are
// not specialized for a dtype they do not do anything.
// Specializations of these templates are
// found in the headers defining each dtype.
////////////////////////////////////////////

template<typename T>
#pragma transformer loop_function
inline T conj(T rhs){
  return rhs;
}

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