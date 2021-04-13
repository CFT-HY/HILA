
#ifndef OPS_H
#define OPS_H

///////////////////////////////////////////
// General operation definitions for dtypes
////////////////////////////////////////////

// complex conjugate
//template <typename T> inline T conj(T rhs) { return rhs; }

// transpose
// template <typename T> inline T trans(T rhs) { return rhs; }

// template <typename T> inline auto squarenorm(T val) { return val.squarenorm(); }

inline double squarenorm(double val) { return val * val; }

inline float squarenorm(float val) { return val * val; }

#endif