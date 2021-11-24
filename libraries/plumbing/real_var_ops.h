#ifndef REAL_VAR_OPS_H_
#define REAL_VAR_OPS_H_

#include "plumbing/defs.h"
#include "plumbing/random.h"



//////////////////////////////////////////////////////////////////////////

// define real(), imag(), conj() -functions for basic arithmetic types
template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline T real(T val) {
    return val;
}

template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline T imag(T val) {
    return (T)0;
}

template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline T conj(T val) {
    return val;
}

// define squarenorm and norm
template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline T squarenorm(T val) {
    return val * val;
}

template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline T norm(T val) {
    return sqrt(val * val);
}

// this defines gaussian_random(T & arg) for all double -> T convertible types
template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline void gaussian_random(T &d, T width = 1.0) {
    d = static_cast<T>(hila::gaussrand()) * width;
}

// this does the same for std. random
template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline void random(T &d) {
    d = static_cast<T>(hila::random());
}

#endif