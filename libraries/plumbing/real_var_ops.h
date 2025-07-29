#ifndef HILA_REAL_VAR_OPS_H_
#define HILA_REAL_VAR_OPS_H_

#include "plumbing/defs.h"
#include "plumbing/random.h"


//////////////////////////////////////////////////////////////////////////

// define abs() -function for floating point
template <typename T, std::enable_if_t<std::is_floating_point<T>::value, int> = 0>
inline T abs(T val) {
    return std::abs(val);
}

// also def min and max, for all arith types
template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
inline T min(T val1, T val2) {
    return std::min(val1, val2);
}

// also def min and max
template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
inline T max(T val1, T val2) {
    return std::max(val1, val2);
}

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

template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline T dagger(T val) {
    return val;
}

// define squarenorm and norm
template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline T squarenorm(T val) {
    return val * val;
}

// for integral types, sqrt -> double, so auto here makes sens
template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline auto norm(T val) {
    return sqrt(val * val);
}

namespace hila {


/// convert to string: separator does nothing, but for compatibility w. other to_strings

template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
std::string to_string(const T v, int prec = 8, char separator = ' ') {
    std::stringstream ss;
    ss.precision(prec);
    ss << v;
    return ss.str();
}

template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
std::string prettyprint(const T v, int prec = 8) {
    return to_string(v, prec);
}


} // namespace hila

#endif