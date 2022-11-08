#ifndef REAL_VAR_OPS_H_
#define REAL_VAR_OPS_H_

#include "plumbing/defs.h"
#include "plumbing/random.h"


//////////////////////////////////////////////////////////////////////////

// define abs() -function for floating point
template <typename T, std::enable_if_t<std::is_floating_point<T>::value, int> = 0>
inline T abs(T val) {
    return std::abs(val);
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

// this defines gaussian_random(T & arg) for all double -> T convertible types
template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline void gaussian_random(T &d, T width = 1.0) {
    static_assert(hila::is_floating_point<T>::value,
                  "gaussian_random() requires floating point type argument");

    d = static_cast<T>(hila::gaussrand()) * width;
}

// this does the same for std. random
template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline void random(T &d) {
    static_assert(hila::is_floating_point<T>::value,
                  "random() requires floating point type argument");

    d = static_cast<T>(hila::random());
}

/// convert to string: separator does nothing, but for compatibility w. other to_strings

namespace hila {
template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
std::string to_string(const T v, int prec = 8, char separator = ' ') {
    std::stringstream ss;
    ss.precision(prec);
    ss << v;
    return ss.str();
}

template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
std::string prettyprint(const T v, int prec = 8) {
    return to_string(v,prec);
}


} // namespace hila

#endif