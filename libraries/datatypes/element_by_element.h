#ifndef HILA_ELEMENT_BY_ELEMENT
#define HILA_ELEMENT_BY_ELEMENT

/**
 * @file element_by_element.h
 * @brief easy element-by-element operations for matrices and vectors, working also with bare
 * numbers
 * @details Provides functions like
 *    hila::elem::exp(m)
 * which does element-by-element exponentiation of m.
 * This also works with scalars and complex numbers as expected, useful for generic template
 * programming.  Thus, hila::elem::exp(3.2) == exp(3.2)
 *
 * For example, if m is of type SquareMatrix<3,double>, then
 *     sqr(m)   // == m*m  matrix product
 *     elements::sqr(m)   // gives matrix with e(i,j) = sqr(m.e(i.j))
 *
 * Function hila::elem::add(a,b)  does element-by-element addition of same shape
 * Matrix/Vector/Array.  If one of the arguments is scalar or Complex, this is added to all
 * elements.  For matrices this is equivalent to (a.asArray() + b.asArray()).asMatrix(),
 * except that for hila::elem::add() the return value is the same as the input matrix type
 * if possible.
 *
 * Provides functions sqr, sqrt, cbrt, exp, log, sin, cos, tan, asin, acos, atan,
 *                    sinh, cosh, tanh, asinh, acosh, atanh,
 *                    pow,   // with real/complex power
 *                    add, sub, mul, div
 *
 */


#include "datatypes/matrix.h"
#include "datatypes/array.h"

namespace hila {
namespace elem {

/// @internal Helper routines to_array and back, for matrix/vector and other types

template <typename T, std::enable_if_t<hila::is_matrix<T>::value, int> = 0>
auto &to_array(const T &m) {
    return m.asArray();
}

template <typename T, typename A, std::enable_if_t<hila::is_matrix<T>::value, int> = 0>
const T &back(const A &m) {
    return m.asMatrix();
}


template <
    typename T,
    std::enable_if_t<hila::is_complex_or_arithmetic<T>::value || hila::is_array<T>::value, int> = 0>
const T &to_array(const T &m) {
    return m;
}

template <
    typename T,
    std::enable_if_t<hila::is_complex_or_arithmetic<T>::value || hila::is_array<T>::value, int> = 0>
const T &back(const T &m) {
    return m;
}

////////////////////////////////////////////////

/// @brief element-by-element square
template <typename T>
inline T sqr(const T &m) {
    return back<T>(::sqr(to_array(m)));
}

/// @brief e-by-e sqrt

template <typename T>
inline T sqrt(const T &m) {
    return back<T>(::sqrt(to_array(m)));
}

/// @brief cube root

template <typename T>
inline T cbrt(const T &m) {
    return back<T>(::cbrt(to_array(m)));
}

/// @brief exp

template <typename T>
inline T exp(const T &m) {
    return back<T>(::exp(to_array(m)));
}

/// @brief log

template <typename T>
inline T log(const T &m) {
    return back<T>(::log(to_array(m)));
}

/// @brief sin

template <typename T>
inline T sin(const T &m) {
    return back<T>(::sin(to_array(m)));
}

/// @brief cos

template <typename T>
inline T cos(const T &m) {
    return back<T>(::cos(to_array(m)));
}

/// @brief tan

template <typename T>
inline T tan(const T &m) {
    return back<T>(::tan(to_array(m)));
}

/// @brief asin

template <typename T>
inline T asin(const T &m) {
    return back<T>(::asin(to_array(m)));
}

/// @brief acos

template <typename T>
inline T acos(const T &m) {
    return back<T>(::acos(to_array(m)));
}

/// @brief atan

template <typename T>
inline T atan(const T &m) {
    return back<T>(::atan(to_array(m)));
}

/// @brief sinh

template <typename T>
inline T sinh(const T &m) {
    return back<T>(::sinh(to_array(m)));
}

/// @brief cosh

template <typename T>
inline T cosh(const T &m) {
    return back<T>(::cosh(to_array(m)));
}

/// @brief tanh

template <typename T>
inline T tanh(const T &m) {
    return back<T>(::tanh(to_array(m)));
}

/// @brief asinh

template <typename T>
inline T asinh(const T &m) {
    return back<T>(::asinh(to_array(m)));
}

/// @brief acosh

template <typename T>
inline T acosh(const T &m) {
    return back<T>(::acosh(to_array(m)));
}

/// @brief atanh

template <typename T>
inline T atanh(const T &m) {
    return back<T>(::atanh(to_array(m)));
}

/// @brief Power pow(m,p) -- accept only scalar powers in this construct
/// @details If more complicated power needed use Arrays

// template <typename T, typename P,
//           std::enable_if_t<hila::is_complex_or_arithmetic<P>::value &&
//                                hila::is_assignable<T &, hila::type_mul<T, P>>::value,
//                            int> = 0>
// inline T pow(const T &a, const P &p) {
//     return back<T>(::pow(to_array(a), p);
// }

template <typename T, typename P, typename R = hila::type_mul<T, P>,
          std::enable_if_t<hila::is_complex_or_arithmetic<P>::value, int> = 0>
inline R pow(const T &a, const P &p) {
    return back<R>(::pow(to_array(a), p));
}

/// @brief addition
/// @details Function hila::elem::add(a,b)  does element-by-element addition of same shape
/// Matrix/Vector/Array.  If one of the arguments is scalar or Complex, this is added to all
/// elements.  For matrices this is equivalent to (a.asArray() + b.asArray()).asMatrix(),
/// except that for hila::elem::add() the return value is the same as the input matrix type
/// if possible.


template <typename A, typename B, typename R = hila::type_plus<A, B>,
          typename V = decltype(to_array(std::declval<A>()) + to_array(std::declval<B>()))>
inline R add(const A &a, const B &b) {
    return back<R>(to_array(a) + to_array(b));
}

/// @brief subtract

template <typename A, typename B, typename R = hila::type_plus<A, B>,
          typename V = decltype(to_array(std::declval<A>()) + to_array(std::declval<B>()))>
inline R sub(const A &a, const B &b) {
    return back<R>(to_array(a) - to_array(b));
}

/// @brief multiply
/// Result is still "type_plus<A,B>", because that gives correct type for
/// element-by-element multiply

template <typename A, typename B, typename R = hila::type_plus<A, B>,
          typename V = decltype(to_array(std::declval<A>()) * to_array(std::declval<B>()))>
inline R mul(const A &a, const B &b) {
    return back<R>(to_array(a) * to_array(b));
}

/// @brief divide
/// Result is still "type_plus<A,B>", because that gives correct type for
/// element-by-element multiply

template <typename A, typename B, typename R = hila::type_plus<A, B>,
          typename V = decltype(to_array(std::declval<A>()) / to_array(std::declval<B>()))>
inline R div(const A &a, const B &b) {
    return back<R>(to_array(a) / to_array(b));
}


} // namespace elem
} // namespace hila


#endif // HILA_ELEMENT_BY_ELEMENT
