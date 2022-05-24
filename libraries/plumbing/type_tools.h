#ifndef TEMPLATE_TOOLS_H
#define TEMPLATE_TOOLS_H

#include <iostream>
#include <assert.h>
#include <sstream>

#include "plumbing/defs.h"

// // Useful c++14 template missing in c++11
// #if defined(PUHTI) && defined(HILAPP)
// namespace std {
// template <bool B, class T = void>
// using enable_if_t = typename std::enable_if<B, T>::type;
// }
// #endif

namespace hila {

/////////////////////////////////////////////////////////////////////////////////
/// Utility for obtaining the numeric base type of a class
/// Use as  "hila::number_type<T>"
/// which returns the underlying arithmetic type used in class T (e.g. float, double
/// etc.)
///
/// Alt form  "typename base_type_struct<T>"

template <typename T, typename Enable = void>
struct base_type_struct {
    /// The base type of the class
    using type = typename T::base_type;
};

template <typename T>
struct base_type_struct<T, typename std::enable_if_t<hila::is_arithmetic<T>::value>> {
    // In this case the base type is just T
    using type = T;
};

template <typename T>
using number_type = typename base_type_struct<T>::type;

//////////////////////////////////////////////////////////////////////////////
/// Get the inner type in nesting templates:
///    inner_type<A<B<C>>>
/// returns type B<C>
/// If type is arithmetic, returns it
//////////////////////////////////////////////////////////////////////////////
template <typename T, typename Enable = void>
struct inner_type_struct {
    using type = typename T::argument_type;
};

template <typename T>
struct inner_type_struct<T, typename std::enable_if_t<hila::is_arithmetic<T>::value>> {
    using type = T;
};

template <typename T>
using inner_type = typename inner_type_struct<T>::type;

//////////////////////////////////////////////////////////////////////////////
/// Define boolean contains_type<A,B>::value, which returns true if
/// type A "wraps" type B at some level
/// Example:
///    hila::contains_type<Vector<4,Complex<double>>>, Complex<double>>::value
/// is true
//////////////////////////////////////////////////////////////////////////////

template <typename A, typename B, typename Enable = void>
struct contains_type
    : std::integral_constant<bool, hila::contains_type<inner_type<A>, B>::value> {};

template <typename A>
struct contains_type<A, A> : std::integral_constant<bool, true> {};

template <typename A, typename B>
struct contains_type<A, B, typename std::enable_if_t<hila::is_arithmetic<A>::value>>
    : std::integral_constant<bool, std::is_same<A, B>::value> {};

//////////////////////////////////////////////////////////////////////////////
/// Helper operations to make generic templates for arithmetic operators
/// e.g. hila::type_plus<A,B> gives the type of the operator a + b, where a is of type A
/// and b type B.
//////////////////////////////////////////////////////////////////////////////
template <typename A, typename B>
using type_plus = decltype(std::declval<A>() + std::declval<B>());
template <typename A, typename B>
using type_minus = decltype(std::declval<A>() - std::declval<B>());
template <typename A, typename B>
using type_mul = decltype(std::declval<A>() * std::declval<B>());
template <typename A, typename B>
using type_div = decltype(std::declval<A>() / std::declval<B>());


//////////////////////////////////////////////////////////////////////////////
/// Access variables as if arrays of number_type numbers
//////////////////////////////////////////////////////////////////////////////

template <typename T>
inline hila::number_type<T> get_number_in_var(const T &var, int i) {
    return *(reinterpret_cast<const hila::number_type<T> *>(&var) + i);
}
template <typename T>
inline void set_number_in_var(T &var, int i, const hila::number_type<T> val ) {
    *(reinterpret_cast<const hila::number_type<T> *>(&var) + i) = val;
}




} // namespace hila

#endif