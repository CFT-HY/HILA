#ifndef TEMPLATE_TOOLS_H
#define TEMPLATE_TOOLS_H

#include <iostream>
#include <assert.h>
#include <sstream>

#include "plumbing/defs.h"

namespace hila {

/////////////////////////////////////////////////////////////////////////////////
/// Utility for obtaining the numeric base type of a class
///   Use as "typename base_type_struct<T>" which returns the "arithmetic type"
/// which is used to construct T

template <typename T, typename Enable = void> struct base_type_struct {
    /// The base type of the class
    using type = typename T::base_type;
};

/// Utility for selecting the numeric base type of a class
template <typename T>
struct base_type_struct<T, typename ::std::enable_if_t<hila::is_arithmetic<T>::value>> {
    // In this case the base type is just T
    using type = T;
};

/// Main utility for obtaining the numeric base type of a hila class.
/// Use as hila::number_type<T>
template <typename T> using number_type = typename base_type_struct<T>::type;

/// struct to return the inner type, i.e.
///    typename inner_type_struct<A<B>>::type
/// returns type B.
//////////////////////////////////////////////////////////////////////////////
template <typename T, typename Enable = void> struct inner_type_struct {
    using type = typename T::argument_type;
};

template <typename T>
struct inner_type_struct<T, typename ::std::enable_if_t<hila::is_arithmetic<T>::value>> {
    using type = T;
};

//////////////////////////////////////////////////////////////////////////////
/// Main interface to inquire inner type in nesting templates:
///    inner_type<A<B<C>>>
/// teturns type B<C>
//////////////////////////////////////////////////////////////////////////////

template <typename T> using inner_type = typename inner_type_struct<T>::type;

//////////////////////////////////////////////////////////////////////////////
/// Define boolean contains_type<A,B>::value, which returns true if
/// type A "wraps" type B at some level
/// Example:
///    hila::contains_type<Vector<4,Complex<double>>>, Complex<double>>::value
/// is true
//////////////////////////////////////////////////////////////////////////////

template <typename A, typename B, typename Enable = void>
struct contains_type
    : ::std::integral_constant<bool, hila::contains_type<inner_type<A>, B>::value> {};

template <typename A>
struct contains_type<A, A> : ::std::integral_constant<bool, true> {};

template <typename A, typename B>
struct contains_type<A, B, typename ::std::enable_if_t<hila::is_arithmetic<A>::value>>
    : ::std::integral_constant<bool, ::std::is_same<A, B>::value> {};

// Useful c++14 template missing in Puhti compilation of hilapp
#if defined(PUHTI) && defined(HILAPP)
namespace std {
template <bool B, class T = void>
using enable_if_t = typename ::std::enable_if<B, T>::type;
}
#endif

//////////////////////////////////////////////////////////////////////////////
/// Helper operations to make generic templates for arithmetic operators
/// e.g. hila::type_plus<A,B> gives the type of the operator a + b, where a is of type A
/// and b type B.
//////////////////////////////////////////////////////////////////////////////
template <typename A, typename B>
using type_plus = decltype(::std::declval<A>() + ::std::declval<B>());
template <typename A, typename B>
using type_minus = decltype(::std::declval<A>() - ::std::declval<B>());
template <typename A, typename B>
using type_mul = decltype(::std::declval<A>() * ::std::declval<B>());
template <typename A, typename B>
using type_div = decltype(::std::declval<A>() / ::std::declval<B>());

} // namespace hila

#endif