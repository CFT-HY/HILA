#ifndef TEMPLATE_TOOLS_H
#define TEMPLATE_TOOLS_H

#include <iostream>
#include <array>
#include <vector>
#include <assert.h>
#include <sstream>

#include "plumbing/defs.h"


/// Utility for selecting the numeric base type of a class
template <class T, class Enable = void> struct base_type_struct {
    /// The base type of the class
    using type = typename T::base_type;
};

/// Utility for selecting the numeric base type of a class
template <typename T>
struct base_type_struct<T, typename std::enable_if_t<is_arithmetic<T>::value>> {
    /// In this case the base type is just T
    using type = T;
};

/// Utility for selecting the numeric base type of a class
/// Use as number_type<T>
template <typename T> using number_type = typename base_type_struct<T>::type;



// Useful c++14 template missing in Puhti compilation of hilapp
#if defined(PUHTI) && defined(HILAPP)
namespace std {
template <bool B, class T = void>
using enable_if_t = typename std::enable_if<B, T>::type;
}
#endif

// These are helpers, to make generic templates
// e.g. type_plus<A,B> gives the type of the operator a + b, where a is of type A and b
// B.
template <typename A, typename B>
using type_plus = decltype(std::declval<A>() + std::declval<B>());
template <typename A, typename B>
using type_minus = decltype(std::declval<A>() - std::declval<B>());
template <typename A, typename B>
using type_mul = decltype(std::declval<A>() * std::declval<B>());
template <typename A, typename B>
using type_div = decltype(std::declval<A>() / std::declval<B>());



#endif