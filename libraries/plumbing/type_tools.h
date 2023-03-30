#ifndef TEMPLATE_TOOLS_H
#define TEMPLATE_TOOLS_H

#include <iostream>
#include <assert.h>
#include <sstream>

#include "plumbing/defs.h"


namespace hila {

/////////////////////////////////////////////////////////////////////////////////
/// is class type defined so that can be used in Field<type> -variables
/// use the existence of type::base_type -element as a proxy for this
/// use: hila::is_field_class_type<T>::value
template <class, class = void>
struct is_field_class_type : std::false_type {};

template <class T>
struct is_field_class_type<T, std::void_t<typename T::base_type>> : std::true_type {};

/////////////////////////////////////////////////////////////////////////////////
/// General routine to determine if type can be used in Field<type> -variables
/// union of is_arithmetic and is_field_class_type
/// use: hila::is_field_type<T>::value
template <class, class = void>
struct is_field_type : std::false_type {};

template <class T>
struct is_field_type<T, typename std::enable_if_t<is_field_class_type<T>::value ||
                                                  is_arithmetic<T>::value>> : std::true_type {
};

/////////////////////////////////////////////////////////////////////////////////
/// Utility for obtaining the numeric base type of a class
/// Use as  "hila::scalar_type<T>"
/// which returns the underlying arithmetic type used in class T (e.g. float, double
/// etc.)  if it is field class type.  For other types returns the input type
/// 
///
/// Alt form  "typename base_type_struct<T>"

template <typename T, typename Enable = void>
struct base_type_struct {
    using type = T;
};

template <typename T>
struct base_type_struct<T, typename std::enable_if_t<is_field_class_type<T>::value>> {
    // In this case the base type is just T
    using type = typename T::base_type;
};

template <typename T>
using scalar_type = typename base_type_struct<T>::type;

//////////////////////////////////////////////////////////////////////////////
/// Get the inner type in nesting templates:
///    inner_type<A<B<C>>>
/// returns type B<C>
/// If type is not field class type, returns input
//////////////////////////////////////////////////////////////////////////////

template <typename T, typename Enable = void>
struct inner_type_struct {
    using type = T;
};

template <typename T>
struct inner_type_struct<T, typename std::enable_if_t<hila::is_field_class_type<T>::value>> {
    using type = typename T::argument_type;
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
struct contains_type : std::integral_constant<bool, std::is_same<A, B>::value> {};

template <typename A, typename B>
struct contains_type<A, B, typename std::enable_if_t<hila::is_field_class_type<A>::value && !std::is_same<A,B>::value>>
    : std::integral_constant<bool, hila::contains_type<hila::inner_type<A>, B>::value> {};

template <typename A>
struct contains_type<A, A> : std::integral_constant<bool, true> {};


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
/// Access variables as if arrays of scalar_type numbers
/// 
//////////////////////////////////////////////////////////////////////////////

template <typename T>
inline hila::scalar_type<T> get_number_in_var(const T &var, int i) {
    return *(reinterpret_cast<const hila::scalar_type<T> *>(&var) + i);
}
template <typename T>
inline void set_number_in_var(T &var, int i, const hila::scalar_type<T> val) {
    *(reinterpret_cast<hila::scalar_type<T> *>(&var) + i) = val;
}


} // namespace hila

#endif