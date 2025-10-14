#ifndef HILA_CMPLX_UTILITIES_H_
#define HILA_CMPLX_UTILITIES_H_

// This header file is meant to be included from cmplx.h
/**
 * @file cmplx_utilities.h
<<<<<<< Updated upstream
 * @brief defines complex template utilities, in namespace hila::
=======
 * @brief defines complex template utilities, in namespace hila
>>>>>>> Stashed changes
 *
 */

namespace hila {
////////////////////////////////////////////////////////////////////////
// Section for adding complex utility templates (in namespace hila)
// hila::is_complex<T>::value
// hila::is_complex_or_arithmetic<T>::value
// hila::contains_complex<T>::value
// hila::number_type<T>   returns Complex or arithmetic type
// hila::ntype_op<A,B>        returns the conventionally upgraded complex or scalar number type
// hila::complex_x_scalar_type<A,B>  type of operation Complex<A> * scalar<B>
// hila::number_type_double<T>  returns double or Complex<double> (compare to number_type)

/**
 * @brief  hila::is_complex<T>::value is true if T is Complex<float> or Complex<double>
 */
template <typename T>
struct is_complex : std::integral_constant<bool, false> {};

template <typename T>
struct is_complex<Complex<T>> : std::integral_constant<bool, true> {};

template <typename T>
struct is_complex<Imaginary_t<T>> : std::integral_constant<bool, true> {};

/***
 * hila::is_complex_or_arithmetic<T>::value is true for complex or arithmetic types
 */
template <typename T>
struct is_complex_or_arithmetic
    : std::integral_constant<bool, hila::is_arithmetic<T>::value || hila::is_complex<T>::value> {};

/////////////////////////////////////////////////////////////////////////
/**
 * @brief Utility to check that the type contains complex numbers (only for hila types)
 * Use as contains_complex<T>::value
 */
template <typename T, typename Enable = void>
struct contains_complex : std::integral_constant<bool, false> {};

template <typename T>
struct contains_complex<T, typename std::enable_if_t<hila::is_field_class_type<T>::value>>
    : std::integral_constant<bool,
                             hila::contains_type<T, Complex<hila::arithmetic_type<T>>>::value> {};

/////////////////////////////////////////////////////////////////////////
template <typename T, typename Enable = void, typename N = hila::arithmetic_type<T>>
struct complex_or_arithmetic_type_struct {
    using type = N;
};

template <typename T>
struct complex_or_arithmetic_type_struct<
    T, typename std::enable_if_t<hila::contains_complex<T>::value>> {
    using type = Complex<hila::arithmetic_type<T>>;
};

/**
 * @brief Utility hila::number_type<T>  returns the number type from which T is constructed,
 * either arithmetic type or Complex
 */
template <typename T>
using number_type = typename complex_or_arithmetic_type_struct<T>::type;

/////////////////////////////////////////////////////////////////////////
// Utility hila::ntype_op<T1,T2>  returns real or complex type,
// "conventionally" defined double/float etc. upgrading.
// Note:  hila::type_plus<T1,T2> gives e.g. Complex<float>,double -> Complex<float>

template <typename T1, typename T2, typename Enable = void>
struct ntype_op_s {
    using type = hila::type_plus<hila::arithmetic_type<T1>, hila::arithmetic_type<T2>>;
};

// if one type is complex?
template <typename T1, typename T2>
struct ntype_op_s<T1, T2,
                  typename std::enable_if_t<(hila::contains_complex<T1>::value ||
                                             hila::contains_complex<T2>::value)>> {
    using type = Complex<hila::type_plus<hila::arithmetic_type<T1>, hila::arithmetic_type<T2>>>;
};

template <typename T1, typename T2>
using ntype_op = typename ntype_op_s<T1, T2>::type;

////////////////////////////////////////////////////////////////////////////
// utility complex_x_scalar, used in Complex<A> + scalar<B> ops
// - if result is convertible to Complex<A>, return Complex<A>
// - otherwise return Complex<type_plus<A,B>>
// This is done in order to keep vector types as vectors

template <typename A, typename B, typename Enable = void>
struct complex_x_scalar_s {
    using type = Complex<hila::type_plus<A, B>>;
};

template <typename A, typename B>
struct complex_x_scalar_s<
    A, B, typename std::enable_if_t<std::is_assignable<A &, hila::type_plus<A, B>>::value>> {
    using type = Complex<A>;
};

template <typename A, typename B>
using complex_x_scalar_type = typename complex_x_scalar_s<A, B>::type;

////////////////////////////////////////////////////////////////////////
// as_complex_array(T var)
// casts the var to a pointer to complex<scalar_type<T>>
// assuming var contains complex type.  This enables access
// of complex elements as
//  as_complex_array(var)[i]
// comment out as hilapp gets confused at the moment
////////////////////////////////////////////////////////////////////////

// template <typename T, std::enable_if_t<hila::contains_complex<T>::value, int> = 0>
// inline const Complex<hila::arithmetic_type<T>> * as_complex_array(const T &var) {
//     return (const Complex<hila::arithmetic_type<T>> *)(void *)&var;
// }

// template <typename T, std::enable_if_t<hila::contains_complex<T>::value, int> = 0>
// inline Complex<hila::arithmetic_type<T>> * as_complex_array(T &var) {
//     return (Complex<hila::arithmetic_type<T>> *)(void *)&var;
// }

} // namespace hila

#endif