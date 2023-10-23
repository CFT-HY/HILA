/**
 * @file cmplx.h
 * @brief Definition of Complex types
 * @details This file contains definitions and methods for Complex numbers and Imaginary type
 *
 */

#ifndef CMPLX_H_
#define CMPLX_H_

// let's not include the std::complex
// #include <complex>
// #include <cmath>

#include "plumbing/defs.h"

/// TEMPORARY location for vector intrinsic analogues -- result obvious

template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline T mul_add(T a, T b, T c) {
    return a * b + c;
}

template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline T mul_sub(T a, T b, T c) {
    return a * b - c;
}

template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline T nmul_add(T a, T b, T c) {
    return c - a * b;
}

/**
 * @brief Complex definition
 * @details Define complex type as a class. This allows Hilapp to replace the internal type with
 * a vector.
 *
 * __NOTE:__ T must be arithmetic and integrable. In the following documentation MyType refers to T,
 * as in an arithmetic and integrable type.
 * @param re Real part
 * @param im Imaginary part
 * @tparam T Arithmetic type
 */
template <typename T>
class Complex {

    static_assert((hila::is_arithmetic<T>::value && !std::is_integral<T>::value),
                  "Complex can be used only with floating point type");
    // This incantation is needed to make Field<Complex<>> vectorized

  public:
    using base_type = hila::arithmetic_type<T>;
    using argument_type = T;

    // and the content of the complex number
    T re, im;

    /**
     * @brief Construct a new Complex object
     * @details The following ways of constructing a Complex object are
     *
     * __Default constructor:__
     *
     * \code{.cpp}
     * Complex<MyType> C;
     * \endcode
     *
     * The default constructor initializes \p Complex#re and \p Complex#im to 0
     *
     * __Complex constructor:__
     *
     * Initialize both real and imaginary element
     *
     * \code{.cpp}
     * MyType a,b;
     * a = hila::random();
     * b = hila::random();
     * Complex<MyType> C(a,b); // C.re = a, C.im = b
     * \endcode
     * __Copy constructor:__
     *
     * Initialize form already existing Complex number
     *
     * \code {.cpp}
     * MyType a,b;
     * a = hila::random();
     * b = hila::random();
     * Complex<MyType> C(a,b);
     * Complex<MyType> B = C; // B.re = a, B.im = b
     * \endcode
     *
     * Equivalent initializing is `Complex<MyType> B(C)`
     *
     * __Real constructor:__
     *
     * Initialize only real element and sets imaginary to 0
     *
     * \code {.cpp}
     * MyType a = hila::random();
     * Complex<MyType> C(a); // C.re = a, C.im = 0
     * \endcode
     *
     * Not equivalent to `Complex<MyType> C = a`
     *
     * __Zero constructor:__
     *
     * Initialize to zero with nullpointer trick
     *
     * \code {.cpp}
     * Complex<MyType> C = 0; // C.re = 0, C.im = 0
     * \endcode
     *
     * Equivalent to `Complex<MyType> C` and `Complex<MyType> C(0)`
     *
     */
    Complex<T>() = default;
    ~Complex<T>() = default;
    Complex<T>(const Complex<T> &a) = default;

    // constructor from single complex --IS THIS NEEDED?
    // template <typename A,
    //           std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0 >
    // constexpr Complex<T>(const Complex<A> a) : re(static_cast<T>(a.re)),
    // im(static_cast<T>(a.im)) {}

    // constructor from single scalar value
    // Remember to mark this explicit, we do not want this to be invoked
    // in automatic conversions (there should be methods)

//#pragma hila loop_function
    template <typename S, std::enable_if_t<hila::is_arithmetic<S>::value, int> = 0>
    explicit constexpr Complex<T>(const S val) : re(val), im(0) {}

    // allow construction from 0 (nullptr)
    constexpr Complex<T>(const std::nullptr_t n) : re(0), im(0) {}

    // constructor c(a,b)
    template <typename A, typename B, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0,
              std::enable_if_t<hila::is_arithmetic<B>::value, int> = 0>
//#pragma hila novector loop_function
    explicit constexpr Complex<T>(const A &a, const B &b) : re(a), im(b) {
    }

    // make also std accessors real() and imag() - don't return reference, because
    // v.real() would not then work inside loops!


    /**
     * @brief Real part of Complex number
     * @details Overloads exist for pass by value and pass by reference for when Complex is defined
     * as non-const or const
     *
     * @return T Real part
     */
    inline T real() const {
        return re;
    }
    /**
     * @brief Imaginary part of Complex number
     * @details Overloads exist for pass by value and pass by reference for when Complex is defined
     * as non-const or const
     *
     * @return T Imaginary part
     */
    inline T imag() const {
        return im;
    }

    inline T &real() const_function {
        return re;
    }

    inline T &imag() const_function {
        return im;
    }
    // automatic casting from Complex<T> -> Complex<A>
    // TODO: ensure this works if A is vector type!
    template <typename A>
#pragma hila loop_function // TODO
    operator Complex<A>() const {
        return Complex<A>(re, im);
    }

    /**
     * @brief Assignment operator
     * @details The following ways of assigning a Complex object are
     *
     * __Assignment for Complex:__
     *
     * \code {.cpp}
     * Complex<MyType> C(hila::random(),hila::random());
     * Complex<MyType> B;
     * B = C; // B.re = C.re, B.im = C.im
     * \endcode
     *
     * __Assignment from real:__
     *
     * Assigns real part and imaginary part to zero
     *
     * \code {.cpp}
     * Complex<MyType> B;
     * MyType a = hila::random;
     * B = a; // B.re = a, B.im = 0
     * \endcode
     *
     * @param s
     * @return Complex<T>&
     */
    inline Complex<T> &operator=(const Complex<T> &s) = default;

    // Assignment from Complex<A>
    template <typename A>
    inline Complex<T> &operator=(const Complex<A> &s) {
        re = s.re;
        im = s.im;
        return *this;
    }

// #pragma hila loop_function
    template <typename S, std::enable_if_t<hila::is_arithmetic<S>::value, int> = 0>
    inline Complex<T> &operator=(S s) {
        re = s;
        im = 0;
        return *this;
    }

    /**
     * @brief Compute square norm of Complex number
     * @details
     * \f{align}{ |z|^2 = \Re(z)^2 + \Im(z)^2\f}
     *
     * @return T Squared norm
     */
    inline T squarenorm() const {
        return re * re + im * im;
    }

    /**
     * @brief Compute absolute value of Complex number
     * @details Essentially sqrt(squarenorm(z)):
     * \f{align}{ |z| = \sqrt{\Re(z)^2 + \Im(z)^2}\f}
     *
     * @return T
     */
    inline T abs() const {
        return sqrt(squarenorm());
    }

    /**
     * @brief Compute argument of Complex number
     * @details
     * \f{align}{ \arg(z) = \arctan2(\Re(z),\Im(z))\f}
     *
     * @return T
     */
    inline T arg() const {
        return atan2(im, re);
    }

    /**
     * @brief Compute conjugate of Complex number
     * @details
     * \f{align}{ z &= x + i\cdot y \\
     * \Rightarrow z^* &= x - i\cdot y\f}
     *
     * @return Complex<T>
     */
    inline Complex<T> conj() const {
        return Complex<T>(re, -im);
    }

    /**
     * @brief Compute dagger of Complex number
     * @details Alias to Complex::conj
     *
     * \f{align}{ z^* = z^\dagger \f}
     * @return Complex<T>
     */
    inline Complex<T> dagger() const {
        return Complex<T>(re, -im);
    }

    /**
     * @brief Stores and returns Complex number given in polar coordinates
     *
     * \f{align}{ z = r\cdot e^{i\theta} \f}
     * \code {.cpp}
     * Complex<double> c;
     * double r = 1;
     * double theta = 3.14159265358979/2; // pi/2
     * c.polar(r,theta); // c.re = 0, c.im = 1
     * \endcode
     *
     *
     * @param r Radius of Complex number
     * @param theta Angle of complex number in radians
     * @return Complex<T> Complex number
     */
    inline Complex<T> polar(const T r, const T theta) out_only {
        re = r * cos(theta);
        im = r * sin(theta);
        return *this;
    }

    /**
     * @brief Assign random values to Complex real and imaginary part
     * @details Uses hila::random for both real and imaginary part
     *
     * @return Complex<T>&
     */
    inline Complex<T> &random() out_only {
        re = hila::random();
        im = hila::random();
        return *this;
    }

    /**
     * @brief Produces complex gaussian random values
     * @details Uses hila::gaussrand2 for both real and imaignary part
     * Assigns same random value for both real and imaginary part
     * @param width gaussian_random
     * @return Complex<T>&
     */
    inline Complex<T> &gaussian_random(double width = 1.0) out_only {
        double d;
        re = hila::gaussrand2(d) * width;
        im = d * width;
        return *this;
    }

    /**
     * @brief Unary + operator
     * @details Leaves Complex number unchanged
     * @return Complex<T>
     */
    inline Complex<T> operator+() const {
        return *this;
    }

    /**
     * @brief Unary - operator
     * @details Negates Complex number
     *
     * @return Complex<T>
     */
    inline Complex<T> operator-() const {
        return Complex<T>(-re, -im);
    }

    // mark += and *= as loop_functions, because these can be used
    // in reductions implicitly

    /**
     * @brief += addition assignment operator
     * @details Addition assignment for Complex numbers can be performed in the following ways
     *
     * __Complex addition assignment:__
     *
     * \code{.cpp}
     * Complex<double> z(0,0);
     * Complex<double> w(1,1);
     * z += w; // z.re = 1, z.im = 1
     * \endcode
     *
     * __Real addition assignment:__
     *
     * Add assign only to real part of Complex number
     *
     * \code {.cpp}
     * Complex<double> z(0,0);
     * z += 1; // z.re = 1, z.im = 0
     * \endcode
     *
     *
     */
//#pragma hila loop_function
    template <typename A>
    inline Complex<T> &operator+=(const Complex<A> &lhs) {
        re += lhs.re;
        im += lhs.im;
        return *this;
    }

    template <typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
    inline Complex<T> &operator+=(const A &a) {
        re += a;
        return *this;
    }

    /**
     * @brief -= subtraction assignment operator
     * @details Subtraction assignment for Complex numbers can be performed in the following ways
     *
     * __Complex subtract assign:__
     *
     * \code{.cpp}
     * Complex<double> z(0,0);
     * Complex<double> w(1,1);
     * z -= w; // z.re = -1, z.im = -1
     * \endcode
     *
     * __Real subtract assign:__
     *
     * Subtract assign only to real part of Complex number
     *
     * \code {.cpp}
     * Complex<double> z(0,0);
     * z -= 1; // z.re = -1, z.im = 0
     * \endcode
     *
     *
     */
    template <typename A>
    inline Complex<T> &operator-=(const Complex<A> &lhs) {
        re -= lhs.re;
        im -= lhs.im;
        return *this;
    }

    template <typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
    inline Complex<T> &operator-=(const A &a) {
        re -= a;
        return *this;
    }

    // inline Complex<T> & operator*= (const Complex<T> & lhs) {
    //   T r = re * lhs.re - im * lhs.im;
    //   im  = im * lhs.re + re * lhs.im;
    //   re = r;
    //   return *this;
    // }

    /**
     * @brief *= multiply assignment operator
     * @details Multiply assignment for Complex numbers can be performed in the following ways
     *
     * __Complex multiply assign:__
     *
     * Standard Complex number multiplication
     *
     * \f{align}{z &= x + iy, w = x' + iy' \\
     * z w &= (x + iy)(x' + iy') = (xx'-yy') + i(xy' + yx')\f}
     * \code{.cpp}
     * Complex<double> z,w;
     * //
     * // z,w get values
     * //
     * z*=w; // z = zw as defined above
     * \endcode
     * __Real multiply assign:__
     *
     * Multiply assign by real number to both components of Complex number
     *
     * \code {.cpp}
     * Complex<double> z(1,1);
     * z *= 2; // z.re = 2, z.im = 2
     * \endcode
     *
     *
     */
//#pragma hila loop_function
    template <typename A>
    inline Complex<T> &operator*=(const Complex<A> &lhs) {
        T r = mul_sub(re, lhs.re, im * lhs.im); // a*b-c
        im = mul_add(im, lhs.re, re * lhs.im);  // a*b+c
        re = r;
        return *this;
    }

    template <typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
    inline Complex<T> &operator*=(const A a) {
        re *= a;
        im *= a;
        return *this;
    }

    // a/b = a b*/|b|^2 = (a.re*b.re + a.im*b.im + i(a.im*b.re - a.re*b.im))/|b|^2
    // inline Complex<T> & operator/= (const Complex<T> & lhs) {
    //   T n = lhs.squarenorm();
    //   T r = (re * lhs.re + im * lhs.im)/n;
    //   im  = (im * lhs.re - re * lhs.im)/n;
    //   re = r;
    //   return *this;
    // }

    /**
     * @brief /= divide assignment operator
     * @details Divide assignment for Complex numbers can be performed in the following ways
     *
     * __Complex divide assign:__
     *
     * Standard Complex number division
     *
     * \f{align}{z &= x + iy, w = x' + iy' \\
     * \frac{z}{w} &= \frac{x + iy}{x' + iy'} = \frac{(xx'+ yy') + i( yx' - xy')}{|w|^2}\f}
     * \code{.cpp}
     * Complex<double> z,w;
     * //
     * // z,w get values
     * //
     * z/=w; // z = z/w as defined above
     * \endcode
     * __Real divide assign:__
     *
     * Divide assign by real number to both components of Complex number
     *
     * \code {.cpp}
     * Complex<double> z(2,2);
     * z /= 2; // z.re = 1, z.im = 1
     * \endcode
     *
     *
     */
    template <typename A>
    inline Complex<T> &operator/=(const Complex<A> &lhs) {
        T n = lhs.squarenorm();
        T r = mul_add(re, lhs.re, im * lhs.im) / n; // a*b+c
        im = mul_sub(im, lhs.re, re * lhs.im) / n;  // a*b-c
        re = r;
        return *this;
    }

    template <typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
    inline Complex<T> &operator/=(const A &a) {
        re /= a;
        im /= a;
        return *this;
    }

    // define also increment and decrement operators (though not so intuitive for complex)

    /**
     * @brief ++ increment operator
     *
     * Increments real part of Complex number
     *
     * \code {.cpp}
     * Complex<double> z(1,1);
     * z++; // z.re = 2, z.im = 1
     * \endcode
     *
     *
     * @return Complex<T>&
     */
    inline Complex<T> &operator++() {
        this->re++;
        return *this;
    }
    /**
     * @brief -- decrement operator
     *
     * Decrement real part of Complex number
     *
     * \code {.cpp}
     * Complex<double> z(1,1);
     * z--; // z.re = 0, z.im = 1
     * \endcode
     *
     *
     * @return Complex<T>&
     */
    inline Complex<T> &operator--() {
        this->re--;
        return *this;
    }

    inline Complex<T> operator++(int) {
        Complex<T> a = *this;
        this->re++;
        return a;
    }

    inline Complex<T> operator--(int) {
        Complex<T> a = *this;
        this->re--;
        return a;
    }

    std::string str(int prec = 8, char separator = ' ') const {
        std::stringstream ss;
        ss.precision(prec);
        ss << re << separator << im;
        return ss.str();
    }


    // Convenience method a.conj_mul(b) == a^* b

    /**
     * @brief Conjugate multiply method
     * @details Conjugate a (*this) and multiply with give argument b: `a.conj_mul(b)` \f$ \equiv
     * a^* \cdot b\f$
     *
     * @tparam A Type for Complex number b
     * @param b number to multiply with
     * @return Complex<T>
     */
    template <typename A>
    inline Complex<T> conj_mul(const Complex<A> &b) const {
        return Complex<T>(re * b.re + im * b.im, re * b.im - im * b.re);
    }

    // Convenience method a.mul_conj(b) == a * b^*

    /**
     * @brief Multiply conjugate method
     * @details Multiply a (*this) with conjugate of given argument b: `a.mul_conj(b)` \f$ \equiv a
     * \cdot b^*\f$
     *
     * @tparam A Type for Complex number b
     * @param b number to conjugate and multiply withv
     * @return Complex<T>
     */
    template <typename A>
    inline Complex<T> mul_conj(const Complex<A> &b) const {
        return Complex<T>(re * b.re + im * b.im, im * b.re - re * b.im);
    }

    // cast to another number type (IS THIS NEEDED FOR COMPLEX?)
    template <typename Ntype>
    Complex<Ntype> cast_to() const {
        Complex<Ntype> res;
        res.re = re;
        res.im = im;
        return res;
    }
};


namespace hila {
////////////////////////////////////////////////////////////////////////
/// Section for adding complex utility templates (in namespace hila)
/// hila::is_complex_or_arithmetic<T>::value
/// hila::contains_complex<T>::value
/// hila::number_type<T>   returns Complex or arithmetic type
/// hila::ntype_op<A,B>        returns the conventionally upgraded complex or scalar number type
/// hila::complex_x_scalar_type<A,B>  type of operation Complex<A> * scalar<B>


////////////////////////////////////////////////////////////////////////
/// Define hila::is_complex<T>::value -template, using specialization
template <typename T>
struct is_complex : std::integral_constant<bool, false> {};

template <typename T>
struct is_complex<Complex<T>> : std::integral_constant<bool, true> {};

/// hila::is_complex_or_arithmetic<T>::value
template <typename T>
struct is_complex_or_arithmetic
    : std::integral_constant<bool, hila::is_arithmetic<T>::value || hila::is_complex<T>::value> {};

/////////////////////////////////////////////////////////////////////////
/// Utility to check that the type contains complex numbers
/// Use as contains_complex<T>::value

template <typename T, typename Enable = void>
struct contains_complex : std::integral_constant<bool, false> {};

template <typename T>
struct contains_complex<T, typename std::enable_if_t<hila::is_field_class_type<T>::value>>
    : std::integral_constant<bool, hila::contains_type<T, Complex<hila::arithmetic_type<T>>>::value> {};

/////////////////////////////////////////////////////////////////////////
/// Utility hila::number_type<T>  returns complex or arithmetic type
/// depending on the type

template <typename T, typename Enable = void, typename N = hila::arithmetic_type<T>>
struct complex_or_arithmetic_type_struct {
    using type = N;
};

template <typename T>
struct complex_or_arithmetic_type_struct<
    T, typename std::enable_if_t<hila::contains_complex<T>::value>> {
    using type = Complex<hila::arithmetic_type<T>>;
};

template <typename T>
using number_type = typename complex_or_arithmetic_type_struct<T>::type;

/////////////////////////////////////////////////////////////////////////
/// Utility hila::ntype_op<T1,T2>  returns real or complex type,
/// "conventionally" defined double/float etc. upgrading.
/// Note:  hila::type_plus<T1,T2> gives e.g. Complex<float>,double -> Complex<float>

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
/// utility complex_x_scalar, used in Complex<A> + scalar<B> ops
/// - if result is convertible to Complex<A>, return Complex<A>
/// - otherwise return Complex<type_plus<A,B>>
/// This is done in order to keep vector types as vectors

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
/// as_complex_array(T var)
/// casts the var to a pointer to complex<scalar_type<T>>
/// assuming var contains complex type.  This enables access
/// of complex elements as
///  as_complex_array(var)[i]
/// comment out as hilapp gets confused at the moment
////////////////////////////////////////////////////////////////////////

// #pragma hila loop_function
// template <typename T, std::enable_if_t<hila::contains_complex<T>::value, int> = 0>
// inline const Complex<hila::arithmetic_type<T>> * as_complex_array(const T &var) {
//     return (const Complex<hila::arithmetic_type<T>> *)(void *)&var;
// }

// #pragma hila loop_function
// template <typename T, std::enable_if_t<hila::contains_complex<T>::value, int> = 0>
// inline Complex<hila::arithmetic_type<T>> * as_complex_array(T &var) {
//     return (Complex<hila::arithmetic_type<T>> *)(void *)&var;
// }

////////////////////////////////////////////////////////////////////////
/// get_complex_element(var,i)
///    returns the i:th complex number embedded in variable v
/// set_complex_element(var,i,value)
///    sets the i:th element in var
////////////////////////////////////////////////////////////////////////

template <typename T, std::enable_if_t<hila::contains_complex<T>::value, int> = 0>
inline Complex<hila::arithmetic_type<T>> get_complex_in_var(const T &var, int i) {
    return *(reinterpret_cast<const Complex<hila::arithmetic_type<T>> *>(&var) + i);
}

template <typename T, std::enable_if_t<hila::contains_complex<T>::value, int> = 0>
inline void set_complex_in_var(T &var, int i, const Complex<hila::arithmetic_type<T>> &val) {
    *(reinterpret_cast<Complex<hila::arithmetic_type<T>> *>(&var) + i) = val;
}

} // namespace hila

// generic Complex constructor - type from arguments
template <typename T, std::enable_if_t<hila::is_floating_point<T>::value, int> = 0>
Complex<T> complex(const T re, const T im) {
    return Complex<T>(re, im);
}


// functions real(), imag()

template <typename T>
inline T real(const Complex<T> &a) {
    return a.re;
}


template <typename T>
inline T imag(const Complex<T> &a) {
    return a.im;
}

template <typename T>
Complex<T> polar(T r, T arg) {
    Complex<T> res(r * cos(arg), r * sin(arg));
    return res;
}


// template <typename T>
// inline Complex<T> operator+(const Complex<T> & a, const Complex<T> & b) {
//   return Complex<T>(a.re + b.re, a.im + b.im);
// }

//#pragma hila loop_function
template <typename T1, typename T2, typename Tr = hila::type_plus<T1, T2>>
inline Complex<Tr> operator+(const Complex<T1> &a, const Complex<T2> &b) {
    return Complex<Tr>(a.re + b.re, a.im + b.im);
}

// TODO: for avx vector too -- #define new template macro
template <typename T, typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline auto operator+(const Complex<T> &c, const A &a) {
    return hila::complex_x_scalar_type<T, A>(c.re + a, c.im);
}

template <typename T, typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline auto operator+(const A &a, const Complex<T> &c) {
    return hila::complex_x_scalar_type<T, A>(c.re + a, c.im);
}

// -
// template <typename T>
// inline Complex<T> operator-(const Complex<T> & a, const Complex<T> & b) {
//   return Complex<T>(a.re - b.re, a.im - b.im);
// }
template <typename T1, typename T2, typename Tr = hila::type_plus<T1, T2>>
inline Complex<Tr> operator-(const Complex<T1> &a, const Complex<T2> &b) {
    return Complex<Tr>(a.re - b.re, a.im - b.im);
}

// complex - scalar
template <typename T, typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline auto operator-(const Complex<T> &c, const A &a) {
    return hila::complex_x_scalar_type<T, A>(c.re - a, c.im);
}

// scalar - complex
template <typename T, typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline auto operator-(const A &a, const Complex<T> &c) {
    return hila::complex_x_scalar_type<T, A>(a - c.re, -c.im);
}

///
/// complex * complex
template <typename T1, typename T2, typename Tr = hila::type_mul<T1, T2>>
inline Complex<Tr> operator*(const Complex<T1> &a, const Complex<T2> &b) {
    return Complex<Tr>(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}

/// complex * scalar
template <typename T, typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline auto operator*(const Complex<T> &c, const A &a) {
    return hila::complex_x_scalar_type<T, A>(c.re * a, c.im * a);
}

template <typename T, typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline auto operator*(const A &a, const Complex<T> &c) {
    return hila::complex_x_scalar_type<T, A>(a * c.re, a * c.im);
}

// /   a/b = ab*/|b|^2
// template <typename T>
// inline Complex<T> operator/(const Complex<T> & a, const Complex<T> & b) {
//   T n = b.squarenorm();
//   return Complex<T>( (a.re*b.re + a.im*b.im)/n, (a.im*b.re - a.re*b.im)/n );
// }
template <typename T1, typename T2, typename Tr = hila::type_mul<T1, T2>>
inline Complex<Tr> operator/(const Complex<T1> &a, const Complex<T2> &b) {
    T2 n = 1.0 / b.squarenorm();
    return Complex<Tr>((a.re * b.re + a.im * b.im) * n, (a.im * b.re - a.re * b.im) * n);
}

template <typename T, typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline auto operator/(const Complex<T> &c, const A &a) {
    return hila::complex_x_scalar_type<T, A>(c.re / a, c.im / a);
}

// a/c = ac*/|c|^2
template <typename T, typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline auto operator/(const A &a, const Complex<T> &c) {
    T n = c.squarenorm();
    return hila::complex_x_scalar_type<T, A>((a * c.re) / n, -(a * c.im) / n);
}

// write also multiply-add directly with complex numbers
template <typename T>
inline Complex<T> mul_add(const Complex<T> &a, const Complex<T> &b, const Complex<T> &c) {
    // a*b + c
    Complex<T> r;
    T t1 = mul_add(a.re, b.re, c.re);
    T t2 = mul_add(a.re, b.im, c.im);
    r.re = mul_add(a.im, b.im, t1); // -a.im*b.im + a.re*b.re + c.re
    r.im = mul_add(a.im, b.re, t2); // a.im*b.re + a.re*b.im + c.im
    return r;
}

/**
 * @brief Compare equality of two complex numbers
 *
 * Two numbers are equal, if the arithmetic values are equal: thus,
 * complex and real comparison  (a + i b) == a
 * is true if b == 0.
 *
 * @tparam A
 * @tparam B
 * @param a
 * @param b
 * @return true if values compare to equal
 */

template <typename A, typename B>
inline bool operator==(const Complex<A> &a, const Complex<B> &b) {
    return (a.re == b.re && a.im == b.im);
}

template <typename A, typename B, std::enable_if_t<hila::is_arithmetic<B>::value, int> = 0>
inline bool operator==(const Complex<A> &a, const B b) {
    return (a.re == b && a.im == 0);
}

template <typename A, typename B, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline bool operator==(const A a, const Complex<B> &b) {
    return b == a;
}

/**
 * @brief Compare non-equality of two complex numbers
 *
 * Negation of operator==()
 *
 * @tparam A
 * @tparam B
 * @param a
 * @param b
 * @return true if values are not arithmetically equal
 */

template <typename A, typename B>
inline bool operator!=(const Complex<A> &a, const Complex<B> &b) {
    return (a.re != b.re || a.im != b.im);
}

template <typename A, typename B, std::enable_if_t<hila::is_arithmetic<B>::value, int> = 0>
inline bool operator!=(const Complex<A> &a, const B b) {
    return (a.re != b || a.im != 0);
}

template <typename A, typename B, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline bool operator!=(const A a, const Complex<B> &b) {
    return b != a;
}


// Cast operators to different number or Complex type
// cast_to<double>(a);
// cast_to<float>(b);
// Cast from number->number, number->Complex, Complex->Complex OK,
//     Complex->number not.

template <typename Ntype, typename T, std::enable_if_t<hila::is_arithmetic<Ntype>::value, int> = 0>
Complex<Ntype> cast_to(const Complex<T> &m) {
    Complex<Ntype> res;
    res = m;
    return res;
}

//////////////////////////////////////////////////////////////////////////////////
// Some operations in function form.  Useful in templates when the arg type is not known

/// abs
template <typename T>
inline T abs(const Complex<T> &a) {
    return a.abs();
}

/// arg
template <typename T>
inline T arg(const Complex<T> &a) {
    return a.arg();
}

/// Conjugate
template <typename T>
inline Complex<T> conj(const Complex<T> &val) {
    return val.conj();
}

/// Alias dagger to conjugate
template <typename T>
inline Complex<T> dagger(const Complex<T> &val) {
    return val.conj();
}

/// norm_squared
template <typename T>
inline auto squarenorm(const Complex<T> &val) {
    return val.squarenorm();
}


//////////////////////////////////////////////////////////////////////////////////
/// Print a complex value as (re,im)
//////////////////////////////////////////////////////////////////////////////////

template <typename T>
std::ostream &operator<<(std::ostream &strm, const Complex<T> &A) {
    return strm << A.real() << ' ' << A.imag();
}

/// //////////////////////////////////////////////////////////////////////////////
/// Function hila::to_string

namespace hila {
template <typename T>
std::string to_string(const Complex<T> &A, int prec = 8, char separator = ' ') {
    return A.str(prec, separator);
}

template <typename T>
std::string prettyprint(const Complex<T> &A, int prec = 8) {
    std::stringstream ss;
    ss.precision(prec);
    ss << "( " << A.real() << ", " << A.imag() << " )";
    return ss.str();
}

} // namespace hila


/////////////////////////////////////////////////////////////////////////////
/// @brief Imaginary type, used to represent purely imaginary numbers
///
/// Useful for reducing multiply operations in im * complex or im * real -ops
/// Derived from Complex class, so generic complex ops should remain valid
/// Defines only operators * and /, others go via Complex class
///
/// Note: Imaginary_t should be used in Field variables
///
/// @tparam T  type of imaginary (float/double)
/////////////////////////////////////////////////////////////////////////////
template <typename T>
class Imaginary_t : public Complex<T> {
  public:
    constexpr Imaginary_t() = default;
    ~Imaginary_t() = default;
    constexpr Imaginary_t(const Imaginary_t &i) = default;

    // construct from scalar
#pragma hila loop_function
    template <typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
    explicit constexpr Imaginary_t(const A v) : Complex<T>(0, v) {}

    constexpr Imaginary_t operator-() const {
        return Imaginary_t<T>(-this->im);
    }

    constexpr Imaginary_t operator+() const {
        return *this;
    }
};

///////////////////////////////////////////////////////////////////////////////////////
/// @brief Imaginary unit I - global variable
///
/// Don't use #define'd I : this will conflict with some headers in rocm
///
/// For some reason it is sufficient to use only __device__
#if defined(CUDA) || defined(HIP)
__device__
#endif
    constexpr Imaginary_t<double>
        I(1.0);

// constexpr Complex<double> I(0,1);
// #define I Imaginary_t<double>(1.0)

///////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Imaginary * object containing complex
 * @param i  Imaginary value
 * @param c  multiplied value which contains complex numbers
 * @return   i * c, of the same type as c
 */
template <typename A, typename T, std::enable_if_t<hila::contains_complex<T>::value, int> = 0>
inline auto operator*(const Imaginary_t<A> &i, const T &c) {
    Complex<hila::arithmetic_type<T>> ca, cb;
    T res;
    constexpr int n_cmplx = sizeof(T) / sizeof(Complex<hila::arithmetic_type<T>>);
    for (int k = 0; k < n_cmplx; k++) {
        ca = hila::get_complex_in_var(c, k);
        cb.re = -ca.im * i.imag();
        cb.im = ca.re * i.imag();
        hila::set_complex_in_var(res, k, cb);
    }
    return res;
}

/// @brief object containing complex * imaginary
template <typename A, typename T, std::enable_if_t<hila::contains_complex<T>::value, int> = 0>
inline auto operator*(const T &c, const Imaginary_t<A> &i) {
    return i * c;
}


/// @brief Imag * scalar, returns imag
/// note: using std::is_arithmetic, not done for vector types
template <typename A, typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
inline Imaginary_t<A> operator*(Imaginary_t<A> i, const T &c) {
    i.imag() *= c;
    return i;
}

/// @brief scalar * imag, returns imag
template <typename A, typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
inline Imaginary_t<A> operator*(const T &c, Imaginary_t<A> i) {
    i.imag() *= c;
    return i;
}

/**
 * @brief  Imaginary * Imaginary, returning real value
 *
 * @return  -(a.im * b.im) real value, of appropriate type
 */
template <typename A, typename B>
inline auto operator*(const Imaginary_t<A> &a, const Imaginary_t<B> &b) {
    return -(a.imag() * b.imag());
}

////////////////////////////


/// @brief Imaginary / real value, returning imaginary
/// Note: for vectorized types the generic version is used
template <typename T, typename A, std::enable_if_t<std::is_arithmetic<A>::value, int> = 0>
inline auto operator/(Imaginary_t<T> i, const A &a) {
    i.imag() /= a;
    return i;
}

/// @brief Imaginary / Imaginary, returning real
template <typename A, typename B>
inline auto operator/(const Imaginary_t<A> &a, const Imaginary_t<B> &b) {
    return a.imag() / b.imag();
}


///////////////////////////////////////////////////////////////////////////////
/// Set of complex functions
///////////////////////////////////////////////////////////////////////////////


/// exp(z)
template <typename T>
inline Complex<T> exp(const Complex<T> z) {
    return exp(z.re) * Complex<T>(cos(z.im), sin(z.im));
}

/// exp(i x)
template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline Complex<T> expi(T a) {
    return Complex<T>(cos(a), sin(a));
}

/// @brief exp(imaginary)
template <typename T>
inline Complex<T> exp(const Imaginary_t<T> im) {
    return expi(im.imag());
}

/// log(z)
template <typename T>
inline Complex<T> log(Complex<T> z) {
    return Complex<T>(static_cast<T>(0.5) * log(z.squarenorm()), z.arg());
}

/// sqrt(z) branch cut at -x axis
template <typename T>
inline Complex<T> sqrt(Complex<T> z) {
    T r = z.squarenorm();
    T a = z.arg();
    return pow(r, 0.25) * expi(0.5 * a);
}

/// cbrt(z)
template <typename T>
inline Complex<T> cbrt(Complex<T> z) {
    T r = z.squarenorm();
    T a = z.arg();
    return pow(r, (1.0 / 6.0)) * expi((1.0 / 3.0) * a);
}

/// pow(z.p) = z^p = exp(p*log(z))
template <typename A, typename B>
inline auto pow(Complex<A> z, Complex<B> p) {
    return exp(p * log(z));
}

/// pow(z.p) with scalar power
template <typename T, typename S, std::enable_if_t<hila::is_arithmetic<S>::value, int> = 0>
inline Complex<T> pow(Complex<T> z, S p) {
    return exp(p * log(z));
}

/// pow(z.p) with scalar base
template <typename T, typename S, std::enable_if_t<hila::is_arithmetic<S>::value, int> = 0>
inline Complex<T> pow(S z, Complex<T> p) {
    return exp(p * log(z));
}

/// sin(z)
/// = sin(re + i im) = sin(re)cos(i im) + cos(re)sin(i im)
/// = sin(re) cosh(im) + i cos(re) sinh(im)
template <typename T>
inline Complex<T> sin(Complex<T> z) {
    return Complex<T>(sin(z.re) * cosh(z.im), cos(z.re) * sinh(z.im));
}

/// cos(z)
/// = cos(re)cos(i im) - sin(re)sin(i im) = cos(re)cosh(im) - i sin(re)sinh(im)
template <typename T>
inline Complex<T> cos(Complex<T> z) {
    return Complex<T>(cos(z.re) * cosh(z.im), -sin(z.re) * sinh(z.im));
}

/// tan(z) - rely on optimizer to simplify
template <typename T>
inline Complex<T> tan(Complex<T> z) {
    return sin(z) / cos(z);
}

/// sinh(z) = sinh(re)cosh(i im) + cosh(re)sinh(i im)
/// = sinh(re)cos(im) + i cosh(re)sin(im)
template <typename T>
inline Complex<T> sinh(Complex<T> z) {
    return Complex<T>(sinh(z.re) * cos(z.im), cosh(z.re) * sin(z.im));
}

/// cosh(z) = cosh(re)cosh(i im) - sinh(re)sinh(i im)
/// = cosh(re)cos(im) - i sinh(re)sin(im)
template <typename T>
inline Complex<T> cosh(Complex<T> z) {
    return Complex<T>(cosh(z.re) * cos(z.im), sinh(z.re) * sin(z.im));
}

/// tanh(z)
template <typename T>
inline Complex<T> tanh(Complex<T> z) {
    return sinh(z) / cosh(z);
}

/// arctan(z)
template <typename T>
inline Complex<T> atan(Complex<T> z) {
    return -0.5 * I * (log((I - z) / (I + z)));
}

/// arcsin(z)
template <typename T>
inline Complex<T> asin(Complex<T> z) {
    return -I * (log(I * z + sqrt(1 - z * z)));
}

/// arccos(z)
template <typename T>
inline Complex<T> acos(Complex<T> z) {
    return -I * (log(z + I * (sqrt(1 - z * z))));
}

/// artanh(z)
template <typename T>
inline Complex<T> atanh(Complex<T> z) {
    return 0.5 * log((1 + z) / (1 - z));
}

/// arsinh(z)
template <typename T>
inline Complex<T> asinh(Complex<T> z) {
    return log(z + sqrt(1 + z * z));
}

/// arcosh(z)
template <typename T>
inline Complex<T> acosh(Complex<T> z) {
    return log(z + sqrt(z * z - 1));
}

//////////////////////////////////////////////////////////////////////////////////
/// Operators to implement imaginary unit 1_i, enablig expressions  3 + 2_i  etc.
/// This is defined as an user-defined literal, which requires an underscore.
////////////////////////////////////////////////////////////////////////////////

constexpr Imaginary_t<double> operator""_i(long double a) {
    return Imaginary_t<double>{a};
}

constexpr Imaginary_t<double> operator""_i(unsigned long long a) {
    return Imaginary_t<double>(a);
}


#endif
