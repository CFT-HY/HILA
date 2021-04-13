#ifndef CMPLX_H_
#define CMPLX_H_

// let's not include the std::complex
//#include <complex>
//#include <cmath>

#include <sstream>
#include <iostream>
#include <math.h>
#include <type_traits>
#include <cmath>

#include "plumbing/defs.h"

/// TEMPORARY location for vector intrinsic analogues -- result obvious

template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
inline T mul_add(T a, T b, T c) {
    return a * b + c;
}

template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
inline T mul_sub(T a, T b, T c) {
    return a * b - c;
}

template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
inline T nmul_add(T a, T b, T c) {
    return c - a * b;
}

///////////////////////////////////////////////////////////////////////////////
/// main cmpx type definition
/// Define complex type as a template. This allows Hilapp to replace the internal
/// type with a vector. The datatype T must be an arithmetic type.
///////////////////////////////////////////////////////////////////////////////

template <typename T = double> struct Complex {

    static_assert(hila::is_arithmetic<T>::value,
                  "Complex can be used only with arithmetic type");
    // This incantation is needed to make Field<Complex<>> vectorized

    using base_type = hila::number_type<T>;
    using argument_type = T;

    // and the content of the complex number
    T re, im;

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

#pragma hila loop_function
    template <typename S, std::enable_if_t<hila::is_arithmetic<S>::value, int> = 0>
    explicit constexpr Complex<T>(const S val) {
        re = val;
        im = 0;
    }

    // allow construction from 0 (nullptr)
    constexpr Complex<T>(const std::nullptr_t n) { re = im = 0; }

    // constructor c(a,b)
    template <typename A, typename B,
              std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0,
              std::enable_if_t<hila::is_arithmetic<B>::value, int> = 0>
#pragma hila loop_function
    explicit constexpr Complex<T>(const A &a, const B &b) {
        re = a;
        im = b;
    }

    // make also std accessors real() and imag() - though no real difference to .re, .im
    inline T real() const { return re; }
    inline T &real() { return re; }

    inline T imag() const { return im; }
    inline T &imag() { return im; }

    // automatic casting from Complex<T> -> Complex<A>
    // TODO: ensure this works if A is vector type!
    template <typename A>
#pragma hila loop_function // TODO
    operator Complex<A>() const {
        return Complex<A>(re, im);
    }

    inline Complex<T> &operator=(const Complex<T> &s) = default;

    // Assignment from Complex<A>
    template <typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
    inline Complex<T> &operator=(const Complex<A> &s) {
        re = s.re;
        im = s.im;
        return *this;
    }

    template <typename S, std::enable_if_t<std::is_arithmetic<S>::value, int> = 0>
    inline Complex<T> &operator=(S s) {
        re = s;
        im = 0;
        return *this;
    }

    inline T squarenorm() const { return re * re + im * im; }

    // TODO: make this work for vector type!
    // currently this gives a compilation error
    inline T abs() const { return sqrt(squarenorm()); }
    inline T arg() const { return atan2(im, re); }

    inline Complex<T> conj() const { return Complex<T>(re, -im); }

    inline Complex<T> polar(const T r, const T theta) {
        return Complex<T>({r * cos(theta), r * sin(theta)});
    }

    inline Complex<T> &random() output_only {
        re = hila::random();
        im = hila::random();
        return *this;
    }

    inline Complex<T> &gaussian() output_only {
        re = hila::gaussian_ran2(im);
        return *this;
    }

    // unary + and -
    inline Complex<T> operator+() const { return *this; }
    inline Complex<T> operator-() const { return Complex<T>(-re, -im); }

    // mark += and *= as loop_functions, because these can be used
    // in reductions implicitly
#pragma hila loop_function
    inline Complex<T> &operator+=(const Complex<T> &lhs) {
        re += lhs.re;
        im += lhs.im;
        return *this;
    }

    template <typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
    inline Complex<T> &operator+=(const A &a) {
        re += a;
        return *this;
    }

    inline Complex<T> &operator-=(const Complex<T> &lhs) {
        re -= lhs.re;
        im -= lhs.im;
        return *this;
    }

    // TODO: for vector too
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

#pragma hila loop_function
    inline Complex<T> &operator*=(const Complex<T> lhs) {
        T r = mul_sub(re, lhs.re, im * lhs.im); // a*b-c
        im = mul_add(im, lhs.re, re * lhs.im);  // a*b+c
        re = r;
        return *this;
    }

    // TODO: for vector too
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
    inline Complex<T> &operator/=(const Complex<T> &lhs) {
        T n = lhs.squarenorm();
        T r = mul_add(re, lhs.re, im * lhs.im) / n; // a*b+c
        im = mul_sub(im, lhs.re, re * lhs.im) / n;  // a*b-c
        re = r;
        return *this;
    }

    // TODO: for vector too
    template <typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
    inline Complex<T> &operator/=(const A &a) {
        re /= a;
        im /= a;
        return *this;
    }

    template <typename A = T, std::enable_if_t<!std::is_arithmetic<A>::value, int> = 0>
    std::string str() const {
        std::string text = "(" + re.str() + "," + im.str() + ")";
        return text;
    }

    template <typename A = T, std::enable_if_t<std::is_arithmetic<A>::value, int> = 0>
    std::string str() const {
        std::string text = "(" + std::to_string(re) + "," + std::to_string(im) + ")";
        return text;
    }
};

// functions real(), imag()

template <typename T> inline T real(const Complex<T> a) { return a.re; }

template <typename T> inline T imag(const Complex<T> a) { return a.im; }

// template <typename T>
// inline Complex<T> operator+(const Complex<T> & a, const Complex<T> & b) {
//   return Complex<T>(a.re + b.re, a.im + b.im);
// }

template <typename T> inline Complex<T> operator+(Complex<T> a, const Complex<T> &b) {
    a += b;
    return a;
}

// TODO: for avx vector too -- #define new template macro
template <typename T, typename A,
          std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline Complex<T> operator+(const Complex<T> &c, const A &a) {
    return Complex<T>(c.re + a, c.im);
}

template <typename T, typename A,
          std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline Complex<T> operator+(const A &a, const Complex<T> &c) {
    return Complex<T>(c.re + a, c.im);
}

// -
// template <typename T>
// inline Complex<T> operator-(const Complex<T> & a, const Complex<T> & b) {
//   return Complex<T>(a.re - b.re, a.im - b.im);
// }
template <typename T> inline Complex<T> operator-(Complex<T> a, const Complex<T> &b) {
    a -= b;
    return a;
}

// TODO: for avx vector too -- #define new template macro
template <typename T, typename A,
          std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline Complex<T> operator-(const Complex<T> &c, const A &a) {
    return Complex<T>(c.re - a, c.im);
}

template <typename T, typename A,
          std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline Complex<T> operator-(const A &a, const Complex<T> &c) {
    return Complex<T>(a - c.re, -c.im);
}

//
// template <typename T>
// inline Complex<T> operator*(const Complex<T> & a, const Complex<T> & b) {
//   return Complex<T>(a.re*b.re - a.im*b.im, a.im*b.re + a.re*b.im);
// }
template <typename T> inline Complex<T> operator*(Complex<T> a, const Complex<T> &b) {
    a *= b;
    return a;
}

template <typename T, typename A,
          std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline Complex<T> operator*(const Complex<T> &c, const A &a) {
    return Complex<T>(c.re * a, c.im * a);
}

template <typename T, typename A,
          std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline Complex<T> operator*(const A &a, const Complex<T> &c) {
    return Complex<T>(a * c.re, a * c.im);
}

// /   a/b = ab*/|b|^2
// template <typename T>
// inline Complex<T> operator/(const Complex<T> & a, const Complex<T> & b) {
//   T n = b.squarenorm();
//   return Complex<T>( (a.re*b.re + a.im*b.im)/n, (a.im*b.re - a.re*b.im)/n );
// }
template <typename T> inline Complex<T> operator/(Complex<T> a, const Complex<T> &b) {
    a /= b;
    return a;
}

template <typename T, typename A,
          std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline Complex<T> operator/(const Complex<T> &c, const A &a) {
    return Complex<T>(c.re / a, c.im / a);
}

// a/c = ac*/|c|^2
template <typename T, typename A,
          std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline Complex<T> operator/(const A &a, const Complex<T> &c) {
    T n = c.squarenorm();
    return Complex<T>((a * c.re) / n, -(a * c.im) / n);
}

// write also multiply-add directly with complex numbers
template <typename T>
inline Complex<T> mul_add(const Complex<T> &a, const Complex<T> &b,
                          const Complex<T> &c) {
    // a*b + c
    Complex<T> r;
    T t1 = mul_add(a.re, b.re, c.re);
    T t2 = mul_add(a.re, b.im, c.im);
    r.re = nmul_add(a.im, b.im, t1); // -a.im*b.im + a.re*b.re + c.re
    r.im = mul_add(a.im, b.re, t2);  // a.im*b.re + a.re*b.im + c.im
    return r;
}

//////////////////////////////////////////////////////////////////////////////////
// Some operations in function form.  Useful in templates when the arg type is not known

/// abs
template <typename T> inline T abs(const Complex<T> &a) { return a.abs(); }

/// arg
template <typename T> inline T arg(const Complex<T> &a) { return a.arg(); }

/// Conjugate
template <typename T> inline Complex<T> conj(const Complex<T> &val) {
    return val.conj();
}

/// norm_squared
template <typename T> inline auto squarenorm(const Complex<T> &val) {
    return val.squarenorm();
}


/// random() : set argument to random vals [0,1]
template <typename T> inline void random(Complex<T> &c) {
    ::random(c.re);
    ::random(c.im);
}

template <typename T> inline void gaussian_random(Complex<T> &c) {
    gaussian_random(c.re);
    gaussian_random(c.im);
}

//////////////////////////////////////////////////////////////////////////////////
/// Print a complex value as (re,im)
//////////////////////////////////////////////////////////////////////////////////

template <typename T>
std::ostream &operator<<(std::ostream &strm, const Complex<T> &A) {
    return strm << "(" << A.re << ", " << A.im << ")";
}

//////////////////////////////////////////////////////////////////////////////////
/// Operators to implement imaginary unit 1_i, enablig expressions  3 + 2_i  etc.
/// This is defined as an user-defined literal, which requires an underscore.
////////////////////////////////////////////////////////////////////////////////

constexpr Complex<double> operator""_i(long double a) {
    return Complex<double>{0.0, a};
}

constexpr Complex<double> operator""_i(unsigned long long a) {
    return Complex<double>(0.0, static_cast<double>(a));
}

//////////////////////////////////////////////////////////////////////////

// define also real(), imag(), conj() -functions for basic arithmetic types
template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline T real(T val) {
    return val;
}

template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline T imag(T val) {
    return 0;
}

template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline T conj(T val) {
    return val;
}

///////////////////////////////////////////////////////////////////////////////
/// Set of complex functions
///////////////////////////////////////////////////////////////////////////////

/// multiply_by_i(z)
/// an auxiliary function, less operations than with 1_i * z == (0,1) * z
template <typename T> inline Complex<T> multiply_by_i(Complex<T> z) {
    return Complex<T>(-z.im, z.re);
}

/// exp(z)
template <typename T> inline Complex<T> exp(const Complex<T> z) {
    return exp(z.re) * Complex<T>(cos(z.im), sin(z.im));
}

/// exp(i x)
template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline Complex<T> expi(T a) {
    return Complex<T>(cos(a), sin(a));
}

/// log(z)
template <typename T> inline Complex<T> log(Complex<T> z) {
    return Complex<T>(static_cast<T>(0.5) * log(z.squarenorm()), z.arg());
}

/// sqrt(z) branch cut at -x axis
template <typename T> inline Complex<T> sqrt(Complex<T> z) {
    T r = z.squarenorm();
    T a = z.arg();
    return pow(r,0.25) * expi(0.5*a);
}

/// cbrt(z)
template <typename T> inline Complex<T> cbrt(Complex<T> z) {
    T r = z.squarenorm();
    T a = z.arg();
    return pow(r,(1.0/6.0)) * expi((1.0/3.0)*a);
}


/// pow(z.p) = z^p = exp(p*log(z))
template <typename T> inline Complex<T> pow(Complex<T> z, Complex<T> p) {
    return exp(p * log(z));
}

/// pow(z.p) with scalar power
template <typename T> inline Complex<T> pow(Complex<T> z, T p) {
    return exp(p * log(z));
}

/// pow(z.p) with scalar base
template <typename T> inline Complex<T> pow(T z, Complex<T> p) {
    return exp(p * log(z));
}


/// sin(z)
/// = sin(re + i im) = sin(re)cos(i im) + cos(re)sin(i im)
/// = sin(re) cosh(im) + i cos(re) sinh(im)
template <typename T> inline Complex<T> sin(Complex<T> z) {
    return Complex<T>(sin(z.re) * cosh(z.im), cos(z.re) * sinh(z.im));
}

/// cos(z)
/// = cos(re)cos(i im) - sin(re)sin(i im) = cos(re)cosh(im) - i sin(re)sinh(im)
template <typename T> inline Complex<T> cos(Complex<T> z) {
    return Complex<T>(cos(z.re) * cosh(z.im), -sin(z.re) * sinh(z.im));
}

/// tan(z) - rely on optimizer to simplify
template <typename T> inline Complex<T> tan(Complex<T> z) { return sin(z) / cos(z); }

/// sinh(z) = sinh(re)cosh(i im) + cosh(re)sinh(i im)
/// = sinh(re)cos(im) + i cosh(re)sin(im)
template <typename T> inline Complex<T> sinh(Complex<T> z) {
    return Complex<T>(sinh(z.re) * cos(z.im), cosh(z.re) * sin(z.im));
}

/// cosh(z) = cosh(re)cosh(i im) - sinh(re)sinh(i im)
/// = cosh(re)cos(im) - i sinh(re)sin(im)
template <typename T> inline Complex<T> cosh(Complex<T> z) {
    return Complex<T>(cosh(z.re) * cos(z.im), sinh(z.re) * sin(z.im));
}

/// tanh(z)
template <typename T> inline Complex<T> tanh(Complex<T> z) {
     return sinh(z) / cosh(z);
}

/// arctan(z)
template <typename T> inline Complex<T> atan(Complex<T> z) {
    return -0.5 * multiply_by_i(log((1_i - z) / (1_i + z)));
}

/// arcsin(z)
template <typename T> inline Complex<T> asin(Complex<T> z) {
    return -multiply_by_i(log(multiply_by_i(z) + sqrt(1 - z * z)));
}

/// arccos(z)
template <typename T> inline Complex<T> acos(Complex<T> z) {
    return -multiply_by_i(log(z + multiply_by_i(sqrt(1 - z * z))));
}

/// artanh(z)
template <typename T> inline Complex<T> atanh(Complex<T> z) {
    return 0.5*log( (1+z)/(1-z) );
}

/// arsinh(z)
template <typename T> inline Complex<T> asinh(Complex<T> z) {
    return log( z + sqrt(1 + z*z));
}

/// arcosh(z)
template <typename T> inline Complex<T> acosh(Complex<T> z) {
    return log( z + sqrt(z*z - 1));
}







namespace hila {
////////////////////////////////////////////////////////////////////////
/// And utility templates
/// Define is_complex<T>::value -template, using specialization
template <typename T> struct is_complex : std::integral_constant<bool, false> {};

template <typename T>
struct is_complex<Complex<T>> : std::integral_constant<bool, true> {};

// and a template is_complex_or_arithmetic<T>::value
template <typename T>
struct is_complex_or_arithmetic
    : std::integral_constant<bool, hila::is_arithmetic<T>::value ||
                                       hila::is_complex<T>::value> {};

/// Utility to check that the type contains complex numbers
/// Use as contains_complex<T>::value
template <typename T>
using contains_complex = hila::contains_type<T, Complex<hila::number_type<T>>>;

} // namespace hila

#endif
