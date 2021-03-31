#ifndef CMPLX_H_
#define CMPLX_H_

// let's not include the std::complex
//#include <complex>
//#include <cmath>

#include <sstream>
#include <iostream>
#include <math.h>
#include <type_traits>

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

    static_assert(is_arithmetic<T>::value,
                  "Complex can be used only with arithmetic type");
    // This incantation is needed to make Field<Complex<>> vectorized
 
    using base_type = number_type<T>;
    constexpr bool complex_base = true;

    // and the content of the complex number
    T re, im;

    Complex<T>() = default;
    ~Complex<T>() = default;
    Complex<T>(const Complex<T> &a) = default;

    // constructor from single complex --IS THIS NEEDED?
    // template <typename A,
    //           std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
    // constexpr Complex<T>(const Complex<A> a) : re(static_cast<T>(a.re)),
    // im(static_cast<T>(a.im)) {}

    // constructor from single scalar value
    // Remember to mark this explicit, we do not want this to be invoked
    // in automatic conversions (there should be methods)

#pragma hila loop_function
    template <typename S, std::enable_if_t<is_arithmetic<S>::value, int> = 0>
    explicit constexpr Complex<T>(const S val) {
        re = val;
        im = 0;
    }

    // allow construction from 0 (nullptr)
    constexpr Complex<T>(const std::nullptr_t n) { re = im = 0; }

    // constructor c(a,b)
    template <typename A, typename B,
              std::enable_if_t<is_arithmetic<A>::value, int> = 0,
              std::enable_if_t<is_arithmetic<B>::value, int> = 0>
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
    template <typename A, std::enable_if_t<is_arithmetic<A>::value, int> = 0>
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

    inline T norm_sq() const { return re * re + im * im; }

    // TODO: make this work for vector type!
    // currently this gives a compilation error
    inline T abs() const { return sqrt(norm_sq()); }
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

    template <typename A, std::enable_if_t<is_arithmetic<A>::value, int> = 0>
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
    template <typename A, std::enable_if_t<is_arithmetic<A>::value, int> = 0>
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
    template <typename A, std::enable_if_t<is_arithmetic<A>::value, int> = 0>
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
        T n = lhs.norm_sq();
        T r = mul_add(re, lhs.re, im * lhs.im) / n; // a*b+c
        im = mul_sub(im, lhs.re, re * lhs.im) / n;  // a*b-c
        re = r;
        return *this;
    }

    // TODO: for vector too
    template <typename A, std::enable_if_t<is_arithmetic<A>::value, int> = 0>
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
template <typename T, typename A, std::enable_if_t<is_arithmetic<A>::value, int> = 0>
inline Complex<T> operator+(const Complex<T> &c, const A &a) {
    return Complex<T>(c.re + a, c.im);
}

template <typename T, typename A, std::enable_if_t<is_arithmetic<A>::value, int> = 0>
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
template <typename T, typename A, std::enable_if_t<is_arithmetic<A>::value, int> = 0>
inline Complex<T> operator-(const Complex<T> &c, const A &a) {
    return Complex<T>(c.re - a, c.im);
}

template <typename T, typename A, std::enable_if_t<is_arithmetic<A>::value, int> = 0>
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

template <typename T, typename A, std::enable_if_t<is_arithmetic<A>::value, int> = 0>
inline Complex<T> operator*(const Complex<T> &c, const A &a) {
    return Complex<T>(c.re * a, c.im * a);
}

template <typename T, typename A, std::enable_if_t<is_arithmetic<A>::value, int> = 0>
inline Complex<T> operator*(const A &a, const Complex<T> &c) {
    return Complex<T>(a * c.re, a * c.im);
}

// /   a/b = ab*/|b|^2
// template <typename T>
// inline Complex<T> operator/(const Complex<T> & a, const Complex<T> & b) {
//   T n = b.norm_sq();
//   return Complex<T>( (a.re*b.re + a.im*b.im)/n, (a.im*b.re - a.re*b.im)/n );
// }
template <typename T> inline Complex<T> operator/(Complex<T> a, const Complex<T> &b) {
    a /= b;
    return a;
}

template <typename T, typename A, std::enable_if_t<is_arithmetic<A>::value, int> = 0>
inline Complex<T> operator/(const Complex<T> &c, const A &a) {
    return Complex<T>(c.re / a, c.im / a);
}

// a/c = ac*/|c|^2
template <typename T, typename A, std::enable_if_t<is_arithmetic<A>::value, int> = 0>
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
template <typename T> inline auto norm_squared(const Complex<T> &val) {
    return val.norm_sq();
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
template <typename T, std::enable_if_t<is_arithmetic<T>::value, int> = 0>
inline T real(T val) {
    return val;
}

template <typename T, std::enable_if_t<is_arithmetic<T>::value, int> = 0>
inline T imag(T val) {
    return 0;
}

template <typename T, std::enable_if_t<is_arithmetic<T>::value, int> = 0>
inline T conj(T val) {
    return val;
}

////////////////////////////////////////////////////////////////////////
/// And utility templates
/// Define is_cmplx<T>::value -template, using specialization
template <typename T> struct is_cmplx : std::integral_constant<bool, false> {};

template <typename T>
struct is_cmplx<Complex<T>> : std::integral_constant<bool, true> {};

// and a template is_cmplx_or_real<T>::value
template <typename T>
struct is_cmplx_or_arithmetic
    : std::integral_constant<bool, is_arithmetic<T>::value || is_cmplx<T>::value> {};

#endif
