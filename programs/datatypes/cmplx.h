#ifndef CMPLX_H
#define CMPLX_H

// let's not include the std::complex
//#include <complex>
//#include <cmath>

#include <iostream>
#include <math.h>
#include <type_traits>
#include "../plumbing/defs.h"

template <typename T = double>
struct cmplx {
  T re,im;
  
  cmplx<T>() = default;
  
  cmplx<T>(const cmplx<T> & a) =default;

  // constructor from single complex 
  template <typename A,
            std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
  #pragma transformer loop_function
  constexpr cmplx<T>(const cmplx<A> a): re(static_cast<T>(a.re)), im(static_cast<T>(a.im)) {}

  // constructor from single scalar value 
  template <typename scalar_t,
            std::enable_if_t<is_arithmetic<scalar_t>::value, int> = 0 >
  #pragma transformer loop_function
  constexpr cmplx<T>(const scalar_t val): re(static_cast<T>(val)), im(static_cast<T>(0)) {}

  // constructor c(a,b)
//   template <typename A, typename B,
//             std::enable_if_t<is_arithmetic<A>::value, int> = 0,
//             std::enable_if_t<is_arithmetic<B>::value, int> = 0 >
//   constexpr cmplx<T>(const A & a, const B & b) {
//     re = static_cast<T>(a);
//     im = static_cast<T>(b);
//   }

  // constructor c(a,b)
  template <typename A, typename B,
            std::enable_if_t<is_arithmetic<A>::value, int> = 0,
            std::enable_if_t<is_arithmetic<B>::value, int> = 0 >
  #pragma transformer loop_function
  constexpr cmplx<T>(const A & a, const B & b): re(static_cast<T>(a)), im(static_cast<T>(b)) {}

  ~cmplx<T>() =default;
  
  // automatic casting from cmplx<T> -> cmplx<A>
  // TODO: ensure this works if A is vector type!
  template <typename A>
  #pragma transformer loop_function
  operator cmplx<A>() const { 
    return cmplx<A>( static_cast<A>(re), static_cast<A>(im) );
  }

  // Assignment from cmplx<A>
  template <typename A,
            std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
  #pragma transformer loop_function
  cmplx<T> & operator=(cmplx<A> s) {
    re = static_cast<T>(s.re);
    im = static_cast<T>(s.im);
    return *this;
  }
  
  template <typename scalar_t,
            std::enable_if_t<is_arithmetic<scalar_t>::value, int> = 0 >
  #pragma transformer loop_function
  cmplx<T> & operator=(scalar_t s) {
    re = static_cast<T>(s);
    im = 0.0;
    return *this;
  }
  
  #pragma transformer loop_function
  T real() const { return re; }
  #pragma transformer loop_function
  T imag() const { return im; }

  #pragma transformer loop_function
  T norm() const { return re*re + im*im; }
  // TODO: make this work for vector type!  Not double
  
  //currently this gives a compilation error
  #pragma transformer loop_function
  double abs() const { return sqrt(static_cast<double>(norm()) ); }
  #pragma transformer loop_function
  double arg() const { return atan2(static_cast<double>(im),static_cast<double>(re)); }


  #pragma transformer loop_function
  cmplx<T> conj() const { return cmplx<T>( { re, -im } ); }

  #pragma transformer loop_function
  cmplx<T> polar(const T r, const T theta) { 
    return cmplx<T>( { r*cos(theta), r*sin(theta) } );
  }

  // unary + and -
  #pragma transformer loop_function
  cmplx<T> operator+() const {return *this;}
  #pragma transformer loop_function
  cmplx<T> operator-() const {return cmplx<T>(-re, -im); }


  #pragma transformer loop_function
  cmplx<T> & operator+= (const cmplx<T> & lhs) {
    re += lhs.re;
    im += lhs.im;
    return *this;
  }

  // TODO: for avx vector too -- #define new template macro
  template <typename A,
            std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
  #pragma transformer loop_function
  cmplx<T> & operator+= (const A & a) {
    re += static_cast<T>(a);
    return *this;
  }

  #pragma transformer loop_function
  cmplx<T> & operator-= (const cmplx<T> & lhs) {
    re -= lhs.re;
    im -= lhs.im;
    return *this;
  }
  
  // TODO: for vector too
  template <typename A,
            std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
  #pragma transformer loop_function
  cmplx<T> & operator-= (const A & a) {
    re -= static_cast<T>(a);
    return *this;
  }
  
  #pragma transformer loop_function
  cmplx<T> & operator*= (const cmplx<T> & lhs) {
    re = re * lhs.re - im * lhs.im;
    im = im * lhs.re + re * lhs.im;
    return *this;
  }
  
  // TODO: for vector too
  template <typename A,
            std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
  #pragma transformer loop_function
  cmplx<T> & operator*= (const A & a) {
    re *= static_cast<T>(a);
    im *= static_cast<T>(a);
    return *this;
  }

  // a/b = a b*/|b|^2 = (a.re*b.re + a.im*b.im + i(a.im*b.re - a.re*b.im))/|b|^2
  #pragma transformer loop_function
  cmplx<T> & operator/= (const cmplx<T> & lhs) {
    T n = lhs.norm();
    re = (re * lhs.re + im * lhs.im)/n;
    im = (im * lhs.re - re * lhs.im)/n;
    return *this;
  }
  
  // TODO: for vector too
  template <typename A,
            std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
  #pragma transformer loop_function
  cmplx<T> & operator/= (const A & a) {
    re /= static_cast<T>(a);
    im /= static_cast<T>(a);
    return *this;
  }

  #pragma transformer loop_function
  auto norm_sq(){
    return re*re + im*im;
  }
};

template <typename T>
#pragma transformer loop_function
cmplx<T> operator+(const cmplx<T> & a, const cmplx<T> & b) {
  return cmplx<T>(a.re + b.re, a.im + b.im);
}

  // TODO: for avx vector too -- #define new template macro
template <typename T, typename A,
          std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
#pragma transformer loop_function
cmplx<T> operator+(const cmplx<T> & c, const A & a) {
  return cmplx<T>(c.re + a, c.im);
}

template <typename T, typename A,
          std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
#pragma transformer loop_function
cmplx<T> operator+(const A &a, const cmplx<T> & c) {
  return cmplx<T>(c.re + a, c.im);
}

// -
template <typename T>
#pragma transformer loop_function
cmplx<T> operator-(const cmplx<T> & a, const cmplx<T> & b) {
  return cmplx<T>(a.re - b.re, a.im - b.im);
}

// TODO: for avx vector too -- #define new template macro
template <typename T, typename A,
          std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
#pragma transformer loop_function
cmplx<T> operator-(const cmplx<T> & c, const A & a) {
  return cmplx<T>(c.re - a, c.im);
}

template <typename T, typename A,
          std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
#pragma transformer loop_function
cmplx<T> operator-(const A &a, const cmplx<T> & c) {
  return cmplx<T>(a - c.re, -c.im);
}


// *
template <typename T>
#pragma transformer loop_function
cmplx<T> operator*(const cmplx<T> & a, const cmplx<T> & b) {
  return cmplx<T>(a.re*b.re - a.im*b.im, a.im*b.re + a.re*b.im);
}

// TODO: for avx vector too -- #define new template macro
template <typename T, typename A,
          std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
#pragma transformer loop_function
cmplx<T> operator*(const cmplx<T> & c, const A & a) {
  return cmplx<T>(c.re * a, c.im * a);
}

template <typename T, typename A,
          std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
#pragma transformer loop_function
cmplx<T> operator*(const A &a, const cmplx<T> & c) {
  return cmplx<T>(a * c.re, a * c.im);
}


// /   a/b = ab*/|b|^2
template <typename T>
#pragma transformer loop_function
cmplx<T> operator/(const cmplx<T> & a, const cmplx<T> & b) {
  T n = b.norm();
  return cmplx<T>( (a.re*b.re + a.im*b.im)/n, (a.im*b.re - a.re*b.im)/n );
}

// TODO: for avx vector too -- #define new template macro
template <typename T, typename A,
          std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
#pragma transformer loop_function
cmplx<T> operator/(const cmplx<T> & c, const A & a) {
  return cmplx<T>(c.re / a, c.im / a);
}

// a/c = ac*/|c|^2
template <typename T, typename A,
          std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
#pragma transformer loop_function
cmplx<T> operator/(const A &a, const cmplx<T> & c) {
  T n = c.norm();
  return cmplx<T>((a * c.re)/n, -(a * c.im)/n);
}


// Operators to implement imaginary unit 1_i, enablig expressions  3 + 2_i  etc.
// Underscore seems to be required here
#pragma transformer loop_function
constexpr cmplx<double> operator""_i(long double a) {
  return cmplx<double>{0.0,a};
}

#pragma transformer loop_function
constexpr cmplx<double> operator""_i(unsigned long long a) {
  return cmplx<double>(0.0,static_cast<double>(a));
}

template <typename T>
#pragma transformer loop_function
std::ostream& operator<<(std::ostream &strm, const cmplx<T> A) {
  return strm << "(" << A.re << ", " << A.im << ")";
}

template<typename Accuracy> 
#pragma transformer loop_function
inline cmplx<Accuracy> conj(cmplx<Accuracy> val){
  return val.conj();
}

template<typename T> 
#pragma transformer loop_function
inline auto norm_sq(cmplx<T> val){
  return val.norm_sq();
}


#endif

