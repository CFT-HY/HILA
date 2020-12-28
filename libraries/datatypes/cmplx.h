#ifndef CMPLX_H
#define CMPLX_H

// let's not include the std::complex
//#include <complex>
//#include <cmath>

#include <sstream>
#include <iostream>
#include <math.h>
#include <type_traits>
// #include "plumbing/defs.h"

/// TEMPORARY location for vector intrinsic analogues -- result obvious

template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0 >
#pragma hila loop_function
inline T mul_add(T a, T b, T c) { return a*b + c; }

template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0 >
#pragma hila loop_function
inline T mul_sub(T a, T b, T c) { return a*b - c; }

template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0 >
#pragma hila loop_function
inline T nmul_add(T a, T b, T c) { return c - a*b; }


/// Define complex type as a template. This allows Hilapp to replace the internal
/// type with a vector.
template <typename T = double>
struct cmplx {
  
  // This incantation is needed to make field<cmplx<>> vectorized 
  using base_type = typename base_type_struct<T>::type;

  T re,im;
  
  cmplx<T>() = default;
  
  cmplx<T>(const cmplx<T> & a) =default;

  // constructor from single complex 
  template <typename A,
            std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
  #pragma hila loop_function
  constexpr cmplx<T>(const cmplx<A> a): re(static_cast<T>(a.re)), im(static_cast<T>(a.im)) {}

  // constructor from single scalar value 
  template <typename scalar_t,
            std::enable_if_t<is_arithmetic<scalar_t>::value, int> = 0 >
  #pragma hila loop_function
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
  #pragma hila loop_function
  constexpr cmplx<T>(const A & a, const B & b): re(static_cast<T>(a)), im(static_cast<T>(b)) {}

  ~cmplx<T>() =default;
  
  // automatic casting from cmplx<T> -> cmplx<A>
  // TODO: ensure this works if A is vector type!
  template <typename A>
  #pragma hila loop_function
  operator cmplx<A>() const { 
    return cmplx<A>( static_cast<A>(re), static_cast<A>(im) );
  }

  // Assignment from cmplx<A>
  template <typename A,
            std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
  #pragma hila loop_function
  inline cmplx<T> & operator=(cmplx<A> s) {
    re = static_cast<T>(s.re);
    im = static_cast<T>(s.im);
    return *this;
  }
  
  template <typename scalar_t,
            std::enable_if_t<is_arithmetic<scalar_t>::value, int> = 0 >
  #pragma hila loop_function
  inline cmplx<T> & operator=(scalar_t s) {
    re = static_cast<T>(s);
    im = static_cast<T>(0);
    return *this;
  }
  
  #pragma hila loop_function
  inline T real() const { return re; }
  #pragma hila loop_function
  inline T imag() const { return im; }

  #pragma hila loop_function
  inline T squarenorm() const { return re*re + im*im; }
  #pragma hila loop_function
  inline T norm_sq() const { return re*re + im*im; }

  // TODO: make this work for vector type!  Not double  
  //currently this gives a compilation error
  #pragma hila loop_function
  inline double abs() const { return sqrt(static_cast<double>(squarenorm()) ); }
  #pragma hila loop_function
  inline double arg() const { return atan2(static_cast<double>(im),static_cast<double>(re)); }


  #pragma hila loop_function
  inline cmplx<T> conj() const { return cmplx<T>( { re, -im } ); }

  #pragma hila loop_function
  inline cmplx<T> polar(const T r, const T theta) { 
    return cmplx<T>( { r*cos(theta), r*sin(theta) } );
  }

  #pragma hila loop_function
  template <typename A=T, std::enable_if_t<std::is_arithmetic<A>::value, int> = 0 > 
  cmplx<A> & random(){
    re = static_cast<T>(hila_random());
    im = static_cast<T>(hila_random());
    return *this;
  }

  // unary + and -
  #pragma hila loop_function
  inline cmplx<T> operator+() const {return *this;}
  #pragma hila loop_function
  inline cmplx<T> operator-() const {return cmplx<T>(-re, -im); }


  #pragma hila loop_function
  inline cmplx<T> & operator+= (const cmplx<T> & lhs) {
    re += lhs.re;
    im += lhs.im;
    return *this;
  }

  // TODO: for avx vector too -- #define new template macro
  template <typename A,
            std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
  #pragma hila loop_function
  inline cmplx<T> & operator+= (const A & a) {
    re += static_cast<T>(a);
    return *this;
  }

  #pragma hila loop_function
  inline cmplx<T> & operator-= (const cmplx<T> & lhs) {
    re -= lhs.re;
    im -= lhs.im;
    return *this;
  }
  
  // TODO: for vector too
  template <typename A,
            std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
  #pragma hila loop_function
  inline cmplx<T> & operator-= (const A & a) {
    re -= static_cast<T>(a);
    return *this;
  }
  
  // #pragma hila loop_function
  // inline cmplx<T> & operator*= (const cmplx<T> & lhs) {
  //   T r = re * lhs.re - im * lhs.im;
  //   im  = im * lhs.re + re * lhs.im;
  //   re = r;
  //   return *this;
  // }
  
  #pragma hila loop_function
  inline cmplx<T> & operator*= (const cmplx<T> lhs) {
    T r = mul_sub(re, lhs.re, im * lhs.im);   // a*b-c
    im  = mul_add(im, lhs.re, re * lhs.im);   // a*b+c
    re = r;
    return *this;
  }
  

  // TODO: for vector too
  template <typename A,
            std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
  #pragma hila loop_function
  inline cmplx<T> & operator*= (const A a) {
    re *= static_cast<T>(a);
    im *= static_cast<T>(a);
    return *this;
  }

  // a/b = a b*/|b|^2 = (a.re*b.re + a.im*b.im + i(a.im*b.re - a.re*b.im))/|b|^2
  // #pragma hila loop_function
  // inline cmplx<T> & operator/= (const cmplx<T> & lhs) {
  //   T n = lhs.squarenorm();
  //   T r = (re * lhs.re + im * lhs.im)/n;
  //   im  = (im * lhs.re - re * lhs.im)/n;
  //   re = r;
  //   return *this;
  // }
  #pragma hila loop_function
  inline cmplx<T> & operator/= (const cmplx<T> & lhs) {
    T n = lhs.squarenorm();
    T r = mul_add(re, lhs.re, im * lhs.im)/n;  // a*b+c
    im  = mul_sub(im, lhs.re, re * lhs.im)/n;  // a*b-c
    re = r;
    return *this;
  }
  
  // TODO: for vector too
  template <typename A,
            std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
  #pragma hila loop_function
  inline cmplx<T> & operator/= (const A & a) {
    re /= static_cast<T>(a);
    im /= static_cast<T>(a);
    return *this;
  }

  template <typename A = T,
            std::enable_if_t<!std::is_arithmetic<A>::value, int> = 0 >
  std::string str() const {
    std::string text = "(" + re.str() + "," + im.str() + ")"; 
    return text;
  }

  template <typename A = T,
            std::enable_if_t<std::is_arithmetic<A>::value, int> = 0 >
  std::string str() const {
    std::string text = "(" + std::to_string(re) + "," + std::to_string(im) + ")"; 
    return text;
  }

};

// template <typename T>
// #pragma hila loop_function
// inline cmplx<T> operator+(const cmplx<T> & a, const cmplx<T> & b) {
//   return cmplx<T>(a.re + b.re, a.im + b.im);
// }

template <typename T>
#pragma hila loop_function
inline cmplx<T> operator+(cmplx<T> a, const cmplx<T> & b) {
  a += b;
  return a;
}

  // TODO: for avx vector too -- #define new template macro
template <typename T, typename A,
          std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
#pragma hila loop_function
inline cmplx<T> operator+(const cmplx<T> & c, const A & a) {
  return cmplx<T>(c.re + a, c.im);
}

template <typename T, typename A,
          std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
#pragma hila loop_function
inline cmplx<T> operator+(const A &a, const cmplx<T> & c) {
  return cmplx<T>(c.re + a, c.im);
}

// -
// template <typename T>
// #pragma hila loop_function
// inline cmplx<T> operator-(const cmplx<T> & a, const cmplx<T> & b) {
//   return cmplx<T>(a.re - b.re, a.im - b.im);
// }
template <typename T>
#pragma hila loop_function
inline cmplx<T> operator-(cmplx<T> a, const cmplx<T> & b) {
  a -= b;
  return a;
}

// TODO: for avx vector too -- #define new template macro
template <typename T, typename A,
          std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
#pragma hila loop_function
inline cmplx<T> operator-(const cmplx<T> & c, const A & a) {
  return cmplx<T>(c.re - a, c.im);
}

template <typename T, typename A,
          std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
#pragma hila loop_function
inline cmplx<T> operator-(const A &a, const cmplx<T> & c) {
  return cmplx<T>(a - c.re, -c.im);
}


// *
// template <typename T>
// #pragma hila loop_function
// inline cmplx<T> operator*(const cmplx<T> & a, const cmplx<T> & b) {
//   return cmplx<T>(a.re*b.re - a.im*b.im, a.im*b.re + a.re*b.im);
// }
template <typename T>
#pragma hila loop_function
inline cmplx<T> operator*(cmplx<T> a, const cmplx<T> & b) {
  a *= b;
  return a;
}

// TODO: for avx vector too -- #define new template macro
template <typename T, typename A,
          std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
#pragma hila loop_function
inline cmplx<T> operator*(const cmplx<T> & c, const A & a) {
  return cmplx<T>(c.re * a, c.im * a);
}

template <typename T, typename A,
          std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
#pragma hila loop_function
inline cmplx<T> operator*(const A &a, const cmplx<T> & c) {
  return cmplx<T>(a * c.re, a * c.im);
}


// /   a/b = ab*/|b|^2
// template <typename T>
// #pragma hila loop_function
// inline cmplx<T> operator/(const cmplx<T> & a, const cmplx<T> & b) {
//   T n = b.squarenorm();
//   return cmplx<T>( (a.re*b.re + a.im*b.im)/n, (a.im*b.re - a.re*b.im)/n );
// }
template <typename T>
#pragma hila loop_function
inline cmplx<T> operator/(cmplx<T> a, const cmplx<T> & b) {
  a /= b;
  return a;
}

// TODO: for avx vector too -- #define new template macro
template <typename T, typename A,
          std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
#pragma hila loop_function
inline cmplx<T> operator/(const cmplx<T> & c, const A & a) {
  return cmplx<T>(c.re / a, c.im / a);
}

// a/c = ac*/|c|^2
template <typename T, typename A,
          std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
#pragma hila loop_function
inline cmplx<T> operator/(const A &a, const cmplx<T> & c) {
  T n = c.squarenorm();
  return cmplx<T>((a * c.re)/n, -(a * c.im)/n);
}

// write also multiply-add directly with complex numbers
template <typename T>
#pragma hila loop_function
inline cmplx<T> mul_add(const cmplx<T> &a, const cmplx<T>  &b, const cmplx<T> &c ) {
  //a*b + c
  cmplx<T> r;
  T t1 = mul_add(a.re , b.re, c.re);
  T t2 = mul_add(a.re , b.im, c.im);
  r.re = nmul_add( a.im, b.im, t1 );  // -a.im*b.im + a.re*b.re + c.re
  r.im = mul_add(  a.im, b.re, t2 ); // a.im*b.re + a.re*b.im + c.im
  return r;
}
 

// Operators to implement imaginary unit 1_i, enablig expressions  3 + 2_i  etc.
// Underscore seems to be required here
#pragma hila loop_function
constexpr cmplx<double> operator""_i(long double a) {
  return cmplx<double>{0.0,a};
}

#pragma hila loop_function
constexpr cmplx<double> operator""_i(unsigned long long a) {
  return cmplx<double>(0.0,static_cast<double>(a));
}

template <typename T>
#pragma hila loop_function
std::ostream& operator<<(std::ostream &strm, const cmplx<T> A) {
  return strm << "(" << A.re << ", " << A.im << ")";
}

template<typename Accuracy> 
#pragma hila loop_function
inline cmplx<Accuracy> conj(cmplx<Accuracy> val){
  return val.conj();
}

template<typename T> 
#pragma hila loop_function
inline auto norm_squared(cmplx<T> val){
  return val.squarenorm();
}


#endif

