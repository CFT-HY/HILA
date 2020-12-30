#ifndef CMPLX_H_
#define CMPLX_H_

// let's not include the std::complex
//#include <complex>
//#include <cmath>

#include <sstream>
#include <iostream>
#include <math.h>
#include <type_traits>
// #include "plumbing/defs.h"

#include "datatypes/zero.h"

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




///////////////////////////////////////////////////////////////////////////////
/// main cmpx type definition
///////////////////////////////////////////////////////////////////////////////

template <typename T = double>
struct cmplx {
  
  static_assert( is_arithmetic<T>::value, "Cmplx can be used only with arithmetic type" );
  // This incantation is needed to make field<cmplx<>> vectorized 
  using base_type = typename base_type_struct<T>::type;

  T re,im;
  
  cmplx<T>() = default;
  ~cmplx<T>() =default;
  cmplx<T>(const cmplx<T> & a) =default;

  // constructor from single complex --IS THIS NEEDED?
  // template <typename A,
  //           std::enable_if_t<is_arithmetic<A>::value, int> = 0 >
  // #pragma hila loop_function
  // constexpr cmplx<T>(const cmplx<A> a) : re(static_cast<T>(a.re)), im(static_cast<T>(a.im)) {}

  // constructor from single scalar value 
  template <typename scalar_t,
            std::enable_if_t<is_arithmetic<scalar_t>::value, int> = 0 >
  #pragma hila loop_function
  constexpr cmplx<T>(const scalar_t val) : re(static_cast<T>(val)), im(static_cast<T>(0)) {}

  // make zero constructor and assignment 
  constexpr cmplx<T>(const Zero z) : re(static_cast<T>(0)), im(static_cast<T>(0)) {}
  constexpr cmplx<T> & operator=(const Zero z) { 
    re = im = static_cast<T>(0);
    return *this;
  }



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


  // make also std accessors real() and imag() - though no real difference to .re, .im
  #pragma hila loop_function
  inline T real() const { return re; }
  #pragma hila loop_function
  inline T& real() { return re; }

  #pragma hila loop_function
  inline T imag() const { return im; }
  #pragma hila loop_function
  inline T& imag() { return im; }


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
  inline T norm_sq() const { return re*re + im*im; }

  // TODO: make this work for vector type!  Not double  
  //currently this gives a compilation error
  #pragma hila loop_function
  inline T abs() const { return sqrt( norm_sq() ); }
  #pragma hila loop_function
  inline T arg() const { return atan2( im, re ); }


  #pragma hila loop_function
  inline cmplx<T> conj() const { return cmplx<T>( { re, -im } ); }

  #pragma hila loop_function
  inline cmplx<T> polar(const T r, const T theta) { 
    return cmplx<T>( { r*cos(theta), r*sin(theta) } );
  }

  #pragma hila loop_function
  inline cmplx<T> & random(){
    re = hila_random();
    im = hila_random();
    return *this;
  }

  #pragma hila loop_function
  inline cmplx<T> & gaussian(){
    re = gaussian_ran2(im);
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
    T n = lhs.norm_sq();
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

// functions real(), imag()

template <typename T>
#pragma hila loop_function
inline T real(const cmplx<T> a) { return a.re; }

template <typename T>
#pragma hila loop_function
inline T imag(const cmplx<T> a) { return a.im; }


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
//   T n = b.norm_sq();
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
  return val.norm_sq();
}

template <typename T>
#pragma hila loop_function
inline void random( cmplx<T> & c ) {
  random(c.re);
  random(c.im); 
}

template <typename T>
#pragma hila loop_function
inline void gaussian_random( cmplx<T> & c ) {
  gaussian_random(c.re);
  gaussian_random(c.im); 
}

//////////////////////////////////////////////////////////////////////////

// define also real(), imag(), conj() -functions for basic arithmetic types
template<typename T, std::enable_if_t<is_arithmetic<T>::value,int> = 0 >
#pragma hila loop_function
inline T real(T val){
  return val;
}

template<typename T, std::enable_if_t<is_arithmetic<T>::value,int> = 0 >
#pragma hila loop_function
inline T imag(T val){
  return static_cast<T>(0);
}

template<typename T, std::enable_if_t<is_arithmetic<T>::value,int> = 0 >
#pragma hila loop_function
inline T conj(T val){
  return val;
}



////////////////////////////////////////////////////////////////////////
/// And utility templates
/// Define is_cmplx<T>::value -template, using specialization
template< typename T>
struct is_cmplx : std::integral_constant<
  bool, false
> {};

template< typename T>
struct is_cmplx<cmplx<T>> : std::integral_constant<
  bool, true
> {};

// and a template is_cmplx_or_real<T>::value
template< typename T >
struct is_cmplx_or_arithmetic : std::integral_constant<
  bool,
  is_arithmetic<T>::value || is_cmplx<T>::value
> {};



#endif

