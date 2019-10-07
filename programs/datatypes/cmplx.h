#ifndef CMPLX_H
#define CMPLX_H

// let's not include the std::complex
//#include <complex>
//#include <cmath>


template <typename T = double>
struct cmplx {
  T re,im;
  
  // assignment is automatically OK, by c-standard
  //   cmplx operator=(cmplx rhs) { 
  //     re = rhs.re; im = rhs.im; 
  //     return *this; 
  //   }
  cmplx<T>() = default;
  
  cmplx<T>(const cmplx<T> & a) =default;

  // constructor from single scalar value 
  template <typename scalar_t,
            std::enable_if_t<std::is_arithmetic<scalar_t>::value, int> = 0 >
  constexpr cmplx<T>(const scalar_t val): re(static_cast<T>(val)), im(static_cast<T>(0)) {}
  
  // constructor c(a,b)
//   template <typename A, typename B,
//             std::enable_if_t<std::is_arithmetic<A>::value, int> = 0,
//             std::enable_if_t<std::is_arithmetic<B>::value, int> = 0 >
//   constexpr cmplx<T>(const A & a, const B & b) {
//     re = static_cast<T>(a);
//     im = static_cast<T>(b);
//   }

  // constructor c(a,b)
  template <typename A, typename B,
            std::enable_if_t<std::is_arithmetic<A>::value, int> = 0,
            std::enable_if_t<std::is_arithmetic<B>::value, int> = 0 >
  constexpr cmplx<T>(const A & a, const B & b): re(static_cast<T>(a)), im(static_cast<T>(b)) {}

  
  ~cmplx<T>() =default;
  
  // automatic casting from cmplx<T> -> cmplx<A>
  // TODO: ensure this works if A is vector type!
  template <typename A>
  operator cmplx<A>() const { 
    return cmplx<A>( { static_cast<A>(re), static_cast<A>(im) });
  }

  constexpr static int base_element_count(){
    if constexpr( std::is_arithmetic<T>::value ) {
      return 2;
    } else {
      return 2*T::base_element_count();
    }
  }

  constexpr static int base_element_size(){
    if constexpr( std::is_arithmetic<T>::value ) {
      return sizeof(T);
    } else {
      return T::base_element_size();
    }
  }

  void set_from_pointers(char ** pointers){
    if constexpr( std::is_arithmetic<T>::value ) {
      re = *((T*) pointers[0]);
      im = *((T*) pointers[1]);
    } else {
      re.set_from_pointers(pointers);
      im.set_from_pointers(pointers + T::base_element_count());
    }
  }

  void set_to_pointers(char ** pointers){
    if constexpr( std::is_arithmetic<T>::value ) {
      *((T *) pointers[0]) = re;
      *((T *) pointers[1]) = im;
    } else {
      re.set_to_pointers(pointers);
      im.set_to_pointers(pointers + T::base_element_count());
    }
  }


//   // assignment from std::complex<A>  TODO: perhaps remove?
//   template <typename A>
//   cmplx<T> & operator=(const std::complex<A> & c) {
//     re = static_cast<T>(c.real()); im = static_cast<T>(c.imag());
//     return *this;
//   }
  
  template <typename scalar_t,
            std::enable_if_t<std::is_arithmetic<scalar_t>::value, int> = 0 >
  cmplx<T> & operator=(scalar_t s) {
    re = static_cast<T>(s);
    im = 0.0;
    return *this;
  }
  
  T real() const { return re; }
  T imag() const { return im; }

  T norm() const { return re*re + im*im; }
  // TODO: make this work for vector type!  Not double
  
  //currently this gives a compilation error
  double abs() const { return sqrt(static_cast<double>(norm()) ); }
  double arg() const { return atan2(static_cast<double>(im),static_cast<double>(re)); }


  cmplx<T> conj() const { return cmplx<T>( { re, -im } ); }

  cmplx<T> polar(const T r, const T theta) { 
    return cmplx<T>( { r*cos(theta), r*sin(theta) } );
  }

  // unary + and -
  cmplx<T> operator+() const {return *this;}
  cmplx<T> operator-() const {return cmplx<T>(-re, -im); }

  
  cmplx<T> & operator+= (const cmplx<T> & lhs) {
    re += lhs.re;
    im += lhs.im;
    return *this;
  }

  // TODO: for avx vector too -- #define new template macro
  template <typename A,
            std::enable_if_t<std::is_arithmetic<A>::value, int> = 0 >
  cmplx<T> & operator+= (const A & a) {
    re += static_cast<T>(a);
    return *this;
  }
  
  cmplx<T> & operator-= (const cmplx<T> & lhs) {
    re -= lhs.re;
    im -= lhs.im;
    return *this;
  }
  
  // TODO: for vector too
  template <typename A,
            std::enable_if_t<std::is_arithmetic<A>::value, int> = 0 >
  cmplx<T> & operator-= (const A & a) {
    re -= static_cast<T>(a);
    return *this;
  }
  
  cmplx<T> & operator*= (const cmplx<T> & lhs) {
    re = re * lhs.re - im * lhs.im;
    im = im * lhs.re + re * lhs.im;
    return *this;
  }
  
  // TODO: for vector too
  template <typename A,
            std::enable_if_t<std::is_arithmetic<A>::value, int> = 0 >
  cmplx<T> & operator*= (const A & a) {
    re *= static_cast<T>(a);
    im *= static_cast<T>(a);
    return *this;
  }

  // a/b = a b*/|b|^2 = (a.re*b.re + a.im*b.im + i(a.im*b.re - a.re*b.im))/|b|^2
  cmplx<T> & operator/= (const cmplx<T> & lhs) {
    T n = lhs.norm();
    re = (re * lhs.re + im * lhs.im)/n;
    im = (im * lhs.re - re * lhs.im)/n;
    return *this;
  }
  
  // TODO: for vector too
  template <typename A,
            std::enable_if_t<std::is_arithmetic<A>::value, int> = 0 >
  cmplx<T> & operator/= (const A & a) {
    re /= static_cast<T>(a);
    im /= static_cast<T>(a);
    return *this;
  }

  
};

template <typename T>
cmplx<T> operator+(const cmplx<T> & a, const cmplx<T> & b) {
  return cmplx<T>(a.re + b.re, a.im + b.im);
}

  // TODO: for avx vector too -- #define new template macro
template <typename T, typename A,
          std::enable_if_t<std::is_arithmetic<A>::value, int> = 0 >
cmplx<T> operator+(const cmplx<T> & c, const A & a) {
  return cmplx<T>(c.re + a, c.im);
}

template <typename T, typename A,
          std::enable_if_t<std::is_arithmetic<A>::value, int> = 0 >
cmplx<T> operator+(const A &a, const cmplx<T> & c) {
  return cmplx<T>(c.re + a, c.im);
}

// -
template <typename T>
cmplx<T> operator-(const cmplx<T> & a, const cmplx<T> & b) {
  return cmplx<T>(a.re - b.re, a.im - b.im);
}

// TODO: for avx vector too -- #define new template macro
template <typename T, typename A,
          std::enable_if_t<std::is_arithmetic<A>::value, int> = 0 >
cmplx<T> operator-(const cmplx<T> & c, const A & a) {
  return cmplx<T>(c.re - a, c.im);
}

template <typename T, typename A,
          std::enable_if_t<std::is_arithmetic<A>::value, int> = 0 >
cmplx<T> operator-(const A &a, const cmplx<T> & c) {
  return cmplx<T>(a - c.re, -c.im);
}


// *
template <typename T>
cmplx<T> operator*(const cmplx<T> & a, const cmplx<T> & b) {
  return cmplx<T>(a.re*b.re - a.im*b.im, a.im*b.re + a.re*b.im);
}

// TODO: for avx vector too -- #define new template macro
template <typename T, typename A,
          std::enable_if_t<std::is_arithmetic<A>::value, int> = 0 >
cmplx<T> operator*(const cmplx<T> & c, const A & a) {
  return cmplx<T>(c.re * a, c.im * a);
}

template <typename T, typename A,
          std::enable_if_t<std::is_arithmetic<A>::value, int> = 0 >
cmplx<T> operator*(const A &a, const cmplx<T> & c) {
  return cmplx<T>(a * c.re, a * c.im);
}


// /   a/b = ab*/|b|^2
template <typename T>
cmplx<T> operator/(const cmplx<T> & a, const cmplx<T> & b) {
  T n = b.norm();
  return cmplx<T>( (a.re*b.re + a.im*b.im)/n, (a.im*b.re - a.re*b.im)/n );
}

// TODO: for avx vector too -- #define new template macro
template <typename T, typename A,
          std::enable_if_t<std::is_arithmetic<A>::value, int> = 0 >
cmplx<T> operator/(const cmplx<T> & c, const A & a) {
  return cmplx<T>(c.re / a, c.im / a);
}

// a/c = ac*/|c|^2 
template <typename T, typename A,
          std::enable_if_t<std::is_arithmetic<A>::value, int> = 0 >
cmplx<T> operator/(const A &a, const cmplx<T> & c) {
  T n = c.norm();
  return cmplx<T>((a * c.re)/n, -(a * c.im)/n);
}


// Operators to implement imaginary unit 1_i, enablig expressions  3 + 2_i  etc.
// Underscore seems to be required here

constexpr cmplx<double> operator""_i(long double a) {
  return cmplx<double>{0.0,a};
}

constexpr cmplx<double> operator""_i(unsigned long long a) {
  return cmplx<double>(0.0,static_cast<double>(a));
}


#endif
