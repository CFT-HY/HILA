#ifndef CMPLX_H
#define CMPLX_H

#include <complex>


template <typename T = double>
struct cmplx {
  T re,im;
  
  // assignment is automatically OK, by c-standard
  //   cmplx operator=(cmplx rhs) { 
  //     re = rhs.re; im = rhs.im; 
  //     return *this; 
  //   }

  cmplx<T> & operator=(std::complex<int> & c) {
    re = c.real(); im = c.imag();
    return *this;
  }
  cmplx<T> & operator=(std::complex<double> & c) {
    re = c.real(); im = c.imag();
    return *this;
  }

  
  template <typename scalar_t,
            std::enable_if_t<std::is_arithmetic<scalar_t>::value, int> = 0 >
  cmplx<T> & operator=(scalar_t s) {
    re = s;
    im = 0.0;
    return *this;
  }
  
  T real() { return re; }
  T imag() { return im; }

  T norm() { return re*re + im*im; }
  double abs()  { return sqrt((double)norm()); }

  cmplx<T> conj() { return cmplx<T>( { re, -im } ); }

  cmplx<T> polar(const T r, const T theta);
  
};

template <typename T>
cmplx<T> operator+(const cmplx<T> & a, const cmplx<T> & b) {
  return cmplx<T>( { a.re + b.re, a.im + b.im } );
}


// Operators to implement  3 + 2i  etc.
/*
cmplx<long double> operator "" i(long double d) {
  return cmplx<long double>(0.0,d);
}
cmplx<long double> operator "" i(unsigned long long int d) {
  return cmplx<long double>(0.0,d);
}
*/



#endif
