#ifndef MATRIX_H
#define MATRIX_H

#include<type_traits>

#include "cmplx.h"


template <int n, typename T = real_t>
struct cmatrix {
  cmplx<T> c[n][n];

  // do nothing constructor
  // cmatrix<n,T>() = default;
  
  //   // copy construct
  //   cmatrix<n,T>(const cmatrix<n,T> &other) {
  //     *this = other;
  //     return *this;
  //   }
  // 
  //   // another copy construct
  //   template <typename scalart,
  //             std::enable_if_t<std::is_arithmetic<scalart>::value, int> = 0 >    
  //   cmatrix<n,T>(const scalart s) {
  //     this->operator=(s);
  //   }
  
  // cmatrix = cmatrix assignment should happen automatically

  template <typename scalart,
            std::enable_if_t<std::is_arithmetic<scalart>::value, int> = 0 >  
  loop_callable cmatrix<n,T> & operator= (const scalart rhs) {
    for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
      if (i == j) c[i][j] = static_cast<T>(rhs);
      else c[i][j] = static_cast<T>(0);
    }
    return *this;
  }
  
  cmatrix<n,T> & transpose() {
    cmatrix<n,T> res;
    for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
      res.c[i][j] =  c[j][i];
    }
    return res;
  }
};



template <int n, typename T>
loop_callable cmatrix<n,T> operator* (const cmatrix<n,T> &A, const cmatrix<n,T> &B) {
  cmatrix<n,T> res;
  // not sure if this order is the best, but at least the j-loop
  // is in contiquous memory
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) res.c[i][j] = static_cast<T>(0);
    for (int k=0; k<n; k++) for (int j=0; j<n; j++) {
      res.c[i][j] +=  A.c[i][k] * B.c[k][j];
    }
  }
  return res;
}

template <int n, typename T>
loop_callable cmatrix<n,T> operator+ (const cmatrix<n,T> &A, const cmatrix<n,T> &B) {
  cmatrix<n,T> res;
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
    res.c[i][j] =  A.c[i][j] + B.c[i][j];
  }
  return res;
}

// multiply by a scalar 
template <int n, typename T, typename scalart,
          std::enable_if_t<std::is_arithmetic<scalart>::value, int> = 0 >
loop_callable cmatrix<n,T> operator*(const cmatrix<n,T> &A, const scalart s) {
  cmatrix<n,T> res;
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
    res.c[i][j] = s * A.c[i][j];
  }
  return res;
}

template <int n, typename T, typename scalart,
          std::enable_if_t<std::is_arithmetic<scalart>::value, int> = 0 >
loop_callable cmatrix<n,T> operator*(const scalart s, const cmatrix<n,T> &A) {
  return operator*(A,s);
}


#endif
