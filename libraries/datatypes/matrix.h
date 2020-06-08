#ifndef MATRIX_H
#define MATRIX_H

#include<type_traits>

#include "general_matrix.h"


template <int n, typename T>
struct squarematrix {
  using base_type = typename base_type_struct<T>::type;
  T c[n][n];

  squarematrix() = default;

  //constructor from general matrix  
  #pragma hila loop_function
  squarematrix(const matrix<n,n,T> rhs) {
    for (int i=0; i<n; i++){
      for (int j=0; j<n; j++) {
        c[i][j] = rhs.c[i][j];
      }
    }
  }

  //constructor from scalar  
  template <typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 >  
  #pragma hila loop_function
  squarematrix(const scalart rhs) {
    for (int i=0; i<n; i++){
      for (int j=0; j<n; j++) {
        c[i][j] = static_cast<T>(0);
      }
      c[i][i] = static_cast<T>(rhs);
    }
  }

  template <typename scalart,
            std::enable_if_t<std::is_arithmetic<scalart>::value, int> = 0 >  
  #pragma hila loop_function
  squarematrix<n,T> & operator= (const scalart rhs) {
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        c[i][j] = static_cast<T>(0);
      }
      c[i][i] = static_cast<T>(rhs);
    }
    return *this;
  }
  
  squarematrix<n,T> transpose() {
    squarematrix<n,T> res;
    for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
      res.c[i][j] = c[j][i];
    }
    return res;
  }

  squarematrix<n,T> conjugate() const {
    squarematrix<n,T> res;
    for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
      res.c[i][j] = conj(c[j][i]);
    }
    return res;
  }


  squarematrix<n,T> operator-(){
    squarematrix<n,T> r;
    for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
      r.c[i][j] = -c[j][i];
    }
    return *r;
  }


  //*=, +=, -= operators
  #pragma hila loop_function
  squarematrix<n,T> & operator+=(const squarematrix<n,T> & rhs){
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++){
      c[i][j] += rhs.c[i][j]; 
    }
    return *this;
  }

  #pragma hila loop_function
  squarematrix<n,T> & operator-=(const squarematrix<n,T> & rhs){
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++){
      c[i][j] -= rhs.c[i][j]; 
    }
    return *this;
  }

  template <typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 >
  #pragma hila loop_function
  squarematrix<n,T> & operator*=(const scalart rhs){
    T val;
    val=rhs;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++){
      c[i][j]*=val;
    }
    return *this;
  }


  #pragma hila loop_function
  squarematrix<n,T> & operator*=(const squarematrix<n,T> & rhs){
    squarematrix<n,T> res;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++){
      res.c[i][j] = (0);
      for (int k = 0; k < n; k++){
        res.c[i][j] += (c[i][k] * rhs.c[k][j]);
      }
    }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++){
      c[i][j] = res.c[i][j];
    }
    return *this;
  }

  T trace() const {
    T result = static_cast<T>(0);
    for (int i = 0; i < n; i++){
      result += c[i][i];
    }
    return result;
  }

  template <typename A=T, std::enable_if_t<is_arithmetic<A>::value, int> = 0 > 
  squarematrix<n,T> & random(){
    for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
      c[i][j] = static_cast<T>(hila_random());
    }
    return *this;
  }

  #pragma hila loop_function
  template <typename A=T, std::enable_if_t<!is_arithmetic<A>::value, int> = 0 > 
  squarematrix<n,T> & random(){
    for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
      c[i][j].random();
    }
    return *this;
  }

  auto norm_sq(){
    auto result = norm_squared(c[0][0]);
    for (int i=0; i<n; i++) for (int j=0; j<n; j++) if(i>0&&j>0) {
      result += norm_squared(c[i][j]);
    }
    return result;
  }

  // find determinant using LU decomposition. Algorithm: numerical Recipes, 2nd ed. p. 47 ff
  template <typename radix, std::enable_if_t<is_arithmetic<radix>::value, int> = 0, std::enable_if_t<std::is_same<T,cmplx<radix>>::value, int> = 0 > 
  cmplx<radix> det_lu(){
    int i, imax, j, k;
    radix big, d, temp, dum;
    cmplx<radix> cdum, csum, ctmp1;
    radix vv[n];
    cmplx<radix> a[n][n];
    cmplx<radix> one;

    one=cmplx<radix>(1,0);
    d=1;
    imax = -1;

    for (i=0; i<n; i++) for(j=0; j<n; j++) a[i][j] = this->c[i][j];
    for (i=0; i<n; i++) {
      big = 0;
      for(j=0; j<n; j++) {
        if ((temp = a[i][j].abs()) > big)
          big = temp;
      }
      asser(big != 0.0 && "Determinant does not exist\n");
      vv[i] = 1.0/big;
    }

    for (j=0; j<n; j++) {
      for (i=0; i<j; i++) {
        csum = a[i][j];
        for (k=0; k<i; k++) {
          csum -= a[i][k]*a[k][j];
        }
        a[i][j] = csum;
      }

      big = 0;
      for (i=j; i<n; i++) {
        csum = a[i][j];
        for (k=0; k<j; k++) {
            csum -= a[i][k]*a[k][j];
        }
        a[i][j] = csum;
        if ((dum = vv[i]*csum.abs()) >= big) {
            big = dum;
            imax = i;
        }
      }

      if (j != imax) {
        for (k=0; k<n; k++) {
          cdum = a[imax][k];
          a[imax][k] = a[j][k];
          a[j][k] = cdum;
        }
        d = -d;
        vv[imax] = vv[j];
      }

      if (a[j][j].abs() == static_cast<radix>(0.0)) 
        a[j][j] = cmplx<radix>(1e-20,0);

      if (j != n-1) {
        cdum = one/a[j][j];
        for (i=j+1; i<n; i++) {
          a[i][j] = a[i][j]*cdum;
        }
      }
    }

    csum = cmplx<radix>(d,0.0);
    for (j=0; j<n; j++) {
      csum = csum*a[j][j];
    }

    return (csum);
  }
  
  std::string str() const {
    std::string text = "";
    for (int i=0; i<n; i++){
      for (int j=0; j<n; j++) {
        text + c[i][j].str() + " "; 
      }
      text + "\n"; 
    }
    return text;
  }
};




template <int n, typename T1, typename T2, typename Tr=decltype(std::declval<T1>() * std::declval<T2>())>
squarematrix<n,Tr>
operator* (const squarematrix<n,T1> &A, const squarematrix<n,T2> &B)
{
  squarematrix<n,Tr> res;
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) res.c[i][j] = static_cast<Tr>(0);
    for (int k=0; k<n; k++) for (int j=0; j<n; j++) {
      res.c[i][j] +=  A.c[i][k] * B.c[k][j];
    }
  }
  return res;
}

template <int n, typename T1, typename T2, typename Tr=decltype(std::declval<T1>() + std::declval<T2>())>
squarematrix<n,Tr>
operator+ (const squarematrix<n,T1> &A, const squarematrix<n,T2> &B)
{
  squarematrix<n,Tr> res;
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
    res.c[i][j] =  A.c[i][j] + B.c[i][j];
  }
  return res;
}

template <int n, typename T>
squarematrix<n,T> operator- (const squarematrix<n,T> &A, const squarematrix<n,T> &B) {
  squarematrix<n,T> res;
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
    res.c[i][j] =  A.c[i][j] - B.c[i][j];
  }
  return res;
}


// multiplication by a scalar
template <int n, typename T, typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 > 
#pragma hila loop_function
squarematrix<n,T> operator* (const squarematrix<n,T> &A, const scalart s) {
  squarematrix<n,T> res;
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
    res.c[i][j] = s * A.c[i][j];
  }
  return res;
}

template <int n, typename T, typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 > 
#pragma hila loop_function
squarematrix<n,T> operator/ (const squarematrix<n,T> &A, const scalart s) {
  squarematrix<n,T> res;
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
    res.c[i][j] = s / A.c[i][j];
  }
  return res;
}

template <int n, typename T, typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 > 
#pragma hila loop_function
squarematrix<n,T> operator*(const scalart s, const squarematrix<n,T> &A) {
  return operator*(A,s);
}


template <int n, typename T, std::enable_if_t<!is_arithmetic<T>::value, int> = 0 >
#pragma hila loop_function
squarematrix<n,T> operator * (const squarematrix<n,T> & A, const T & B) {
  squarematrix<n,T> res;
  for (int i = 0; i < n; i++) for (int j = 0; j < n; j++){
    res.c[i][j] = A.c[i][j] * B;
  }
  return res;
}

template <int n, typename T, std::enable_if_t<!is_arithmetic<T>::value, int> = 0 >
#pragma hila loop_function
squarematrix<n,T> operator * (const T & A, const squarematrix<n,T> & B) {
  squarematrix<n,T> res;
  for (int i = 0; i < n; i++) for (int j = 0; j < n; j++){
    res.c[i][j] = B * A.c[i][j];
  }
  return res;
}


template <int n, typename T> 
#pragma hila loop_function
squarematrix<n,T> operator ^ (const squarematrix<n,T> & A, const int pow) {
  squarematrix<n,T> res;
  res = 1;
  for (int i = 0; i < pow; i++){
    res *= A;
  }
  return res;
}


template <int n, typename T> 
#pragma hila loop_function
std::ostream& operator<<(std::ostream &strm, const squarematrix<n,T> &A) {
  for (int i=0; i<n; i++){
    strm << "\n"; 
    for (int j=0; j<n; j++) {
      strm << " " << A.c[i][j] << " "; 
    }
    strm << "\n"; 
  }
  strm << "\n"; 
  return strm;
}


template <int n, typename T> 
inline auto norm_squared(squarematrix<n,T> & rhs){
  auto result = norm_squared(rhs.c[0][0]);
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) if(i>0&&j>0) {
    result += norm_squared(rhs.c[i][j]);
  }
  return result;
}



#endif
