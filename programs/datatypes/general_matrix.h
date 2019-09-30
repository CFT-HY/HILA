#ifndef MATRIX_H
#define MATRIX_H
#include<type_traits>
#include "cmplx.h"

//--- conjugation for general type T -> template specialized for complex types

template<typename T>
inline T typeConj(T val){
  return val; 
}

template<typename Accuracy>
inline cmplx<Accuracy> typeConj(cmplx<Accuracy> val){
  return val.conj();
}

//---

template <int n, int m, typename T>
class matrix {
  public:
  T c[n][m];

  matrix() = default;

  template <typename scalart, std::enable_if_t<std::is_arithmetic<scalart>::value, int> = 0 >  
  matrix<n,m,T> & operator= (const scalart rhs) {
    static_assert(n==m, "rowdim != coldim : cannot assign diagonal from scalar!");
    for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
      if (i == j) c[i][j] = (rhs);
      else c[i][j] = (0);
    }
    return *this;
  }

  //copy constructor from scalar  
  template <typename scalart, std::enable_if_t<std::is_arithmetic<scalart>::value, int> = 0 >  
  matrix(const scalart rhs) {
    static_assert(n==m, "rowdim != coldim : cannot assign diagonal from scalar!");
    for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
      if (i == j) c[i][j] = (rhs);
      else c[i][j] = (0);
    }
  }

  //*=, +=, -= operators
  matrix<n,m,T> & operator+=(const matrix<n,m,T> & rhs){
    for (int i = 0; i < n; i++) for (int j = 0; j < m; j++){
      c[i][j] += rhs.c[i][j]; 
    }
    return *this;
  }

  matrix<n,m,T> & operator-=(const matrix<n,m,T> & rhs){
    for (int i = 0; i < n; i++) for (int j = 0; j < m; j++){
      c[i][j] -= rhs.c[i][j]; 
    }
    return *this;
  }

  template <typename scalart, std::enable_if_t<std::is_arithmetic<scalart>::value, int> = 0 >
  matrix<n,m,T> & operator*=(const scalart rhs){
    T val;
    val=rhs;
    for (int i = 0; i < n; i++) for (int j = 0; j < m; j++){
      c[i][j]*=val;
    }
    return *this;
  }

  template<int p>
  matrix<n,m,T> & operator*=(const matrix<m,p,T> & rhs){
    static_assert(m==p, "can't assign result of *= to matrix A, because doing so would change it's dimensions");
    matrix<m,m,T> rhsTrans = rhs.transpose();
    matrix<n,m,T> res;
    for (int i = 0; i < n; i++) for (int j = 0; j < m; j++){
      res.c[i][j] = (0);
      for (int k = 0; k < m; k++){
        res.c[i][j] += (c[i][k] * rhsTrans.c[j][k]);
      }
    }
    for (int i = 0; i < n; i++) for (int j = 0; j < m; j++){
      c[i][j] = res.c[i][j];
    }
    return *this;
  }

  //numpy style matrix fill 
  template <typename scalart, std::enable_if_t<std::is_arithmetic<scalart>::value, int> = 0 > 
  matrix<n,m,T> & fill(const scalart rhs) {
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++){
      c[i][j] = (rhs);
    }
    return *this;
  }
  
  //return copy of transpose of this matrix
  matrix<m,n,T> transpose() const {
    matrix<m,n,T> res;
    for (int i=0; i<m; i++) for (int j=0; j<n; j++) {
      res.c[i][j] =  c[j][i];
    }
    return res;
  }

  //return copy of complex conjugate of this matrix
  matrix<m,n,T> conjugate() const {
    matrix<m,n,T> res;
    for (int i=0; i<m; i++) for (int j=0; j<n; j++) {
      res.c[i][j] =  typeConj(c[j][i]);
    }
    return res;
  }

  T trace() const {
    static_assert(n==m, "trace not defined for non square matrices!");
    T result = static_cast<T>(0);
    for (int i = 0; i < n; i++){
      result += c[i][i];
    }
    return result;
  }
};

//templates needed for naive calculation of determinants

template<int n, int m, typename T>
matrix<n - 1, m - 1, T> Minor(const matrix<n, m, T> & bigger, int i, int j){
  matrix<n - 1, m - 1, T> result;
  int index = 0;
  for (int p = 0; p < n; p++) for (int l = 0; l < m; l++){
    if (p==i || l==j) continue;
    *(*(result.c) + index) = bigger.c[p][l];
    index++;
  }
  return result;
}

template<int n, int m, typename T>
T det(const matrix<n, m, T> & mat){
  static_assert(n==m, "determinants defined only for square matrices");
  T result = 1.0;
  T parity = 1.0; //assumes that copy constructor from scalar has been defined for T 
  T opposite = -1.0; 
  for (int i = 0; i < n; i++){
    matrix<n - 1, m - 1, T> minor = Minor(mat, 0, i);
    result += parity*det(minor)*mat.c[0][i];
    parity*=opposite;
  }
  return result;
}

template<typename T>
T det(const matrix<2,2,T> & mat){
  return mat.c[0][0]*mat.c[1][1] - mat.c[1][0]*mat.c[0][1];
}

//matrix multiplication for 2 by 2 matrices ; 
template<typename T>
matrix<2,2,T> operator* (const matrix<2,2,T> &A, const matrix<2,2,T> &B) {
  matrix<2,2,T> res = 1;
  res.c[0][0] = A.c[0][0]*B.c[0][0] + A.c[0][1]*B.c[1][0];
  res.c[0][1] = A.c[0][0]*B.c[0][1] + A.c[0][1]*B.c[1][1];
  res.c[1][1] = A.c[1][0]*B.c[0][1] + A.c[1][1]*B.c[1][1];
  res.c[1][0] = A.c[1][0]*B.c[0][0] + A.c[1][1]*B.c[1][0];
  return res;
}

//general naive matrix multiplication 
template <int n, int m, int p, typename T>
matrix<n,p,T> operator* (const matrix<n,m,T> &A, const matrix<m,p,T> &B) {
  matrix<p,m,T> Btrans = B.transpose(); //do matrix multiplication on rows of transpose matrix (should reduce cache misses)
  matrix<n,p,T> res;
  for (int i = 0; i < n; i++) for (int j = 0; j < p; j++){
    res.c[i][j] = (0);
    for (int k = 0; k < m; k++){
      res.c[i][j] += (A.c[i][k] * Btrans.c[j][k]);
    }
  }
  return res;
}

//dot product definitions for vectors
template<int n, typename T>
T operator* (const matrix<1, n, T> & vecA, const matrix<n, 1, T> & vecB) {
  T result = (0.0);
  for (int i = 0; i < n; i++){
    result += vecA.c[0][i]*(typeConj(vecB.c[i][0]));
  }
  return result;
}

template<int n, typename T>
T operator* (const matrix<1, n, T> & vecA, const matrix<1, n, T> & vecB) {
  T result = (0.0);
  for (int i = 0; i < n; i++){
    result += vecA.c[0][i]*(typeConj(vecB.c[0][i]));
  }
  return result;
}

template<int n, typename T>
T operator* (const matrix<n, 1, T> & vecA, const matrix<n, 1, T> & vecB) {
  T result; 
  result=0;
  for (int i = 0; i < n; i++){
    result += vecA.c[i][0]*(typeConj(vecB.c[i][0]));
  }
  return result;
}

//component wise addition
template <int n, int m, typename T>
matrix<n,m,T> operator+ (const matrix<n,m,T> &A, const matrix<n,m,T> &B) {
  matrix<n,m,T> res;
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
    res.c[i][j] =  A.c[i][j] + B.c[i][j];
  }
  return res;
}

//component wise subtraction
template <int n, int m, typename T>
matrix<n,m,T> operator- (const matrix<n,m,T> &A, const matrix<n,m,T> &B) {
  matrix<n,m,T> res;
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
    res.c[i][j] =  A.c[i][j] - B.c[i][j];
  }
  return res;
}

// multiplication by a scalar
template <int n, int m, typename T, typename scalart, std::enable_if_t<std::is_arithmetic<scalart>::value, int> = 0 >
matrix<n,m,T> operator* (const matrix<n,m,T> &A, const scalart s) {
  matrix<n,m,T> res;
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
    res.c[i][j] = s * A.c[i][j];
  }
  return res;
}

template <int n, int m, typename T, typename scalart, std::enable_if_t<std::is_arithmetic<scalart>::value, int> = 0 >
matrix<n,m,T> operator/ (const matrix<n,m,T> &A, const scalart s) {
  matrix<n,m,T> res;
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
    res.c[i][j] = s / A.c[i][j];
  }
  return res;
}

template <int n, int m, typename T, typename scalart, std::enable_if_t<std::is_arithmetic<scalart>::value, int> = 0 >
matrix<n,m,T> operator*(const scalart s, const matrix<n,m,T> &A) {
  return operator*(A,s);
}

#endif
