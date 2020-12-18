#ifndef GENERAL_MATRIX_H_
#define GENERAL_MATRIX_H_
#include<type_traits>
#include "operations.h"
#include "datatypes/cmplx.h"


// Do this macro here to ease "switching" between using
// mul_add operation and the normal sum
// THE mul_add METHOD SEEMS TO BE SLOWER?  
#define MUL_SUM(a, b, c) c += a*b
// #define MUL_SUM(a, b, c) c = mul_add(a, b, c)

// do here transpose() and adjoint() functions for numbers and cmplxes

template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0 >
inline T transpose(const T v) { return v; }

template <typename T>
inline cmplx<T> transpose(const cmplx<T> v) { return v; }

template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0 >
inline T adjoint(const T v) { return v; }

template <typename T>
inline cmplx<T> adjoint(const cmplx<T> v) { return v.conj(); }




template <const int n, const int m, typename T>
class matrix {
  private:
    T c[n*m];

  public:
    // std incantation for field types
    using base_type = typename base_type_struct<T>::type;

    // define these to ensure std::is_trivial
    matrix() = default;
    ~matrix() = default,
    matrix(const matrix<n,m,T> & v) = default;

    // standard access ops m.e(i,j) - assume T is small, as it usually is
    #pragma hila loop_function
    inline T  e(const int i, const int j) const { return c[i*m + j]; }
    #pragma hila loop_function
    inline T& e(const int i, const int j) { return c[i*m + j]; }

    // unary -
    matrix<n,m,T> operator-() const { return -1*(*this); }

    // Assign from "scalar" for square matrix
    template <typename S, std::enable_if_t<is_assignable<T&,S>::value, int> = 0 >  
    #pragma hila loop_function
    matrix<n,m,T> & operator= (const S rhs) {
      static_assert( n==m, "rows != columns : assigning a scalar only possible for a square matrix");
      for (int i=0; i<n*m; i++) c[i] = zero;
      for (int i=0; i<n; i++) e(i,i) = rhs;
      return *this;
    }

    //copy constructor from scalar
    template <typename S, std::enable_if_t<is_assignable<T&,S>::value, int> = 0 >  
    #pragma hila loop_function
    matrix(const S rhs) {
      static_assert(n==m, "rows != columns : scalar assignment possible for square matrices only!");
      for (int i=0; i<n*m; i++) c[i] = zero;
      for (int i=0; i<n; i++) e(i,i) = rhs;
    }

    // assign and construct from zero
    #pragma hila loop_function
    matrix(const Zero z) {
      for (int i=0; i<n*m; i++) c[i] = zero;
    }

    #pragma hila loop_function
    matrix<n,m,T> & operator= (const Zero z) {
      for (int i=0; i<n*m; i++) c[i] = zero;
      return *this;
    }
   

    //*=, +=, -= operators
    #pragma hila loop_function
    template <typename S, 
              std::enable_if_t<std::is_convertible<type_plus<T,S>,T>::value, int> = 0 >
    matrix<n,m,T> & operator+=(const matrix<n,m,S> & rhs){
      for (int i = 0; i < n; i++) for (int j = 0; j < m; j++){
        e(i,j) += rhs.e(i,j);
      }
      return *this;
    }

    #pragma hila loop_function
    template <typename S,
              std::enable_if_t<std::is_convertible<type_minus<T,S>,T>::value, int> = 0 >
    matrix<n,m,T> & operator-=(const matrix<n,m,T> & rhs){
      for (int i = 0; i < n; i++) for (int j = 0; j < m; j++){
        e(i,j) -= rhs.e(i,j);
      }
      return *this;
    }

    #pragma hila loop_function
    template <typename S,
              std::enable_if_t<std::is_convertible<type_mul<T,S>,T>::value, int> = 0 >
    matrix<n,m,T> & operator*=(const S rhs){
      T val;
      val=rhs;
      for (int i=0; i < n*m; i++) {
        c[i] *= val;
      }
      return *this;
    }

    #pragma hila loop_function
    template <int p, typename S,
              std::enable_if_t<std::is_convertible<type_mul<T,S>,T>::value, int> = 0 >
    matrix<n,m,T> & operator*=(const matrix<m,p,S> & rhs){
      static_assert(m==p, "can't assign result of *= to lhs matrix, because doing so would change it's dimensions");
      *this = *this * res;
      return *this;
    }

    //numpy style matrix fill 
    template <typename S, std::enable_if_t<is_assignable<T&,S>::value, int> = 0 >  
    #pragma hila loop_function
    matrix<n,m,T> & fill(const S rhs) {
      T t = rhs;
      for (int i=0; i<n*m; i++) c[i] = t;
      return *this;
    }
  
    //return copy of transpose of this matrix
    #pragma hila loop_function
    matrix<m,n,T> transpose() const { return transpose(*this); }

    //return copy of complex conjugate = adjoint of this matrix
    #pragma hila loop_function
    matrix<m,n,T> adjoint() const { return adjoint(*this); }

    #pragma hila loop_function
    T trace() const {
      static_assert(n==m, "trace not defined for non square matrices!");
      T result = e(0,0);
      for (int i = 1; i < n; i++){
        result += e(i,i);
      }
      return result;
    }

    #pragma hila loop_function
    template <typename A=T, std::enable_if_t<is_arithmetic<A>::value, int> = 0 > 
    matrix<n, m, A> & random(){
    for (int i=0; i<n; i++) for (int j=0; j<m; j++) {
      c[i][j] = static_cast<T>(hila_random());
    }
    return *this;
  }

  #pragma hila loop_function
  template <typename A=T, std::enable_if_t<!is_arithmetic<A>::value, int> = 0 > 
  matrix<n, m, A> & random(){
    for (int i=0; i<n; i++) for (int j=0; j<m; j++) {
      c[i][j].random();
    }
    return *this;
  }

  auto norm_sq(){
    auto result = norm_squared(c[0][0]);
    for (int i=0; i<n; i++) for (int j=0; j<m; j++) if(i>0&&j>0) {
      result += norm_squared(c[i][j]);
    }
    return result;
  }

  inline T dot(const matrix<n, m, T> &rhs) const {
    T r = (0.0);
    for (int i=0; i<n; i++) for (int j=0; j<m; j++) {
      r += conj(c[i][j])*rhs.c[i][j];
    }
    return r;
  }

  std::string str() const {
    std::string text = "";
    for (int i=0; i<n; i++){
      for (int j=0; j<m; j++) {
        text + c[i][j].str() + " "; 
      }
      text + "\n"; 
    }
    return text;
  }
};


// do transpose and adjoint functions here
template <const int n, const int m, typename T> 
inline matrix<n,m,T> transpose(const matrix<m,n,T> & rhs) {
  matrix<n,m,T> res;
  for (int i=0; i<n; i++) for (int j=0; j<m, j++) {
    res.e(i,j) = transpose(rhs.e(j,i));
  }
  return res;
}
// and adjoint function
template <const int n, const int m, typename T> 
inline matrix<n,m,T> adjoint(const matrix<m,n,T> & rhs) {
  matrix<n,m,T> res;
  for (int i=0; i<n; i++) for (int j=0; j<m, j++) {
    res.e(i,j) = adjoint(rhs.e(j,i));
  }
  return res;
}

//templates needed for naive calculation of determinants

template<int n, int m, typename T> 
#pragma hila loop_function
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

//determinant -> use LU factorization later 
template<int n, int m, typename T> 
#pragma hila loop_function
T det(const matrix<n, m, T> & mat){
  static_assert(n==m, "determinants defined only for square matrices");
  T result = 0.0;
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
#pragma hila loop_function
T det(const matrix<2,2,T> & mat){
  return mat.c[0][0]*mat.c[1][1] - mat.c[1][0]*mat.c[0][1];
}

//matrix multiplication for 2 by 2 matrices ; 
template<typename T> 
#pragma hila loop_function
matrix<2,2,T> operator* (const matrix<2,2,T> &A, const matrix<2,2,T> &B) {
  matrix<2,2,T> res = 1;
  res.c[0][0] = A.c[0][0]*B.c[0][0] + A.c[0][1]*B.c[1][0];
  res.c[0][1] = A.c[0][0]*B.c[0][1] + A.c[0][1]*B.c[1][1];
  res.c[1][1] = A.c[1][0]*B.c[0][1] + A.c[1][1]*B.c[1][1];
  res.c[1][0] = A.c[1][0]*B.c[0][0] + A.c[1][1]*B.c[1][0];
  return res;
}

//matrix power 
template <int n, int m, typename T> 
#pragma hila loop_function
matrix<n,m,T> operator ^ (const matrix<n,m,T> & A, const int pow) {
  matrix<n,m,T> res;
  res = 1;
  for (int i = 0; i < pow; i++){
    res *= A;
  }
  return res;
}

//matrix * scalar 
template <int n, int m, typename T> 
#pragma hila loop_function
matrix<n,m,T> operator * (const matrix<n,m,T> & A, const T & B) {
  matrix<n,m,T> res;
  for (int i = 0; i < n; i++) for (int j = 0; j < m; j++){
    res.c[i][j] = A.c[i][j] * B;
  }
  return res;
}



//general matrix * matrix multiplication 
template <int n, int m, int p, typename T> 
#pragma hila loop_function
matrix<n,p,T> operator * (const matrix<n,m,T> &A, const matrix<m,p,T> &B) {
  matrix<n,p,T> res;
  for (int i = 0; i < n; i++) for (int j = 0; j < p; j++){
    res.c[i][j] = 0;
    for (int k = 0; k < m; k++){
      // res.c[i][j] += (A.c[i][k] * B.c[k][j]);
      MUL_SUM( A.c[i][k] , B.c[k][j] , res.c[i][j] );
    }
  }
  return res;
}

//multiplication for matrix * transpose matrix  
template <int n, int m, int p, typename T> 
#pragma hila loop_function
matrix<n,p,T> operator * (const matrix<n,m,T> & A, const transposeMatrix<p,m,T> & B) {
  matrix<n,p,T> res;
  for (int i = 0; i < n; i++) for (int j = 0; j < p; j++){
    res.c[i][j] = 0;
    for (int k = 0; k < m; k++){
      // res.c[i][j] += (A.c[i][k]*B.ref.c[j][k]);
      MUL_SUM(A.c[i][k] , B.ref.c[j][k], res.c[i][j]);
    }
  }
  return res;
}

//multiplication for transpose * matrix  
template <int n, int m, int p, typename T> 
#pragma hila loop_function
matrix<n,p,T> operator * (const transposeMatrix<m,n,T> & A, const matrix<m,p,T> & B) {
  matrix<n,p,T> res;
  for (int i = 0; i < n; i++) for (int j = 0; j < p; j++){
    res.c[i][j] = 0;
    for (int k = 0; k < m; k++){
      //res.c[i][j] += (A.ref.c[k][i]*B.c[k][j]);
      MUL_SUM(A.ref.c[k][i] , B.c[k][j], res.c[i][j]);
    }
  }
  return res;
}


//transpose * transpose
template <int n, int m, int p, typename T> 
#pragma hila loop_function
matrix<n,p,T> operator * (const transposeMatrix<m,n,T> & A, const transposeMatrix<p,m,T> & B) {
  matrix<n,p,T> res;
  for (int i = 0; i < n; i++) for (int j = 0; j < p; j++){
    res.c[i][j] = 0;
    for (int k = 0; k < m; k++){
      // res.c[i][j] += (A.ref.c[k][i]*B.ref.c[j][k]);
      MUL_SUM(A.ref.c[k][i] , B.ref.c[j][k], res.c[i][j]);
    }
  }
  return res;
}

//multiplication for matrix * conjugate matrix  
template <int n, int m, int p, typename T> 
#pragma hila loop_function
matrix<n,p,T> operator * (const matrix<n,m,T> & A, const conjugateMatrix<p,m,T> & B) {
  matrix<n,p,T> res;
  for (int i = 0; i < n; i++) for (int j = 0; j < p; j++){
    res.c[i][j] = 0;
    for (int k = 0; k < m; k++){
      // res.c[i][j] += (A.c[i][k]*conj(B.ref.c[j][k]));
      MUL_SUM( A.c[i][k] , conj(B.ref.c[j][k]), res.c[i][j] );
    }
  }
  return res;
}

//multiplication for conjugate * matrix  
template <int n, int m, int p, typename T> 
#pragma hila loop_function
matrix<n,p,T> operator * (const conjugateMatrix<m,n,T> & A, const matrix<m,p,T> & B) {
  matrix<n,p,T> res;
  for (int i = 0; i < n; i++) for (int j = 0; j < p; j++){
    res.c[i][j] = 0;
    for (int k = 0; k < m; k++){
      // res.c[i][j] += (conj(A.ref.c[k][i])*B.c[k][j]);
      MUL_SUM( conj(A.ref.c[k][i]) , B.c[k][j], res.c[i][j] );
    }
  }
  return res;
}


//conjugate * conjugate
template <int n, int m, int p, typename T> 
#pragma hila loop_function
matrix<n,m,T> operator * (const conjugateMatrix<m,n,T> & A, const conjugateMatrix<p,m,T> & B) {
  matrix<n,p,T> res;
  for (int i = 0; i < n; i++) for (int j = 0; j < p; j++){
    res.c[i][j] = 0;
    for (int k = 0; k < m; k++){
      res.c[i][j] += (conj(A.ref.c[k][i])*conj(B.ref.c[j][k]));
    }
  }
  return res;
}

//addition for matrix + conjugate matrix  
template <int n, int m, typename T> 
#pragma hila loop_function
matrix<n,m,T> operator + (const matrix<n,m,T> & A, const conjugateMatrix<m,n,T> & B) {
  matrix<n,m,T> res;
  for (int i = 0; i < n; i++) for (int j = 0; j < m; j++){
    res.c[i][j] = A.c[i][j] + conj(B.ref.c[j][i]);
  }
  return res;
}

//addition for conjugate + matrix  
template <int n, int m, typename T> 
#pragma hila loop_function
matrix<n,m,T> operator + (const conjugateMatrix<m,n,T> & A, const matrix<n,m,T> & B) {
  matrix<n,m,T> res;
  for (int i = 0; i < n; i++) for (int j = 0; j < m; j++){
    res.c[i][j] = conj(A.ref.c[j][i]) + B.c[i][j];
  }
  return res;
}

//addition for conjugate + conjugate
template <int n, int m, typename T> 
#pragma hila loop_function
matrix<n,m,T> operator + (const conjugateMatrix<m,n,T> & A, const conjugateMatrix<n,m,T> & B) {
  matrix<n,m,T> res;
  for (int i = 0; i < n; i++) for (int j = 0; j < m; j++){
      res.c[i][j] = conj(A.ref.c[j][i]) + conj(B.ref.c[j][i]);
  }
  return res;
}

//subtraction for matrix - conjugate matrix  
template <int n, int m, typename T> 
#pragma hila loop_function
matrix<n,m,T> operator - (const matrix<n,m,T> & A, const conjugateMatrix<m,n,T> & B) {
  matrix<n,m,T> res;
  for (int i = 0; i < n; i++) for (int j = 0; j < m; j++){
    res.c[i][j] = A.c[i][j] - conj(B.ref.c[j][i]);
  }
  return res;
}

//subtraction for conjugate - matrix  
template <int n, int m, typename T> 
#pragma hila loop_function
matrix<n,m,T> operator - (const conjugateMatrix<m,n,T> & A, const matrix<n,m,T> & B) {
  matrix<n,m,T> res;
  for (int i = 0; i < n; i++) for (int j = 0; j < m; j++){
    res.c[i][j] = conj(A.ref.c[j][i]) - B.c[i][j];
  }
  return res;
}

//subtraction for conjugate - conjugate
template <int n, int m, typename T> 
#pragma hila loop_function
matrix<n,m,T> operator - (const conjugateMatrix<m,n,T> & A, const conjugateMatrix<n,m,T> & B) {
  matrix<n,m,T> res;
  for (int i = 0; i < n; i++) for (int j = 0; j < m; j++){
      res.c[i][j] = conj(A.ref.c[j][i]) - conj(B.ref.c[j][i]);
  }
  return res;
}


//component wise addition
template <int n, int m, typename T> 
#pragma hila loop_function
matrix<n,m,T> operator+ (const matrix<n,m,T> &A, const matrix<n,m,T> &B) {
  matrix<n,m,T> res;
  for (int i=0; i<n; i++) for (int j=0; j<m; j++) {
    res.c[i][j] =  A.c[i][j] + B.c[i][j];
  }
  return res;
}

//component wise subtraction
template <int n, int m, typename T> 
#pragma hila loop_function
matrix<n,m,T> operator- (const matrix<n,m,T> &A, const matrix<n,m,T> &B) {
  matrix<n,m,T> res;
  for (int i=0; i<n; i++) for (int j=0; j<m; j++) {
    res.c[i][j] =  A.c[i][j] - B.c[i][j];
  }
  return res;
}

// multiplication by a scalar
template <int n, int m, typename T, typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 > 
#pragma hila loop_function
matrix<n,m,T> operator* (const matrix<n,m,T> &A, const scalart s) {
  matrix<n,m,T> res;
  for (int i=0; i<n; i++) for (int j=0; j<m; j++) {
    res.c[i][j] = s * A.c[i][j];
  }
  return res;
}

template <int n, int m, typename T, typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 > 
#pragma hila loop_function
matrix<n,m,T> operator/ (const matrix<n,m,T> &A, const scalart s) {
  matrix<n,m,T> res;
  for (int i=0; i<n; i++) for (int j=0; j<m; j++) {
    res.c[i][j] = s / A.c[i][j];
  }
  return res;
}

template <int n, int m, typename T, typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 > 
#pragma hila loop_function
matrix<n,m,T> operator*(const scalart s, const matrix<n,m,T> &A) {
  return operator*(A,s);
}


template <int n, int m, typename T>
#pragma hila loop_function
std::ostream& operator<<(std::ostream &strm, const matrix<n,m,T> &A) {
  for (int i=0; i<n; i++){
    strm << "\n"; 
    for (int j=0; j<m; j++) {
      strm << " " << A.c[i][j] << " "; 
    }
    strm << "\n"; 
  }
  strm << "\n";
  return strm;
}

template<int n, int m, typename T> 
inline transposeMatrix<n,m,T> trans(matrix<n,m,T> & ref){
  transposeMatrix<n,m,T> result(ref);
  return result;
}

template<int n, int m, typename T> 
inline conjugateMatrix<n,m,T> conj(matrix<n,m,T> & ref){
  conjugateMatrix<n,m,T> result(ref);
  return result;
}

template<int n, int m, typename T>
inline auto norm_squared(matrix<n,m,T> & rhs){
  auto result = norm_squared(rhs.c[0][0]);
  for (int i=0; i<n; i++) for (int j=0; j<m; j++) if(i>0&&j>0) {
    result += norm_squared(rhs.c[i][j]);
  }
  return result;
}

#endif
