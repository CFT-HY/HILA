#ifndef MATRIX_H_
#define MATRIX_H_
#include<type_traits>
#include "operations.h"
#include "datatypes/cmplx.h"


// Do this macro here to ease "switching" between using
// mul_add operation and the normal sum
// THE mul_add METHOD SEEMS TO BE SLOWER?
#define MUL_SUM(a, b, c) c += a*b
// #define MUL_SUM(a, b, c) c = mul_add(a, b, c)

// do here transpose() and adjoint() functions for numbers and cmplxes

template <typename T, std::enable_if_t<is_arithmetic<T>::value, int> = 0 >
inline T transpose(const T v) { return v; }

template <typename T>
inline cmplx<T> transpose(const cmplx<T> v) { return v; }

template <typename T, std::enable_if_t<is_arithmetic<T>::value, int> = 0 >
inline T adjoint(const T v) { return v; }

template <typename T>
inline cmplx<T> adjoint(const cmplx<T> v) { return v.conj(); }

// these fwd declarations seem to be needed
template <const int n, const int m, typename T>
class Matrix;
template <const int n, const int m, typename T>
Matrix<n,m,T> transpose(const Matrix<m,n,T> & rhs);
template <const int n, const int m, typename T>
Matrix<n,m,T> adjoint(const Matrix<m,n,T> & rhs);




template <const int n, const int m, typename T>
class Matrix {
  private:
    T c[n*m];

  public:
    // std incantation for field types
    using base_type = typename base_type_struct<T>::type;

    // define these to ensure std::is_trivial
    Matrix() = default;
    ~Matrix() = default;
    Matrix(const Matrix<n,m,T> & v) = default;

    // standard access ops m.e(i,j) - assume T is small, as it usually is
    #pragma hila loop_function
    inline T  e(const int i, const int j) const { return c[i*m + j]; }
    #pragma hila loop_function
    inline T& e(const int i, const int j) { return c[i*m + j]; }

    // declare single e here too in case we have a vector
    // (one size == 1)
    #pragma hila loop_function
    template <int q=n, int p=m,
              std::enable_if_t< (q == 1 || p == 1), int> = 0 >
    inline T e(const int i) const { return c[i]; }

    #pragma hila loop_function
    template <int q=n, int p=m,
              std::enable_if_t< (q == 1 || p == 1), int> = 0 >
    inline T& e(const int i) { return c[i]; }


    // casting from one Matrix (number) type to another: do not do this automatically.
    // but require an explicit cast operator.  This makes it easier to write code.
    // or should it be automatic?  keep/remove explicit?
    // TODO: CHECK AVX CONVERSIONS

    template <typename S,
              std::enable_if_t<std::is_convertible<T,S>::value, int> = 0 >
    explicit operator Matrix<n,m,S>() {
      Matrix<n,m,S> res;
      for (int i=0; i<n*m; i++) res.c[i] = static_cast<S>(c[i]);
      return res;
    }


    // unary -
    Matrix<n,m,T> operator-() const { return -1*(*this); }

    // Assign from "scalar" for square matrix
    template <typename S, std::enable_if_t<is_assignable<T&,S>::value, int> = 0 >
    #pragma hila loop_function
    Matrix<n,m,T> & operator= (const S rhs) {
      static_assert( n==m, "rows != columns : assigning a scalar only possible for a square Matrix");
      for (int i=0; i<n*m; i++) c[i] = zero;
      for (int i=0; i<n; i++) e(i,i) = rhs;
      return *this;
    }

    //copy constructor from scalar
    template <typename S, std::enable_if_t<is_assignable<T&,S>::value, int> = 0 >
    #pragma hila loop_function
    Matrix(const S rhs) {
      static_assert(n==m, "rows != columns : scalar assignment possible for square matrices only!");
      for (int i=0; i<n*m; i++) c[i] = zero;
      for (int i=0; i<n; i++) e(i,i) = rhs;
    }

    // assign and construct from zero
    #pragma hila loop_function
    Matrix(const Zero z) {
      for (int i=0; i<n*m; i++) c[i] = zero;
    }

    #pragma hila loop_function
    Matrix<n,m,T> & operator= (const Zero z) {
      for (int i=0; i<n*m; i++) c[i] = zero;
      return *this;
    }


    // +=, -= operators for Matrix and scalar
    #pragma hila loop_function
    Matrix<n,m,T> & operator+=(const Matrix<n,m,T> & rhs){
      for (int i=0; i < n*m; i++) {
        c[i] += rhs.c[i];
      }
      return *this;
    }

    #pragma hila loop_function
    Matrix<n,m,T> & operator-=(const Matrix<n,m,T> & rhs){
      for (int i=0; i < n*m; i++) {
        c[i] -= rhs.c[i];
      }
      return *this;
    }

    #pragma hila loop_function
    template <typename S,
              std::enable_if_t<std::is_convertible<type_plus<T,S>,T>::value, int> = 0 >
    Matrix<n,m,T> & operator+=(const S rhs){
      static_assert(n==m, "rows != columns : scalar addition possible for square matrices only!");
      for (int i=0; i < n; i++) {
        e(i,i) += rhs;
      }
      return *this;
    }

    #pragma hila loop_function
    template <typename S,
              std::enable_if_t<std::is_convertible<type_minus<T,S>,T>::value, int> = 0 >
    Matrix<n,m,T> & operator-=(const S rhs){
      static_assert(n==m, "rows != columns : scalar subtraction possible for square matrices only!");
      for (int i=0; i < n; i++) {
        e(i,i) -= rhs;
      }
      return *this;
    }

    // *= operator
    #pragma hila loop_function
    template <int p>
    Matrix<n,m,T> & operator*=(const Matrix<m,p,T> & rhs){
      static_assert(m==p, "can't assign result of *= to lhs Matrix, because doing so would change it's dimensions");
      *this = *this * rhs;
      return *this;
    }

    // *= for scalar
    #pragma hila loop_function
    template <typename S,
              std::enable_if_t<std::is_convertible<type_mul<T,S>,T>::value, int> = 0 >
    Matrix<n,m,T> & operator*=(const S rhs) {
      for (int i=0; i<n*m; i++) {
        c[i] *= rhs;
      }
      return *this;
    }

    // and operator /= for a scalar rhs
    #pragma hila loop_function
    template <typename S,
              std::enable_if_t<std::is_convertible<type_div<T,S>,T>::value, int> = 0 >
    Matrix<n,m,T> & operator/=(const S rhs) {
      for (int i=0; i<n*m; i++) {
        c[i] /= rhs;
      }
      return *this;
    }

    //numpy style matrix fill
    template <typename S, std::enable_if_t<is_assignable<T&,S>::value, int> = 0 >
    #pragma hila loop_function
    Matrix<n,m,T> & fill(const S rhs) {
      T t = static_cast<T>(rhs);
      for (int i=0; i<n*m; i++) c[i] = t;
      return *this;
    }

    //return copy of transpose of this matrix
    #pragma hila loop_function
    Matrix<m,n,T> transpose() const { return ::transpose(*this); }

    //return copy of complex conjugate = adjoint of this Matrix
    #pragma hila loop_function
    Matrix<m,n,T> adjoint() const { return ::adjoint(*this); }

    #pragma hila loop_function
    T trace() const {
      static_assert(n==m, "trace not defined for non square matrices!");
      T result = e(0,0);
      for (int i = 1; i < n; i++){
        result += e(i,i);
      }
      return result;
    }


    auto norm_sq() const {
      auto result = norm_squared(c[0]);
      for (int i=1; i<n*m; i++) {
        result += norm_squared(c[i]);
      }
      return result;
    }


    // enable dot only for vectors (confusing otherwise)
    auto dot(const Matrix<n,m,T> &rhs) const {
      static_assert( n==1 || m==1, "dot product only for vectors");
      auto r = adjoint(c[0]) * rhs.c[0];
      for (int i=1; i<n*m; i++) {
        r += adjoint(c[i]) * rhs.c[i];
      }
      return r;
    }

    #pragma hila loop_function
    Matrix<n, m, T> & random() {
      for (int i=0; i<n*m; i++) {
        ::random(c[i]);
      }
      return *this;
    }

    #pragma hila loop_function
    inline Matrix<n, m, T> & gaussian(){ 
      for (int i = 0; i < n*m; i++) {
        ::gaussian_random(c[i]);
      }
      return *this;
    }



    std::string str() const {
      std::stringstream text;
      for (int i=0; i<n; i++){
        for (int j=0; j<m; j++) {
          text << e(i,j) << " ";
        }
        text << '\n';
      }
      return text.str();
    }
};


// do transpose and adjoint functions here
template <const int n, const int m, typename T>
inline Matrix<n,m,T> transpose(const Matrix<m,n,T> & rhs) {
  Matrix<n,m,T> res;
  for (int i=0; i<n; i++) for (int j=0; j<m; j++) {
    res.e(i,j) = transpose(rhs.e(j,i));
  }
  return res;
}
// and adjoint function
template <const int n, const int m, typename T>
inline Matrix<n,m,T> adjoint(const Matrix<m,n,T> & rhs) {
  Matrix<n,m,T> res;
  for (int i=0; i<n; i++) for (int j=0; j<m; j++) {
    res.e(i,j) = adjoint(rhs.e(j,i));
  }
  return res;
}

//templates needed for naive calculation of determinants

template<int n, int m, typename T>
#pragma hila loop_function
Matrix<n - 1, m - 1, T> Minor(const Matrix<n, m, T> & bigger, int i, int j){
  Matrix<n - 1, m - 1, T> result;
  int ri = 0, bi = 0;
  for (int p = 0; p < n; p++) for (int l = 0; l < m; l++, bi++){
    if (p != i && l != j) {
      result.c[ri] = bigger.c[bi];
      ri++;
    }
  }
  return result;
}

//determinant -> use LU factorization later
template<int n, int m, typename T>
#pragma hila loop_function
T det(const Matrix<n, m, T> & mat){
  static_assert(n==m, "determinants defined only for square matrices");
  T result = 0.0;  // I'll assume this works!
  number_type<T> parity = 1, opposite = -1;
  for (int i = 0; i < n; i++){
    Matrix<n - 1, m - 1, T> minor = Minor(mat, 0, i);
    result += parity*det(minor)*mat.e(0,i);
    parity *= opposite;
  }
  return result;
}

// and 2x2 and 1x1 matrix dets, too...
template<typename T>
#pragma hila loop_function
T det(const Matrix<2,2,T> & mat){
  return mat.e(0,0)*mat.e(1,1) - mat.e(1,0)*mat.e(0,1);
}

template<typename T>
#pragma hila loop_function
T det(const Matrix<1,1,T> & mat) {
  return mat.e(0,0);
}


//Now matrix additions: matrix + matrix

#pragma hila loop_function
template <int n, int m, typename T>
inline Matrix<n,m,T> operator+(Matrix<n,m,T> a, const Matrix<n,m,T> & b){
  a += b;
  return a;
}

#pragma hila loop_function
template <int n, int m, typename T>
inline Matrix<n,m,T> operator-(Matrix<n,m,T> a, const Matrix<n,m,T> & b){
  a -= b;
  return a;
}

// Matrix + scalar
#pragma hila loop_function
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_plus<T,S>,T>::value, int> = 0 >
inline Matrix<n,m,T> operator+(Matrix<n,m,T> a, const S b){
  a += b;
  return a;
}

// scalar + matrix
#pragma hila loop_function
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_plus<T,S>,T>::value, int> = 0 >
inline Matrix<n,m,T> operator+(const S b, Matrix<n,m,T> a){
  a += b;
  return a;
}

// matrix - scalar
#pragma hila loop_function
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_minus<T,S>,T>::value, int> = 0 >
Matrix<n,m,T> operator-(Matrix<n,m,T> a, const S b){
  a -= b;
  return a;
}

// scalar - matrix
#pragma hila loop_function
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_minus<S,T>,T>::value, int> = 0 >
inline Matrix<n,m,T> operator-(const S b, Matrix<n,m,T> a){
  static_assert(n==m, "rows != columns : scalar subtraction possible for square matrices only!");
  for (int i=0; i<n; i++) a.e(i,i) = static_cast<T>(b) - a.e(i,i);
  return a;
}

////////////////////////////////////////
/// matrix * matrix is the crucial bit here

template <int n, int m, int p, typename T>
#pragma hila loop_function
Matrix<n,p,T> operator*(const Matrix<n,m,T> & A, const Matrix<m,p,T> & B) {
  Matrix<n,p,T> res(zero);

  if constexpr (n > 1 && p > 1) {
    // normal matrix*matrix
    for (int i=0; i<n; i++) for (int j=0; j<p; j++) {
      // res.e(i,j) = zero;
      for (int k=0; k<m; k++) {
        res.e(i,j) += A.e(i,k)*B.e(k,j);
      }
    }
  } else if constexpr ( p == 1 ) {
    // matrix * vector
    for (int i=0; i<n; i++) {
      // res.e(i) = zero;
      for (int k=0; k<m; k++) {
        res.e(i) += A.e(i,k)*B.e(k);
      }
    }
  } else if constexpr ( n == 1 ) {
    // horiz. vector * matrix
    for (int j=0; j<p; j++) {
      //res.e(j) = zero;
      for (int k=0; k<m; k++) {
        res.e(j) += A.e(k)*B.e(k,j);
      }
    }
  }
  return res;
}


// matrix * scalar
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_mul<T,S>,T>::value, int> = 0 >
Matrix<n,m,T> operator*(Matrix<n,m,T> mat, const S rhs) {
  mat *= rhs;
  return mat;
}

// scalar * matrix
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_mul<S,T>,T>::value, int> = 0 >
Matrix<n,m,T> operator*(const S rhs, Matrix<n,m,T> mat) {
  mat *= rhs;         // assume here commutativity!  Fails if the "scalar" is also a matrix
  return mat;
}

// matrix / scalar
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_div<T,S>,T>::value, int> = 0 >
Matrix<n,m,T> operator/(Matrix<n,m,T> mat, const S rhs) {
  mat /= rhs;
  return mat;
}


template <int n, int m, typename T>
std::ostream& operator<<(std::ostream &strm, const Matrix<n,m,T> &A) {
  for (int i=0; i<n; i++){
    strm << "\n";
    for (int j=0; j<m; j++) {
      strm << " " << A.e(i,j) << " ";
    }
    strm << "\n";
  }
  strm << "\n";
  return strm;
}


template<int n, int m, typename T>
inline auto norm_squared(Matrix<n,m,T> & rhs){
  auto result = norm_squared(rhs.c[0]);
  for (int i=0; i<n*m; i++) {
    result += norm_squared(rhs.c[i]);
  }
  return result;
}

template<int n, int m, typename T>
inline void random(Matrix<n,m,T> & mat) {
  mat.random();
}

template<int n, int m, typename T>
inline void gaussian_random(Matrix<n,m,T> & mat) {
  mat.gaussian();
}


// find determinant using LU decomposition. Algorithm: numerical Recipes, 2nd ed. p. 47 ff
template <int n, typename T, typename radix=number_type<T>,
          std::enable_if_t<is_arithmetic<radix>::value, int> = 0,
          std::enable_if_t<std::is_same<T,cmplx<radix>>::value, int> = 0 >
cmplx<radix> det_lu( const Matrix<n,n,T> & mat) {

  int i, imax, j, k;
  radix big, d, temp, dum;
  cmplx<radix> cdum, csum, ctmp1;
  radix vv[n];
  cmplx<radix> a[n][n];
  cmplx<radix> one;

  one = cmplx<radix>(1,0);
  d = 1;
  imax = -1;

  for (i=0; i<n; i++) for(j=0; j<n; j++) a[i][j] = mat.e(i,j);
  for (i=0; i<n; i++) {
    big = 0;
    for(j=0; j<n; j++) {
      if ((temp = a[i][j].abs()) > big)
        big = temp;
    }
    assert(big != 0.0 && "Determinant does not exist\n");
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


template <int n, typename T, std::enable_if_t<is_arithmetic<T>::value, int> = 0>
T det_lu(const Matrix<n,n,T> & mat) {
  int i, imax, j, k;
  T big, d, temp, dum;
  T cdum, csum, ctmp1;
  T vv[n];
  T a[n][n];
  d=1;
  imax = -1;

  for (i=0; i<n; i++) for(j=0; j<n; j++) a[i][j] = mat.e(i,j);
  for (i=0; i<n; i++) {
    big = 0;
    for(j=0; j<n; j++) {
      if ((temp = a[i][j]) > big)
        big = temp;
    }
    assert(big != 0.0 && "Determinant does not exist\n");
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
      if ((dum = vv[i]*csum) >= big) {
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

    if (a[j][j] == static_cast<T>(0.0))
      a[j][j] = 1e-20;

    if (j != n-1) {
      cdum = 1.0/a[j][j];
      for (i=j+1; i<n; i++) {
        a[i][j] = a[i][j]*cdum;
      }
    }
  }

  csum = d;
  for (j=0; j<n; j++) {
    csum = csum*a[j][j];
  }

  return (csum);
}

///////////////////////////////////////////////////////////////
/// Finally, define Vector, HorizontalVector and SquareMatrix
///////////////////////////////////////////////////////////////

template <int n,typename T>
using Vector = Matrix<n,1,T>;

template <int n,typename T>
using HorizontalVector = Matrix<1,n,T>;

template <int n,typename T>
using SquareMatrix = Matrix<n,n,T>;








#endif
