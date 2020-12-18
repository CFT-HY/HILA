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
class matrix;
template <const int n, const int m, typename T>
matrix<n,m,T> transpose(const matrix<m,n,T> & rhs);
template <const int n, const int m, typename T>
matrix<n,m,T> adjoint(const matrix<m,n,T> & rhs);




template <const int n, const int m, typename T>
class matrix {
  private:
    T c[n*m];

  public:
    // std incantation for field types
    using base_type = typename base_type_struct<T>::type;

    // define these to ensure std::is_trivial
    matrix() = default;
    ~matrix() = default;
    matrix(const matrix<n,m,T> & v) = default;

    // standard access ops m.e(i,j) - assume T is small, as it usually is
    #pragma hila loop_function
    inline T  e(const int i, const int j) const { return c[i*m + j]; }
    #pragma hila loop_function
    inline T& e(const int i, const int j) { return c[i*m + j]; }

    // casting from one matrix (number) type to another: do not do this automatically.
    // but require an explicit cast operator.  This makes it easier to write code.
    // or should it be automatic?  keep/remove explicit?
    // TODO: CHECK AVX CONVERSIONS

    template <typename S, 
              std::enable_if_t<std::is_convertible<T,S>::value, int> = 0 >
    explicit operator matrix<n,m,S>() {
      matrix<n,m,S> res;
      for (int i=0; i<n*m; i++) res.c[i] = static_cast<S>(c[i]);
      return res;
    }


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
   

    // +=, -= operators for matrix and scalar
    #pragma hila loop_function
    matrix<n,m,T> & operator+=(const matrix<n,m,T> & rhs){
      for (int i=0; i < n*m; i++) {
        c[i] += rhs.c[i];
      }
      return *this;
    }

    #pragma hila loop_function
    matrix<n,m,T> & operator-=(const matrix<n,m,T> & rhs){
      for (int i=0; i < n*m; i++) {
        c[i] -= rhs.c[i];
      }
      return *this;
    }

    #pragma hila loop_function
    template <typename S,
              std::enable_if_t<std::is_convertible<type_plus<T,S>,T>::value, int> = 0 >
    matrix<n,m,T> & operator+=(const S rhs){
      static_assert(n==m, "rows != columns : scalar addition possible for square matrices only!");
      for (int i=0; i < n; i++) {
        e(i,i) += rhs;
      }
      return *this;
    }

    #pragma hila loop_function
    template <typename S,
              std::enable_if_t<std::is_convertible<type_minus<T,S>,T>::value, int> = 0 >
    matrix<n,m,T> & operator-=(const S rhs){
      static_assert(n==m, "rows != columns : scalar subtraction possible for square matrices only!");
      for (int i=0; i < n; i++) {
        e(i,i) -= rhs;
      }
      return *this;
    }

    // *= operator
    #pragma hila loop_function
    template <int p>
    matrix<n,m,T> & operator*=(const matrix<m,p,T> & rhs){
      static_assert(m==p, "can't assign result of *= to lhs matrix, because doing so would change it's dimensions");
      *this = *this * rhs;
      return *this;
    }

    // *= for scalar
    #pragma hila loop_function
    template <typename S,
              std::enable_if_t<std::is_convertible<type_mul<T,S>,T>::value, int> = 0 >
    matrix<n,m,T> & operator*=(const S rhs) {
      for (int i=0; i<n*m; i++) {
        c[i] *= rhs;
      }
      return *this;
    }

    // and operator /= for a scalar rhs
    #pragma hila loop_function
    template <typename S,
              std::enable_if_t<std::is_convertible<type_div<T,S>,T>::value, int> = 0 >
    matrix<n,m,T> & operator/=(const S rhs) {
      for (int i=0; i<n*m; i++) {
        c[i] /= rhs;
      }
      return *this;
    }

    //numpy style matrix fill 
    template <typename S, std::enable_if_t<is_assignable<T&,S>::value, int> = 0 >  
    #pragma hila loop_function
    matrix<n,m,T> & fill(const S rhs) {
      T t = static_cast<T>(rhs);
      for (int i=0; i<n*m; i++) c[i] = t;
      return *this;
    }
  
    //return copy of transpose of this matrix
    #pragma hila loop_function
    matrix<m,n,T> transpose() const { return ::transpose(*this); }

    //return copy of complex conjugate = adjoint of this matrix
    #pragma hila loop_function
    matrix<m,n,T> adjoint() const { return ::adjoint(*this); }

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
        e(i,j) = static_cast<T>(hila_random());
      }
      return *this; 
    }

    #pragma hila loop_function
    template <typename A=T, std::enable_if_t<!is_arithmetic<A>::value, int> = 0 > 
    matrix<n, m, A> & random() {
      for (int i=0; i<n*m; i++) {
        c[i].random();
      }
      return *this;
    }

    auto norm_sq() const {
      auto result = norm_squared(c[0]);
      for (int i=1; i<n*m; i++) {
        result += norm_squared(c[i]);
      }
      return result;
    }

    inline T dot(const matrix<n,m,T> &rhs) const {
      T r = adjoint(c[0]) * rhs.c[0];
      for (int i=1; i<n*m; i++) {
        r += adjoint(c[i]) * rhs.c[i];
      }
      return r;
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
inline matrix<n,m,T> transpose(const matrix<m,n,T> & rhs) {
  matrix<n,m,T> res;
  for (int i=0; i<n; i++) for (int j=0; j<m; j++) {
    res.e(i,j) = transpose(rhs.e(j,i));
  }
  return res;
}
// and adjoint function
template <const int n, const int m, typename T> 
inline matrix<n,m,T> adjoint(const matrix<m,n,T> & rhs) {
  matrix<n,m,T> res;
  for (int i=0; i<n; i++) for (int j=0; j<m; j++) {
    res.e(i,j) = adjoint(rhs.e(j,i));
  }
  return res;
}

//templates needed for naive calculation of determinants

template<int n, int m, typename T> 
#pragma hila loop_function
matrix<n - 1, m - 1, T> Minor(const matrix<n, m, T> & bigger, int i, int j){
  matrix<n - 1, m - 1, T> result;
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
T det(const matrix<n, m, T> & mat){
  static_assert(n==m, "determinants defined only for square matrices");
  T result = 0.0;  // I'll assume this works!
  number_type<T> parity = 1, opposite = -1;
  for (int i = 0; i < n; i++){
    matrix<n - 1, m - 1, T> minor = Minor(mat, 0, i);
    result += parity*det(minor)*mat.e(0,i);
    parity *= opposite;
  }
  return result;
}

// and 2x2 and 1x1 matrix dets, too...
template<typename T> 
#pragma hila loop_function
T det(const matrix<2,2,T> & mat){
  return mat.e(0,0)*mat.e(1,1) - mat.e(1,0)*mat.e(0,1);
}

template<typename T> 
#pragma hila loop_function
T det(const matrix<1,1,T> & mat) {
  return mat.e(0,0);
}


//Now matrix additions: matrix + matrix

#pragma hila loop_function
template <int n, int m, typename T> 
inline matrix<n,m,T> operator+(matrix<n,m,T> a, const matrix<n,m,T> & b){
  a += b;
  return a;
}

#pragma hila loop_function
template <int n, int m, typename T> 
inline matrix<n,m,T> operator-(matrix<n,m,T> a, const matrix<n,m,T> & b){
  a -= b;
  return a;
}

// matrix + scalar
#pragma hila loop_function
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_plus<T,S>,T>::value, int> = 0 >
inline matrix<n,m,T> operator+(matrix<n,m,T> a, const S b){
  a += b;
  return a;
}

// scalar + matrix
#pragma hila loop_function
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_plus<T,S>,T>::value, int> = 0 >
inline matrix<n,m,T> operator+(const S b, matrix<n,m,T> a){
  a += b;
  return a;
}

// matrix - scalar
#pragma hila loop_function
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_minus<T,S>,T>::value, int> = 0 >
matrix<n,m,T> operator-(matrix<n,m,T> a, const S b){
  a -= b;
  return a;
}

// scalar - matrix
#pragma hila loop_function
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_minus<S,T>,T>::value, int> = 0 >
inline matrix<n,m,T> operator-(const S b, matrix<n,m,T> a){
  static_assert(n==m, "rows != columns : scalar subtraction possible for square matrices only!");
  for (int i=0; i<n; i++) a.e(i,i) = static_cast<T>(b) - a.e(i,i);
  return a;
}

//////////
/// matrix * matrix

template <int n, int m, int p, typename T> 
#pragma hila loop_function
matrix<n,p,T> operator*(const matrix<n,m,T> & A, const matrix<m,p,T> & B) {
  matrix<n,p,T> res;
  for (int i=0; i<n; i++) {
    for (int j=0; j<p; j++)   res.e(i,j)  = A.e(i,0)*B.e(0,j);
    for (int k=1; k<m; k++) {
      for (int j=0; j<p; j++) res.e(i,j) += A.e(i,k)*B.e(k,j);
    } 
  }
  return res;
}


// matrix * scalar
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_mul<T,S>,T>::value, int> = 0 >
matrix<n,m,T> operator*(matrix<n,m,T> mat, const S rhs) {
  mat *= rhs;
  return mat;
}

// scalar * matrix
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_mul<S,T>,T>::value, int> = 0 >
matrix<n,m,T> operator*(const S rhs, matrix<n,m,T> mat) {
  mat *= rhs;         // assume here commutativity!  Fails if the "scalar" is also a matrix
  return mat;
}

// matrix / scalar
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_div<T,S>,T>::value, int> = 0 >
matrix<n,m,T> operator/(matrix<n,m,T> mat, const S rhs) {
  mat /= rhs;
  return mat;
}


template <int n, int m, typename T>
std::ostream& operator<<(std::ostream &strm, const matrix<n,m,T> &A) {
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
inline auto norm_squared(matrix<n,m,T> & rhs){
  auto result = norm_squared(rhs.c[0][0]);
  for (int i=0; i<n; i++) for (int j=0; j<m; j++) if(i>0 || j>0) {
    result += norm_squared(rhs.c[i][j]);
  }
  return result;
}

#endif
