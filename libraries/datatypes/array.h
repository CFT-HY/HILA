#ifndef ARRAY_H_
#define ARRAY_H_

#include "datatypes/matrix.h"

////////////////////////////////////////////////////////////////
/// nxm Array type
////////////////////////////////////////////////////////////////
template <const int n, const int m, typename T>
class Array {
  public:
    static_assert(is_cmplx_or_arithmetic<T>::value, "Array requires Cmplx or arithmetic type");

    /// Same data as Matrix, a one dimensional array
    T c[n*m];

    /// std incantation for field types
    using base_type = typename base_type_struct<T>::type;

    /// define default constructors to ensure std::is_trivial
    Array() = default;
    ~Array() = default;
    Array(const Array<n,m,T> & v) = default;


    /// constructor from scalar - make this also explicit, consistency
    template <typename S, std::enable_if_t<is_assignable<T&,S>::value, int> = 0 >
    explicit inline Array(const S rhs) {
      for (int i=0; i<n*m; i++) {
        this->c[i] = rhs;
      }
    }

    // and make non-explicit constructor from 0
    inline Array(const std::nullptr_t & z) {
      for (int i=0; i<n*m; i++) c[i] = static_cast<T>(0);
    }


    /// Define constant methods rows(), columns() - may be useful in template code
    constexpr int rows() const { return n; }
    constexpr int columns() const { return m; }

    // define also method size() for 'vector arrays' and square arrays only!
    template <int q=n, int p=m, std::enable_if_t< q==1, int> = 0 >
    constexpr int size() const { return p; }

    template <int q=n, int p=m, std::enable_if_t< p==1, int> = 0 >
    constexpr int size() const { return q; }

    template <int q=n, int p=m, std::enable_if_t< q==p, int> = 0 >
    constexpr int size() const { return q; }

    /// access operators .e(i,j) and .e(i) from Matrix
    inline T  e(const int i, const int j) const { return c[i*m + j]; }
    /// standard access ops m.e(i,j) - assume T is small, as it should
    inline T& e(const int i, const int j) { return c[i*m + j]; }
    
    /// declare single e here too in case we have a vector
    /// (one size == 1)
    template <int q=n, int p=m,
              std::enable_if_t< (q == 1 || p == 1), int> = 0 >
    inline T e(const int i) const {
      return c[i]; 
    }

    template <int q=n, int p=m,
              std::enable_if_t< (q == 1 || p == 1), int> = 0 >
    inline T& e(const int i) { 
      return c[i]; 
    }


    // cast from array to matrix
    Matrix<n,m,T> & asMatrix() {
      return *reinterpret_cast<Matrix<n,m,T>*>(this);
    }
    const Matrix<n,m,T> asMatrix() const {
      return *reinterpret_cast<const Matrix<n,m,T>*>(this);
    }


    /// casting from one Array (number) type to another   
    /// TODO: CHECK AVX CONVERSIONS
    template <typename S, std::enable_if_t<is_assignable<S&,T>::value, int> = 0 >
    operator Array<n,m,S>() {
      Array<n,m,S> res;
      for (int i=0; i<n*m; i++) {
        res.c[i] = c[i];
      }
      return res;
    }

    /// unary -
    inline Array<n,m,T> operator-() const {
      Array<n,m,T> res;
      for (int i=0; i<n*m; i++) {
        res.c[i] = -c[i];
      }
      return res;
    }

    /// unary +
    inline Array<n,m,T> operator+() const {
      return *this;
    }

    /// Assign from scalar to array
    template <typename S, std::enable_if_t<is_assignable<T&,S>::value, int> = 0 >
    inline Array<n,m,T> & operator= (const S rhs) {
      for (int i=0; i<n*m; i++) {
        c[i] = rhs;
      }
      return *this;
    }

    /// add assign an Array
    template <typename S,
              std::enable_if_t<std::is_convertible<S,T>::value, int> = 0 >
    Array<n,m,T> & operator+=(const Array<n,m,S> & rhs){
      for (int i=0; i < n*m; i++) {
        c[i] += rhs.c[i];
      }
      return *this;
    }

    /// subtract assign an Array
    template <typename S,
              std::enable_if_t<std::is_convertible<S,T>::value, int> = 0 >
    Array<n,m,T> & operator-=(const Array<n,m,S> & rhs){
      for (int i=0; i < n*m; i++) {
        c[i] -= rhs.c[i];
      }
      return *this;
    }

    /// add assign type T and convertible
    template <typename S,
              std::enable_if_t<std::is_convertible<type_plus<T,S>,T>::value, int> = 0 >
    Array<n,m,T> & operator+=(const S rhs){
      for (int i=0; i<n*m; i++) {
        c[i] += rhs;
      }
      return *this;
    }

    /// subtract assign type T and convertible
    template <typename S,
              std::enable_if_t<std::is_convertible<type_minus<T,S>,T>::value, int> = 0 >
    Array<n,m,T> & operator-=(const S rhs){
      for (int i=0; i<n*m; i++) {
        c[i] -= rhs;
      }
      return *this;
    }

    /// multiply assign with Array
    template <typename S,
              std::enable_if_t<std::is_convertible<type_mul<T,S>,T>::value, int> = 0 >
    Array<n,m,T> & operator*=(const Array<n,m,S> & rhs){
      for (int i=0; i<n*m; i++) {
        c[i] *= rhs.c[i];
      }
      return *this;
    }

    /// multiply assign with scalar
    template <typename S,
              std::enable_if_t<std::is_convertible<type_mul<T,S>,T>::value, int> = 0 >
    Array<n,m,T> & operator*=(const S rhs) {
      for (int i=0; i<n*m; i++) {
        c[i] *= rhs;
      }
      return *this;
    }

    /// divide assign by Array
    template <typename S,
              std::enable_if_t<std::is_convertible<type_div<T,S>,T>::value, int> = 0 >
    Array<n,m,T> & operator/=(const Array<n,m,S> & rhs){
      for (int i=0; i<n*m; i++) {
        c[i] /= rhs.c[i];
      }
      return *this;
    }

    /// divide assign with scalar
    template <typename S,
              std::enable_if_t<std::is_convertible<type_div<T,S>,T>::value, int> = 0 >
    Array<n,m,T> & operator/=(const S rhs) {
      for (int i=0; i<n*m; i++) {
        c[i] /= rhs;
      }
      return *this;
    }

    /// complex conjugate
    inline Array<n,m,T> conj() const { 
      Array<n,m,T> res;
      for (int i=0; i<n*m; i++) {
        res.c[i] = ::conj(c[i]);
      }
      return res;
    }
        
    /// return real part
    inline Array<n,m,number_type<T>> real() const { 
      Array<n,m,number_type<T>> res;
      for (int i=0; i<m*n; i++) {
        res.c[i] = real(c[i]);
      }
      return res;
    }

    /// return imaginary part
    inline Array<n,m,number_type<T>> imag() const { 
      Array<n,m,number_type<T>> res;
      for (int i=0; i<m*n; i++) {
        res.c[i] = imag(c[i]);
      }
      return res;
    }

    /// calculate square norm - sum of squared elements
    number_type<T> norm_sq() const {
      number_type<T> result = 0;
      for (int i=0; i<n*m; i++) {
        result += norm_squared(c[i]);
      }
      return result;
    }

 
    /// Generate random elements
    Array<n, m, T> & random() {
      for (int i=0; i<n*m; i++) {
        ::random(c[i]);
      }
      return *this;
    }

    /// Generate gaussian random elements
    inline Array<n, m, T> & gaussian(){ 
      for (int i = 0; i < n*m; i++) {
        ::gaussian_random(c[i]);
      }
      return *this;
    }

    /// Convert to string for printing
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

/// conjugate
template <const int n, const int m, typename T>
inline Array<n,m,T> conj(const Array<n,m,T> & arg) {
  return arg.conj();
}
/// real part
template <const int n, const int m, typename T>
inline Array<n,m,number_type<T>> real(const Array<n,m,T> & arg) {
  return arg.real();
}
/// imaginary part
template <const int n, const int m, typename T>
inline Array<n,m,number_type<T>> imag(const Array<n,m,T> & arg) {
  return arg.imag();
}


/// Now Array additions: Array + Array
template <int n, int m, typename T>
inline Array<n,m,T> operator+(Array<n,m,T> a, const Array<n,m,T> & b){
  a += b;
  return a;
}

/// Array subtract
template <int n, int m, typename T>
inline Array<n,m,T> operator-(Array<n,m,T> a, const Array<n,m,T> & b){
  a -= b;
  return a;
}

/// Array + scalar
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_plus<T,S>,T>::value, int> = 0 >
inline Array<n,m,T> operator+(Array<n,m,T> a, const S b){
  a += b;
  return a;
}

/// scalar + Array
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_plus<T,S>,T>::value, int> = 0 >
inline Array<n,m,T> operator+(const S b, Array<n,m,T> a){
  a += b;
  return a;
}

/// Array - scalar
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_minus<T,S>,T>::value, int> = 0 >
Array<n,m,T> operator-(Array<n,m,T> a, const S b){
  a -= b;
  return a;
}

/// scalar - Array
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_minus<S,T>,T>::value, int> = 0 >
inline Array<n,m,T> operator-(const S b, Array<n,m,T> a){
  for (int i=0; i<n*m; i++) a.c[i] = static_cast<T>(b) - a.c[i];
  return a;
}

/// and Array*Array
template <int n, int m, typename T>
inline Array<n,m,T> operator*(Array<n,m,T> a, const Array<n,m,T> & b){
  a *= b;
  return a;
}

/// and Array/Array
template <int n, int m, typename T>
inline Array<n,m,T> operator/(Array<n,m,T> a, const Array<n,m,T> & b){
  a /= b;
  return a;
}

/// Array * scalar
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_mul<T,S>,T>::value, int> = 0 >
inline Array<n,m,T> operator*(Array<n,m,T> a, const S b){
  a *= b;
  return a;
}

/// scalar * Array
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_mul<T,S>,T>::value, int> = 0 >
inline Array<n,m,T> operator*(const S b, Array<n,m,T> a){
  a *= b;
  return a;
}

/// Array / scalar
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_div<T,S>,T>::value, int> = 0 >
inline Array<n,m,T> operator/(Array<n,m,T> a, const S b){
  a /= b;
  return a;
}

/// scalar / Array 
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_div<S,T>,T>::value, int> = 0 >
inline Array<n,m,T> operator/(const S b, Array<n,m,T> a){
  for (int i=0; i<n*m; i++) a.c[i] = b/a.c[i];
  return a;
}

/// Stream operator
template <int n, int m, typename T>
std::ostream& operator<<(std::ostream &strm, const Array<n,m,T> &A) {
  return operator<<( strm, A.asMatrix() );
}

/// Norm squared function
template<int n, int m, typename T>
inline number_type<T> norm_squared(const Array<n,m,T> & rhs){
  return rhs.norm_sq();
}

/// Function that calls random()-method
template<int n, int m, typename T>
inline void random(Array<n,m,T> & mat) {
  mat.random();
}

/// Function that calls the gaussian()-method
template<int n, int m, typename T>
inline void gaussian_random(Array<n,m,T> & mat) {
  mat.gaussian();
}


template <int n, typename T>
using Array1d = Array<n,1,T>;



#endif