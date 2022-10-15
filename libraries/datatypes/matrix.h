#ifndef MATRIX_H_
#define MATRIX_H_

#include <type_traits>
#include <sstream>
#include "plumbing/defs.h"
#include "datatypes/cmplx.h"

template <const int n, const int m, typename T, typename Mtype>
class Matrix_t;

template <const int n, const int m, typename T = double>
class Array;

///////////////////////////////////////////////////////////////
/// Generic Matrix class fwd definition
///////////////////////////////////////////////////////////////

template <int n, int m, typename T>
class Matrix;

///////////////////////////////////////////////////////////////
/// Define Vector and HorizontalVector as a special case of Matrix
///////////////////////////////////////////////////////////////

template <int n, typename T>
using Vector = Matrix<n, 1, T>;

template <int n, typename T>
using HorizontalVector = Matrix<1, n, T>;

///////////////////////////////////////////////////////////////
/// Define SquareMatrix as an alias
///////////////////////////////////////////////////////////////

template <int n, typename T>
using SquareMatrix = Matrix<n, n, T>;

// template <const int n, const int m, typename T>
// class DaggerMatrix;

// Special case - m.column(), column of a matrix (not used now)
// #include "matrix_column.h"

////////////////////////////////////////////////////////////////
/// The main nxm matrix type template Matrix_t
/// This is a root type, and "useful" types are derived from this type
///
/// Uses curiously recurring template pattern (CRTP), where
/// the last template parameter is the template itself
/// Example: the matrix type below is defined as
///   template <int n, int m, typename T>
///   class Matrix : public Matrix_t<n, m, T, Matrix<n, m, T>> { .. }
///
/// Used because stupid c++ makes it complicated to write
/// generic code, in this case derived functions to return derived type
///////////////////////////////////////////////////////////////

template <const int n, const int m, typename T, typename Mtype>
class Matrix_t {
  public:
    static_assert(hila::is_complex_or_arithmetic<T>::value,
                  "Matrix requires Complex or arithmetic type");

    /// std incantation for field types
    using base_type = hila::number_type<T>;
    using argument_type = T;

    // help for templates, can use T::is_matrix()
    static constexpr bool is_matrix() {
        return true;
    }

    /// The data as a one dimensional array
    // union {
    T c[n * m];
    //    T elem[n][m];
    //};

    /// define default constructors to ensure std::is_trivial
    Matrix_t() = default;
    ~Matrix_t() = default;
    Matrix_t(const Matrix_t &v) = default;

    /// constructor from scalar -- keep it explicit!  Not good for auto use
    template <typename S, int nn = n, int mm = m,
              std::enable_if_t<(hila::is_assignable<T &, S>::value && nn == mm), int> = 0>
    explicit inline Matrix_t(const S rhs) {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                if (i == j)
                    e(i, j) = rhs;
                else
                    e(i, j) = 0;
            }
    }

    /// construct from 0
    inline Matrix_t(const std::nullptr_t &z) {
        for (int i = 0; i < n * m; i++) {
            c[i] = 0;
        }
    }

    /// Construct matrix automatically from right-size initializer list
    /// This does not seem to be dangerous, so keep non-explicit

    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    inline Matrix_t(std::initializer_list<S> rhs) {
        assert(rhs.size() == n * m &&
               "Matrix/Vector initializer list size must match variable size");
        int i = 0;
        for (auto it = rhs.begin(); it != rhs.end(); it++, i++) {
            c[i] = *it;
        }
    }

    // cast to curious type
#pragma hila loop_function
    template <typename Tm = Mtype,
              std::enable_if_t<!std::is_same<Tm, Matrix<n, m, T>>::value, int> = 0>
    inline operator Mtype &() {
        return *reinterpret_cast<Mtype *>(this);
    }

#pragma hila loop_function
    template <typename Tm = Mtype,
              std::enable_if_t<!std::is_same<Tm, Matrix<n, m, T>>::value, int> = 0>
    inline operator const Mtype &() const {
        return *reinterpret_cast<const Mtype *>(this);
    }

    /// automatically cast to generic matrix
#pragma hila loop_function
    inline operator Matrix<n, m, T> &() {
        return *reinterpret_cast<Matrix<n, m, T> *>(this);
    }
#pragma hila loop_function
    inline operator const Matrix<n, m, T> &() const {
        return *reinterpret_cast<const Matrix<n, m, T> *>(this);
    }

    /// Define constant methods rows(), columns() - may be useful in template code
    static constexpr int rows() {
        return n;
    }
    static constexpr int columns() {
        return m;
    }

    // define also method size() for vectors and square matrices only!
    template <int q = n, int p = m, std::enable_if_t<q == 1, int> = 0>
    static constexpr int size() {
        return p;
    }

    template <int q = n, int p = m, std::enable_if_t<p == 1, int> = 0>
    static constexpr int size() {
        return q;
    }

    template <int q = n, int p = m, std::enable_if_t<q == p, int> = 0>
    static constexpr int size() {
        return q;
    }

    /// standard access ops m.e(i,j) - assume T is small, as it should
    inline T e(const int i, const int j) const {
        // return elem[i][j];
        return c[i * m + j];
    }
    /// standard access ops m.e(i,j) - assume T is small, as it should
    inline T &e(const int i, const int j) const_function {
        // return elem[i][j];
        return c[i * m + j];
    }

    /// declare single e here too in case we have a vector
    /// (one size == 1)
    template <int q = n, int p = m, std::enable_if_t<(q == 1 || p == 1), int> = 0>
    inline T e(const int i) const {
        return c[i];
    }

    template <int q = n, int p = m, std::enable_if_t<(q == 1 || p == 1), int> = 0>
    inline T &e(const int i) const_function {
        return c[i];
    }

    /// And also [] for vectors (not matrices!)
    /// (one size == 1)
    template <int q = n, int p = m, std::enable_if_t<(q == 1 || p == 1), int> = 0>
    inline T operator[](const int i) const {
        return c[i];
    }

    template <int q = n, int p = m, std::enable_if_t<(q == 1 || p == 1), int> = 0>
    inline T &operator[](const int i) const_function {
        return c[i];
    }

    // /// get row of a matrix
    // HorizontalVector<m, T> row(int r) const {
    //     HorizontalVector<m, T> v;
    //     for (int i = 0; i < m; i++)
    //         v[i] = e(r, i);
    //     return v;
    // }

    /// return reference to row in a matrix
    const HorizontalVector<m, T> &row(int r) const {
        return *(reinterpret_cast<const HorizontalVector<m, T> *>(this) + r);
    }

    /// return reference to row in a matrix
    HorizontalVector<m, T> &row(int r) {
        return *(reinterpret_cast<HorizontalVector<m, T> *>(this) + r);
    }

    /// return reference to row in a matrix
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    void set_row(int r, const HorizontalVector<m, S> &v) {
        for (int i = 0; i < m; i++)
            e(r, i) = v[i];
    }

    /// get column of a matrix
    const Vector<n, T> column(int c) const {
        Vector<n, T> v;
        for (int i = 0; i < n; i++)
            v[i] = e(i, c);
        return v;
    }

    /// get column of a matrix
    // hila_matrix_column_t<n, T, Mtype> column(int c) {
    //     return hila_matrix_column_t<n, T, Mtype>(*this, c);
    // }

    /// set column of a matrix
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    void set_column(int c, const Vector<n, S> &v) {
        for (int i = 0; i < n; i++)
            e(i, c) = v[i];
    }

    /// return diagonal of a square matrix as a vector
    Vector<n, T> diagonal() {
        static_assert(n == m, "diagonal() method defined only for square matrices");
        Vector<n, T> res;
        for (int i = 0; i < n; i++)
            res.e(i) = (*this).e(i, i);
        return res;
    }

    /// return diagonal of a square matrix as a vector
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    void set_diagonal(const Vector<n,S> & v) {
        static_assert(n == m, "set_diagonal() method defined only for square matrices");
        for (int i = 0; i < n; i++)
            (*this).e(i, i) = v.e(i);
    }


    /// interpret Matrix as Array -  for array ops
    Array<n, m, T> &asArray() const_function {
        return *reinterpret_cast<Array<n, m, T> *>(this);
    }
    const Array<n, m, T> &asArray() const {
        return *reinterpret_cast<const Array<n, m, T> *>(this);
    }

    /// casting from one Matrix (number) type to another: do not do this automatically.
    /// but require an explicit cast operator.  This makes it easier to write code.
    /// or should it be automatic?  keep/remove explicit?
    /// TODO: CHECK AVX CONVERSIONS

    // template <typename S, typename Rtype,
    //           std::enable_if_t<
    //               Rtype::is_matrix() && Rtype::rows() == n && Rtype::columns() == m
    //               &&
    //                   hila::is_assignable<typename (Rtype::argument_type) &,
    //                   T>::value,
    //               int> = 0>
    // explicit operator Rtype() const {
    //     Rtype res;
    //     for (int i = 0; i < n * m; i++)
    //         res.c[i] = c[i];
    //     return res;
    // }

    /// unary -
    inline Mtype operator-() const {
        Mtype res;
        for (int i = 0; i < n * m; i++) {
            res.c[i] = -c[i];
        }
        return res;
    }

    /// unary +  - there's an automatic cast here
    inline const Mtype &operator+() const {
        return *this;
    }

    /// assign from 0
#pragma hila loop_function
    inline Mtype &operator=(const std::nullptr_t &z) out_only {
        for (int i = 0; i < n * m; i++) {
            c[i] = 0;
        }
        return *this;
    }

    bool operator==(const Matrix<n, m, T> &rhs) const {
        hila::number_type<T> epsilon = 0;
        return ((*this) - rhs).squarenorm() <= epsilon;
    }

    /// Assign from different type matrix
#pragma hila loop_function
    template <typename S, typename MT,
              std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    inline Mtype &operator=(const Matrix_t<n, m, S, MT> &rhs) out_only {
        for (int i = 0; i < n * m; i++) {
            c[i] = rhs.c[i];
        }
        return *this;
    }

    /// Assign from "scalar" for square matrix
#pragma hila loop_function
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value && n == m, int> = 0>
    inline Mtype &operator=(const S rhs) out_only {

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                if (i == j)
                    e(i, j) = rhs;
                else
                    e(i, j) = 0;
            }
        return *this;
    }

    /// Assign from initializer list
#pragma hila loop_function
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    inline Mtype &operator=(std::initializer_list<S> rhs) out_only {
        assert(rhs.size() == n * m && "Initializer list has a wrong size in assignment");
        int i = 0;
        for (auto it = rhs.begin(); it != rhs.end(); it++, i++) {
            c[i] = *it;
        }
        return *this;
    }

    /// add assign a Matrix
#pragma hila loop_function
    template <typename S, typename MT,
              std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    Mtype &operator+=(const Matrix_t<n, m, S, MT> &rhs) {
        for (int i = 0; i < n * m; i++) {
            c[i] += rhs.c[i];
        }
        return *this;
    }

    /// subtract assign a Matrix_t
#pragma hila loop_function
    template <typename S, typename MT,
              std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    Mtype &operator-=(const Matrix_t<n, m, S, MT> &rhs) {
        for (int i = 0; i < n * m; i++) {
            c[i] -= rhs.c[i];
        }
        return *this;
    }

    /// add assign type T and convertible
#pragma hila loop_function
    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_plus<T, S>>::value, int> = 0>
    Mtype &operator+=(const S &rhs) {

        static_assert(n == m, "rows != columns : scalar addition possible for square matrix only!");

        for (int i = 0; i < n; i++) {
            e(i, i) += rhs;
        }
        return *this;
    }

    /// subtract assign type T and convertible
#pragma hila loop_function
    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_minus<T, S>>::value, int> = 0>
    Mtype &operator-=(const S rhs) {
        static_assert(n == m,
                      "rows != columns : scalar subtraction possible for square matrix only!");
        for (int i = 0; i < n; i++) {
            e(i, i) -= rhs;
        }
        return *this;
    }

    /// multiply assign with matrix
#pragma hila loop_function
    template <int p, typename S, typename MT,
              std::enable_if_t<hila::is_assignable<T &, hila::type_mul<T, S>>::value, int> = 0>
    Mtype &operator*=(const Matrix_t<m, p, S, MT> &rhs) {
        static_assert(m == p, "can't assign result of *= to lhs Matrix, because doing so "
                              "would change it's dimensions");
        *this = *this * rhs;
        return *this;
    }

    /// multiply assign with scalar
#pragma hila loop_function
    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_mul<T, S>>::value, int> = 0>
    Mtype &operator*=(const S rhs) {
        for (int i = 0; i < n * m; i++) {
            c[i] *= rhs;
        }
        return *this;
    }

    /// divide assign with scalar
#pragma hila loop_function
    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_div<T, S>>::value, int> = 0>
    Mtype &operator/=(const S rhs) {
        for (int i = 0; i < n * m; i++) {
            c[i] /= rhs;
        }
        return *this;
    }

    /// numpy style matrix fill
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    const Mtype &fill(const S rhs) out_only {
        for (int i = 0; i < n * m; i++)
            c[i] = rhs;
        return *this;
    }

    ///////////////////////////////////////////////////////////////////
    /// Transpose of a matrix: For square the return type is the
    /// same as input, for non-square general Matrix<m,n>
    template <int mm = m,
              typename Rtype = typename std::conditional<n == m, Mtype, Matrix<m, n, T>>::type,
              std::enable_if_t<(mm != 1), int> = 0>
    inline const Rtype transpose() const {
        Rtype res;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) {
                res.e(j, i) = e(i, j);
            }
        return res;
    }

    /// Transpose of a vector: just return a ref
    template <int mm = m, std::enable_if_t<mm == 1, int> = 0>
    inline const HorizontalVector<n, T> &transpose() const {
        return *reinterpret_cast<const HorizontalVector<n, T> *>(this);
    }

    ////////////////////////////////////////////////////////////////////
    /// and complex conjugate - same type
    inline const Mtype conj() const {
        Mtype res;
        for (int i = 0; i < n * m; i++) {
            res.c[i] = ::conj(c[i]);
        }
        return res;
    }

    ///////////////////////////////////////////////////////////////////
    /// Hermitean conjugate of a matrix: For square the return type is the
    /// same as input, for non-square general Matrix<m,n>
    template <typename Rtype = typename std::conditional<n == m, Mtype, Matrix<m, n, T>>::type>
    inline const Rtype dagger() const {
        Rtype res;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) {
                res.e(j, i) = ::conj(e(i, j));
            }
        return res;
    }

    /// alias dagger() to adjoint()
    template <typename Rtype = typename std::conditional<n == m, Mtype, Matrix<m, n, T>>::type>
    inline Rtype adjoint() const {
        return dagger();
    }

    /// abs() - give absolute value of elements (real/complex)
    /// With complex type changes!
    template <typename M = T, std::enable_if_t<!hila::contains_complex<M>::value, int> = 0>
    Mtype abs() const {
        Mtype res;
        for (int i = 0; i < n * m; i++) {
            // need to use ::abs for generic real vars
            res.c[i] = ::abs(c[i]);
        }
        return res;
    }

    template <typename M = T, std::enable_if_t<hila::contains_complex<M>::value, int> = 0>
    auto abs() const {
        Matrix<n, m, hila::number_type<T>> res;
        for (int i = 0; i < n * m; i++) {
            res.c[i] = c[i].abs();
        }
        return res;
    }


    // It seems that using special "Dagger" type makes the code slower!
    // Disable it now
    // inline const DaggerMatrix<m,n,T> & dagger() const {
    //     return *reinterpret_cast<const DaggerMatrix<m,n,T>*>(this);
    // }

    /// return real part
    inline Matrix<n, m, hila::number_type<T>> real() const {
        Matrix<n, m, hila::number_type<T>> res;
        for (int i = 0; i < m * n; i++) {
            res.c[i] = ::real(c[i]);
        }
        return res;
    }

    /// return imaginary part
    inline Matrix<n, m, hila::number_type<T>> imag() const {
        Matrix<n, m, hila::number_type<T>> res;
        for (int i = 0; i < m * n; i++) {
            res.c[i] = ::imag(c[i]);
        }
        return res;
    }

    /// matrix trace (return type T)
    T trace() const {
        static_assert(n == m, "trace not defined for non square matrices!");
        T result(0);
        for (int i = 0; i < n; i++) {
            result += e(i, i);
        }
        return result;
    }

    /// Multiply this (nxm) matrix with mxn matrix and return trace
    /// Cheaper than explicit (a*b).trace()
    template <int p, int q, typename S, typename MT>
    hila::type_mul<T, S> mul_trace(const Matrix_t<p, q, S, MT> &rm) const {

        static_assert(p == m && q == n, "mul_trace(): argument matrix size mismatch");

        hila::type_mul<T, S> res = 0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) {
                res += e(i, j) * rm.e(j, i);
            }
        return res;
    }

    /// calculate square norm - sum of squared elements
    hila::number_type<T> squarenorm() const {
        hila::number_type<T> result(0);
        for (int i = 0; i < n * m; i++) {
            result += ::squarenorm(c[i]);
        }
        return result;
    }

    /// calculate vector norm - sqrt of squarenorm
    template <typename S = T,
              std::enable_if_t<hila::is_floating_point<hila::number_type<S>>::value, int> = 0>
    hila::number_type<T> norm() const {
        return sqrt(squarenorm());
    }

    template <typename S = T,
              std::enable_if_t<!hila::is_floating_point<hila::number_type<S>>::value, int> = 0>
    double norm() const {
        return sqrt(static_cast<double>(squarenorm()));
    }


    /// dot product - (*this).dagger() * rhs
    /// could be done as well by writing the operation as above!
    template <int p, int q, typename S, typename R = hila::type_mul<T, S>>
    inline R dot(const Matrix<p, q, S> &rhs) const {
        static_assert(m == 1 && q == 1 && p == n,
                      "dot() product only for vectors of the same length");

        R r = 0;
        for (int i = 0; i < n; i++) {
            r += ::conj(c[i]) * rhs.c[i];
        }
        return r;
    }

    /// outer product - (*this) * rhs.dagger(), sizes (n,1) and (p,1)
    /// gives n * p matrix
    /// could be done as well by the above operation!
    template <int p, int q, typename S, typename R = hila::type_mul<T, S>>
    inline Matrix<n, p, R> outer_product(const Matrix<p, q, S> &rhs) const {
        static_assert(m == 1 && q == 1, "outer_product() only for vectors");

        Matrix<n, p, R> res;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < p; j++) {
                res.e(i, j) = c[i] * ::conj(rhs.c[j]);
            }
        }
        return res;
    }


    /// dot with matrix - matrix
    // template <int p, std::enable_if_t<(p > 1 || m > 1), int> = 0>
    // inline Matrix<m, p, T> dot(const Matrix<n, p, T> &rhs) const {
    //     Matrix<m, p, T> res;
    //     for (int i = 0; i < m; i++)
    //         for (int j = 0; j < p; j++) {
    //             res.e(i, j) = 0;
    //             for (int k = 0; k < n; j++)
    //                 res.e(i, j) += ::conj(e(k, i)) * rhs.e(k, j);
    //         }
    // }

    /// Generate random elements
    Mtype &random() out_only {

        static_assert(hila::is_floating_point<hila::number_type<T>>::value,
                      "Matrix/Vector random() requires non-integral type elements");

        for (int i = 0; i < n * m; i++) {
            ::random(c[i]);
        }
        return *this;
    }

    /// Generate gaussian random elements
    inline Mtype &gaussian_random(double width = 1.0) out_only {

        static_assert(hila::is_floating_point<hila::number_type<T>>::value,
                      "Matrix/Vector gaussian_random() requires non-integral type elements");

        // for Complex numbers gaussian_random fills re and im efficiently
        if constexpr (hila::is_complex<T>::value) {
            for (int i = 0; i < n * m; i++) {
                c[i].gaussian_random(width);
            }
        } else {
            // now not complex matrix
            // if n*m even, max i in loop below is n*m-2.
            // if n*m odd, max i is n*m-3
            double gr;
            for (int i = 0; i < n * m - 1; i += 2) {
                c[i] = hila::gaussrand2(gr) * width;
                c[i + 1] = gr * width;
            }
            if ((n * m) % 2 > 0) {
                c[n * m - 1] = hila::gaussrand() * width;
            }
        }
        return *this;
    }

    /// Convert to string for printing
    std::string str() const {
        std::stringstream text;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                text << e(i, j) << " ";
            }
            text << '\n';
        }
        return text.str();
    }
};

//////////////////////////////////////////////////////////////////////////
/// The matrix class definition here
//////////////////////////////////////////////////////////////////////////

template <int n, int m, typename T>
class Matrix : public Matrix_t<n, m, T, Matrix<n, m, T>> {

  public:
    /// std incantation for field types
    using base_type = hila::number_type<T>;
    using argument_type = T;

    // use the Base::Base -trick to inherit constructors and assignments
    using Matrix_t<n, m, T, Matrix<n, m, T>>::Matrix_t;
    using Matrix_t<n, m, T, Matrix<n, m, T>>::operator=;
    using Matrix_t<n, m, T, Matrix<n, m, T>>::operator+=;
    using Matrix_t<n, m, T, Matrix<n, m, T>>::operator-=;
    using Matrix_t<n, m, T, Matrix<n, m, T>>::operator*=;
    using Matrix_t<n, m, T, Matrix<n, m, T>>::operator/=;
};

namespace hila {

//////////////////////////////////////////////////////////////////////////
// Tool to get "right" result type for matrix (T1) + (T2) -op, where
// T1 and T2 are either complex or arithmetic matrices
// Use: hila::mat_x_mat_type<T1,T2>
//   - T1 == T2 returns T1
//   - If T1 or T2 contains complex, returns Matrix<rows,cols,Complex<typeof T1+T2>>
//   - otherwise returns the "larger" type

template <typename T1, typename T2>
struct matrix_op_type_s {
    using type = Matrix<T1::rows(), T1::columns(), hila::ntype_op<T1, T2>>;
};

// if types are the same
template <typename T>
struct matrix_op_type_s<T, T> {
    using type = T;
};

// useful format mat_x_mat_type<T1,T2>
template <typename T1, typename T2>
using mat_x_mat_type = typename matrix_op_type_s<T1, T2>::type;

////////////////////////////////////////////////////////////////////////////////
// Matrix + scalar result type:
// hila::mat_scalar_type<Mt,S>
//  - if result is convertible to Mt, return Mt
//  - if Mt is not complex and S is, return
//  Matrix<Complex<type_sum(number_type(Mt),number_type(S))>>
//  - otherwise return Matrix<type_sum>

template <typename Mt, typename S, typename Enable = void>
struct matrix_scalar_op_s {
    using type = Matrix<Mt::rows(), Mt::columns(),
                        Complex<hila::type_plus<hila::number_type<Mt>, hila::number_type<S>>>>;
};

template <typename Mt, typename S>
struct matrix_scalar_op_s<
    Mt, S,
    typename std::enable_if_t<std::is_convertible<hila::type_plus<hila::underlying_type<Mt>, S>,
                                                  hila::underlying_type<Mt>>::value>> {
    // using type = Mt;
    using type = typename std::conditional<
        hila::is_floating_point<hila::number_type<Mt>>::value, Mt,
        Matrix<Mt::rows(), Mt::columns(),
               hila::type_plus<hila::number_type<Mt>, hila::number_type<S>>>>::type;
};

template <typename Mt, typename S>
using mat_scalar_type = typename matrix_scalar_op_s<Mt, S>::type;


} // namespace hila

//////////////////////////////////////////////////////////////////////////

/// do transpose and adjoint functions here
template <typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0>
inline auto transpose(const Mtype &arg) {
    return arg.transpose();
}

/// conjugate
template <typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0>
inline auto conj(const Mtype &arg) {
    return arg.conj();
}

/// adjoint (=dagger)
template <typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0>
inline auto adjoint(const Mtype &arg) {
    return arg.adjoint();
}

/// dagger (=adjoint)
template <typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0>
inline auto dagger(const Mtype &arg) {
    return arg.adjoint();
}

/// absolute value
template <typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0>
inline auto abs(const Mtype &arg) {
    return arg.abs();
}

/// trace
template <typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0>
inline auto trace(const Mtype &arg) {
    return arg.trace();
}

/// real part
template <typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0>
inline auto real(const Mtype &arg) {
    return arg.real();
}

/// imaginary part
template <typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0>
inline auto imag(const Mtype &arg) {
    return arg.imag();
}

/// templates needed for naive calculation of determinants
template <
    typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0,
    typename Rtype = Matrix<Mtype::rows() - 1, Mtype::columns() - 1, hila::underlying_type<Mtype>>>
Rtype Minor(const Mtype &bigger, int row, int col) {
    constexpr int n = Mtype::rows();
    constexpr int m = Mtype::columns();

    Rtype result;
    int ri = 0;
    for (int i = 0; i < n; i++)
        if (i != row) {
            int ci = 0;
            for (int j = 0; j < m; j++) {
                if (j != col) {
                    result.e(ri, ci) = bigger.e(i, j);
                    ci++;
                }
            }
            ri++;
        }

    return result;
}

/// determinant -> use LU factorization later
template <typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0,
          typename T = hila::underlying_type<Mtype>, int n = Mtype::rows(),
          int m = Mtype::columns(), std::enable_if_t<(n > 3), int> = 0>
T det(const Mtype &mat) {
    static_assert(n == m, "determinants defined only for square matrices");
    T result(0);
    hila::number_type<T> parity = 1, opposite = -1;
    for (int i = 0; i < n; i++) {
        Matrix<n - 1, m - 1, T> minor = Minor(mat, 0, i);
        result += parity * det(minor) * mat.e(0, i);
        parity *= opposite;
    }
    return result;
}

/// 1x1 matrix det (trivial)
template <
    typename Mtype,
    std::enable_if_t<Mtype::is_matrix() && Mtype::rows() == 1 && Mtype::columns() == 1, int> = 0>
auto det(const Mtype &mat) {
    return mat.e(0, 0);
}

/// 2x2 matrix det
template <
    typename Mtype,
    std::enable_if_t<Mtype::is_matrix() && Mtype::rows() == 2 && Mtype::columns() == 2, int> = 0>
auto det(const Mtype &mat) {
    return mat.e(0, 0) * mat.e(1, 1) - mat.e(1, 0) * mat.e(0, 1);
}

/// 3x3 matrix det
template <
    typename Mtype,
    std::enable_if_t<Mtype::is_matrix() && Mtype::rows() == 3 && Mtype::columns() == 3, int> = 0>
auto det(const Mtype &mat) {
    return mat.e(0, 0) * (mat.e(1, 1) * mat.e(2, 2) - mat.e(2, 1) * mat.e(1, 2)) -
           mat.e(0, 1) * (mat.e(1, 0) * mat.e(2, 2) - mat.e(1, 2) * mat.e(2, 0)) +
           mat.e(0, 2) * (mat.e(1, 0) * mat.e(2, 1) - mat.e(1, 1) * mat.e(2, 0));
}

///////////////////////////////////////////////////////////////////////////
/// Now matrix additions: matrix + matrix

template <typename Mtype1, typename Mtype2,
          std::enable_if_t<Mtype1::is_matrix() && Mtype2::is_matrix(), int> = 0,
          typename Rtype = hila::mat_x_mat_type<Mtype1, Mtype2>,
          std::enable_if_t<!std::is_same<Mtype1, Rtype>::value, int> = 0>
inline Rtype operator+(const Mtype1 &a, const Mtype2 &b) {

    constexpr int n = Mtype1::rows();
    constexpr int m = Mtype1::columns();

    static_assert(n == Mtype2::rows() && m == Mtype2::columns(), "Matrix sizes do not match");

    Rtype r;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            r.e(i, j) = a.e(i, j) + b.e(i, j);
    return r;
}

// Real micro-optimization wrt. above - no extra creation of variable and copy.

template <typename Mtype1, typename Mtype2,
          std::enable_if_t<Mtype1::is_matrix() && Mtype2::is_matrix(), int> = 0,
          typename Rtype = hila::mat_x_mat_type<Mtype1, Mtype2>,
          std::enable_if_t<std::is_same<Mtype1, Rtype>::value, int> = 0>
inline Rtype operator+(Mtype1 a, const Mtype2 &b) {

    constexpr int n = Mtype1::rows();
    constexpr int m = Mtype1::columns();

    static_assert(n == Mtype2::rows() && m == Mtype2::columns(), "Matrix sizes do not match");

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            a.e(i, j) += b.e(i, j);
    return a;
}


/// Matrix - matrix
template <typename Mtype1, typename Mtype2,
          std::enable_if_t<Mtype1::is_matrix() && Mtype2::is_matrix(), int> = 0,
          typename Rtype = hila::mat_x_mat_type<Mtype1, Mtype2>>
inline Rtype operator-(const Mtype1 &a, const Mtype2 &b) {

    constexpr int n = Mtype1::rows();
    constexpr int m = Mtype1::columns();

    static_assert(n == Mtype2::rows() && m == Mtype2::columns(), "Matrix sizes do not match");

    Rtype r;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            r.e(i, j) = a.e(i, j) - b.e(i, j);
    return r;
}

/// Matrix + scalar
template <typename Mtype, typename S,
          std::enable_if_t<Mtype::is_matrix() && hila::is_complex_or_arithmetic<S>::value, int> = 0,
          typename Rtype = hila::mat_scalar_type<Mtype, S>>
inline Rtype operator+(const Mtype &a, const S &b) {

    static_assert(Mtype::rows() == Mtype::columns(),
                  "Matrix + scalar possible only for square matrix");

    Rtype r;
    for (int i = 0; i < Rtype::rows(); i++)
        for (int j = 0; j < Rtype::columns(); j++) {
            r.e(i, j) = a.e(i, j);
            if (i == j)
                r.e(i, j) += b;
        }
    return r;
}

/// scalar + matrix
template <typename Mtype, typename S,
          std::enable_if_t<Mtype::is_matrix() && hila::is_complex_or_arithmetic<S>::value, int> = 0,
          typename Rtype = hila::mat_scalar_type<Mtype, S>>
inline Rtype operator+(const S &b, const Mtype &a) {

    static_assert(Mtype::rows() == Mtype::columns(),
                  "Matrix + scalar possible only for square matrix");
    return a + b;
}

/// matrix - scalar
template <typename Mtype, typename S,
          std::enable_if_t<Mtype::is_matrix() && hila::is_complex_or_arithmetic<S>::value, int> = 0,
          typename Rtype = hila::mat_scalar_type<Mtype, S>>
inline Rtype operator-(const Mtype &a, const S &b) {

    static_assert(Mtype::rows() == Mtype::columns(),
                  "Matrix - scalar possible only for square matrix");

    Rtype r;
    for (int i = 0; i < Rtype::rows(); i++)
        for (int j = 0; j < Rtype::columns(); j++) {
            r.e(i, j) = a.e(i, j);
            if (i == j)
                r.e(i, j) -= b;
        }
    return r;
}

/// scalar - matrix
template <typename Mtype, typename S,
          std::enable_if_t<Mtype::is_matrix() && hila::is_complex_or_arithmetic<S>::value, int> = 0,
          typename Rtype = hila::mat_scalar_type<Mtype, S>>
inline Rtype operator-(const S &b, const Mtype &a) {

    static_assert(Mtype::rows() == Mtype::columns(),
                  "Scalar - matrix possible only for square matrix");

    Rtype r;
    for (int i = 0; i < Rtype::rows(); i++)
        for (int j = 0; j < Rtype::columns(); j++) {
            r.e(i, j) = -a.e(i, j);
            if (i == j)
                r.e(i, j) += b;
        }
    return r;
}

////////////////////////////////////////
/// matrix * matrix is the most important bit

// same type square matrices:
template <typename Mt, std::enable_if_t<Mt::is_matrix() && Mt::rows() == Mt::columns(), int> = 0>
inline Mt operator*(const Mt &A, const Mt &B) {

    constexpr int n = Mt::rows();

    Mt res;

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            res.e(i, j) = A.e(i, 0) * B.e(0, j);
            for (int k = 1; k < n; k++) {
                res.e(i, j) += A.e(i, k) * B.e(k, j);
            }
        }
    return res;
}

// different type matrices - return generic
template <typename Mt1, typename Mt2,
          std::enable_if_t<Mt1::is_matrix() && Mt2::is_matrix() && !std::is_same<Mt1, Mt2>::value,
                           int> = 0,
          typename R = hila::ntype_op<hila::underlying_type<Mt1>, hila::underlying_type<Mt2>>,
          int n = Mt1::rows(), int m = Mt2::columns()>
inline Matrix<n, m, R> operator*(const Mt1 &A, const Mt2 &B) {

    constexpr int p = Mt1::columns();
    static_assert(p == Mt2::rows(), "Matrix size: LHS columns != RHS rows");

    Matrix<n, m, R> res;

    if constexpr (n > 1 && m > 1) {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) {
                res.e(i, j) = A.e(i, 0) * B.e(0, j);
                for (int k = 1; k < p; k++) {
                    res.e(i, j) += A.e(i, k) * B.e(k, j);
                }
            }
    } else if constexpr (m == 1) {
        // matrix x vector
        for (int i = 0; i < n; i++) {
            res.e(i) = A.e(i, 0) * B.e(0);
            for (int k = 1; k < p; k++) {
                res.e(i) += A.e(i, k) * B.e(k);
            }
        }
    } else if constexpr (n == 1) {
        // vector x matrix
        for (int j = 0; j < m; j++) {
            res.e(j) = A.e(0) * B.e(0, j);
            for (int k = 1; k < p; k++) {
                res.e(j) += A.e(k) * B.e(k, j);
            }
        }
    }

    return res;
}

/// and treat separately horiz. vector * vector
template <int m, int n, typename T1, typename T2, typename Rt = hila::ntype_op<T1, T2>>
Rt operator*(const Matrix<1, m, T1> &A, const Matrix<n, 1, T2> &B) {

    static_assert(m == n, "Vector lengths do not match");

    Rt res;
    res = A.e(0) * B.e(0);
    for (int i = 1; i < m; i++) {
        res += A.e(i) * B.e(i);
    }
    return res;
}

///////////////////////////////////////////////////////////////////////////

/// matrix * scalar
template <typename Mt, typename S,
          std::enable_if_t<(Mt::is_matrix() && hila::is_complex_or_arithmetic<S>::value), int> = 0,
          typename Rt = hila::mat_scalar_type<Mt, S>>
inline Rt operator*(const Mt &mat, const S &rhs) {
    Rt res;
    for (int i = 0; i < Rt::rows() * Rt::columns(); i++) {
        res.c[i] = mat.c[i] * rhs;
    }
    return res;
}

/// scalar * matrix
template <typename Mt, typename S,
          std::enable_if_t<(Mt::is_matrix() && hila::is_complex_or_arithmetic<S>::value), int> = 0,
          typename Rt = hila::mat_scalar_type<Mt, S>>
inline Rt operator*(const S &rhs, const Mt &mat) {
    return mat * rhs; // assume commutes
}

/// matrix / scalar
template <typename Mt, typename S,
          std::enable_if_t<(Mt::is_matrix() && hila::is_complex_or_arithmetic<S>::value), int> = 0,
          typename Rt = hila::mat_scalar_type<Mt, S>>
inline Rt operator/(const Mt &mat, const S &rhs) {
    Rt res;
    for (int i = 0; i < Rt::rows() * Rt::columns(); i++) {
        res.c[i] = mat.c[i] / rhs;
    }
    return res;
}


/// mul_trace(a,b) - multiply matrices a and b and take trace
template <typename Mtype1, typename Mtype2,
          std::enable_if_t<Mtype1::is_matrix() && Mtype2::is_matrix(), int> = 0>
inline auto mul_trace(const Mtype1 &a, const Mtype2 &b) {

    static_assert(Mtype1::columns() == Mtype2::rows() && Mtype1::rows() == Mtype2::columns(),
                  "mul_trace(a,b): matrix sizes are not compatible");
    return a.mul_trace(b);
}

//////////////////////////////////////////////////////////////////////////////////


/// Stream operator
template <int n, int m, typename T, typename MT>
std::ostream &operator<<(std::ostream &strm, const Matrix_t<n, m, T, MT> &A) {
    if constexpr (n == 1 || m == 1) {
        // print a vector, horizontally
        strm << '[';
        for (int i = 0; i < n * m; i++)
            strm << ' ' << A.e(i);
        strm << " ]";
        // if (n > 1)
        //    strm << "^T";
    } else {
        // now a matrix - split the matrix on lines.
        // do it so that columns are equal width...
        std::vector<std::string> lines, columns;
        lines.resize(n);
        columns.resize(n);

        for (int i = 0; i < n; i++)
            lines[i] = "[ ";
        for (int j = 0; j < m; j++) {
            int size = 0;
            for (int i = 0; i < n; i++) {
                std::stringstream item;
                item << A.e(i, j);
                columns[i] = item.str();
                if (columns[i].size() > size)
                    size = columns[i].size();
            }
            for (int i = 0; i < n; i++) {
                lines[i].append(size - columns[i].size(), ' ');
                lines[i].append(columns[i]);
                lines[i].append(1, ' ');
            }
        }
        for (int i = 0; i < n; i++) {
            strm << lines[i] << "]\n";
        }
    }
    return strm;
}

/// Norm squared function
template <typename Mt, std::enable_if_t<Mt::is_matrix(), int> = 0>
inline auto squarenorm(const Mt &rhs) {
    return rhs.squarenorm();
}

/// Vector norm - sqrt of squarenorm()
template <typename Mt, std::enable_if_t<Mt::is_matrix(), int> = 0>
inline auto norm(const Mt &rhs) {
    return rhs.norm();
}

/// Function that calls random()-method
template <typename Mt,
          std::enable_if_t<Mt::is_matrix() && hila::is_floating_point<hila::number_type<Mt>>::value,
                           int> = 0>
inline void random(out_only Mt &mat) {
    mat.random();
}

/// Function that calls the gaussian_random()-method
template <typename Mt,
          std::enable_if_t<Mt::is_matrix() && hila::is_floating_point<hila::number_type<Mt>>::value,
                           int> = 0>
inline void gaussian_random(out_only Mt &mat, double width = 1.0) {
    mat.gaussian_random(width);
}

/// find determinant using LU decomposition. Algorithm: numerical Recipes, 2nd ed.
/// p. 47 ff
template <int n, int m, typename T, typename MT, typename radix = hila::number_type<T>,
          std::enable_if_t<hila::is_complex<T>::value, int> = 0>
Complex<radix> det_lu(const Matrix_t<n, m, T, MT> &mat) {

    static_assert(n == m, "Det_lu only for square matrix");

    int i, imax, j, k;
    radix big, d, temp, dum;
    Complex<radix> cdum, csum, ctmp1;
    radix vv[n];
    Complex<radix> a[n][n];
    Complex<radix> one;

    one = Complex<radix>(1, 0);
    d = 1;
    imax = -1;

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            a[i][j] = mat.e(i, j);
    for (i = 0; i < n; i++) {
        big = 0;
        for (j = 0; j < n; j++) {
            if ((temp = a[i][j].abs()) > big)
                big = temp;
        }
        assert(big != 0.0 && "Determinant does not exist\n");
        vv[i] = 1.0 / big;
    }

    for (j = 0; j < n; j++) {
        for (i = 0; i < j; i++) {
            csum = a[i][j];
            for (k = 0; k < i; k++) {
                csum -= a[i][k] * a[k][j];
            }
            a[i][j] = csum;
        }

        big = 0;
        for (i = j; i < n; i++) {
            csum = a[i][j];
            for (k = 0; k < j; k++) {
                csum -= a[i][k] * a[k][j];
            }
            a[i][j] = csum;
            if ((dum = vv[i] * csum.abs()) >= big) {
                big = dum;
                imax = i;
            }
        }

        if (j != imax) {
            for (k = 0; k < n; k++) {
                cdum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = cdum;
            }
            d = -d;
            vv[imax] = vv[j];
        }

        if (a[j][j].abs() == static_cast<radix>(0.0))
            a[j][j] = Complex<radix>(1e-20, 0);

        if (j != n - 1) {
            cdum = one / a[j][j];
            for (i = j + 1; i < n; i++) {
                a[i][j] = a[i][j] * cdum;
            }
        }
    }

    csum = Complex<radix>(d, 0.0);
    for (j = 0; j < n; j++) {
        csum = csum * a[j][j];
    }

    return (csum);
}

/// find determinant using LU decomposition. Algorithm: numerical Recipes, 2nd ed.
/// p. 47 ff
template <int n, int m, typename T, typename MT,
          std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
T det_lu(const Matrix_t<n, m, T, MT> &mat) {
    static_assert(n == m, "Det_lu only for square matrix");

    int i, imax, j, k;
    T big, d, temp, dum;
    T cdum, csum, ctmp1;
    T vv[n];
    T a[n][n];
    d = 1;
    imax = -1;

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            a[i][j] = mat.e(i, j);
    for (i = 0; i < n; i++) {
        big = 0;
        for (j = 0; j < n; j++) {
            if ((temp = a[i][j]) > big)
                big = temp;
        }
        assert(big != 0.0 && "Determinant does not exist\n");
        vv[i] = 1.0 / big;
    }

    for (j = 0; j < n; j++) {
        for (i = 0; i < j; i++) {
            csum = a[i][j];
            for (k = 0; k < i; k++) {
                csum -= a[i][k] * a[k][j];
            }
            a[i][j] = csum;
        }

        big = 0;
        for (i = j; i < n; i++) {
            csum = a[i][j];
            for (k = 0; k < j; k++) {
                csum -= a[i][k] * a[k][j];
            }
            a[i][j] = csum;
            if ((dum = vv[i] * csum) >= big) {
                big = dum;
                imax = i;
            }
        }

        if (j != imax) {
            for (k = 0; k < n; k++) {
                cdum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = cdum;
            }
            d = -d;
            vv[imax] = vv[j];
        }

        if (a[j][j] == static_cast<T>(0.0))
            a[j][j] = 1e-20;

        if (j != n - 1) {
            cdum = 1.0 / a[j][j];
            for (i = j + 1; i < n; i++) {
                a[i][j] = a[i][j] * cdum;
            }
        }
    }

    csum = d;
    for (j = 0; j < n; j++) {
        csum = csum * a[j][j];
    }

    return (csum);
}

// Cast operators to different number type
// cast_to<double>(a);

template <typename Ntype, typename T, int n, int m,
          std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
Matrix<n, m, Ntype> cast_to(const Matrix<n, m, T> &mat) {
    Matrix<n, m, Ntype> res;
    for (int i = 0; i < n * m; i++)
        res.c[i] = mat.c[i];
    return res;
}

template <typename Ntype, typename T, int n, int m,
          std::enable_if_t<hila::is_complex<T>::value, int> = 0>
Matrix<n, m, Complex<Ntype>> cast_to(const Matrix<n, m, T> &mat) {
    Matrix<n, m, Complex<Ntype>> res;
    for (int i = 0; i < n * m; i++)
        res.c[i] = cast_to<Ntype>(mat.c[i]);
    return res;
}

///  Calculate exp of a square matrix
///  Go to  order ORDER in the exponential of the matrices
///  matrices, since unitarity seems important.
///  Evaluation is done as:
/// 	exp(H) = 1 + H + H^2/2! + H^3/3! ..-
///	           = 1 + H*( 1 + (H/2)*( 1 + (H/3)*( ... )))
///  Do it backwards in order to reduce accumulation of errors

template <int n, int m, typename T, typename MT>
inline Matrix_t<n, m, T, MT> exp(const Matrix_t<n, m, T, MT> &mat, const int order = 20) {
    static_assert(n == m, "exp() only for square matrices");

    Matrix_t<n, m, T, MT> r;
    hila::number_type<T> one = 1.0;

    r = mat * (one / order) + one;

    // doing the loop step-by-step should reduce temporaries
    for (int k = order - 1; k > 1; k--) {
        r *= mat;
        r *= (one / k);
        r += one;
    }
    r *= mat;
    r += one;

    return r;
}

#include "datatypes/array.h"

// #include "datatypes/dagger.h"

#endif
