/**
 * @file matrix.h
 * @brief Definition of Matrix types
 * @details This file contains base matrix type Matrix_t which defines all general matrix type
 opirations Matrix types are Matrix, #Vector, #RowVector, #SquareMatrix of which Matrix is
 defined as a class and the rest are special cases of the Matrix class.
 *
 */

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

template <int n, int m, typename T>
class Matrix;

/**
 * @brief Vector is defined as 1-column Matrix
 */
template <int n, typename T>
using Vector = Matrix<n, 1, T>;

/**
 * @brief RowVector is a 1-row Matrix
 */
template <int n, typename T>
using RowVector = Matrix<1, n, T>;

/**
 * @brief Square matrix is defined as alias with special case of Matrix<n,n,T>
 */
template <int n, typename T>
using SquareMatrix = Matrix<n, n, T>;

// template <const int n, const int m, typename T>
// class DaggerMatrix;

// Special case - m.column(), column of a matrix (not used now)
// #include "matrix_column.h"


/**
 * @brief The main \f$ n \times m \f$ matrix type template Matrix_t. This is a root type, and
 * "useful" types are derived from this class
 *
 * @details Uses curiously recurring template pattern (CRTP), where the last template parameter is
 * the template itself
 *
 * Example: the matrix type below is defined as
 * @code{.cpp}
 * template <int n, int m, typename T>
 * class Matrix : public Matrix_t<n, m, T, Matrix<n, m, T>> { .. }
 * @endcode
 *
 * Used because stupid c++ makes it complicated to write generic code, in this case derived
 * functions to return derived type
 * @tparam n row length
 * @tparam m column length
 * @tparam T Data type Matrix
 * @tparam Mtype Specific "Matrix" type for CRTP
 */
template <const int n, const int m, typename T, typename Mtype>
class Matrix_t {

  public:
    /// The data as a one dimensional array
    T c[n * m];

  public:
    static_assert(hila::is_complex_or_arithmetic<T>::value,
                  "Matrix requires Complex or arithmetic type");

    // std incantation for field types
    using base_type = hila::scalar_type<T>;
    using argument_type = T;

    // help for templates, can use T::is_matrix()
    // Not very useful outside template parameters
    static constexpr bool is_matrix() {
        return true;
    }

    // is_vector and is_square bools
    static constexpr bool is_vector() {
        return (n == 1 || m == 1);
    }

    static constexpr bool is_square() {
        return (n == m);
    }

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

    /// Construct from a different type matrix
    // template <typename S, typename MT,
    //           std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    // Matrix_t(const Matrix_t<n, m, S, MT> &rhs) out_only {
    //     for (int i = 0; i < n * m; i++) {
    //         c[i] = rhs.c[i];
    //     }
    // }

    /// construct from 0
    inline Matrix_t(const std::nullptr_t &z) {
        for (int i = 0; i < n * m; i++) {
            c[i] = 0;
        }
    }

    // Construct matrix automatically from right-size initializer list
    // This does not seem to be dangerous, so keep non-explicit
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
    /**
     *
     */

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
    // RowVector<m, T> row(int r) const {
    //     RowVector<m, T> v;
    //     for (int i = 0; i < m; i++)
    //         v[i] = e(r, i);
    //     return v;
    // }

    /// return reference to row in a matrix
    const RowVector<m, T> &row(int r) const {
        return *(reinterpret_cast<const RowVector<m, T> *>(this) + r);
    }

    /// return reference to row in a matrix
    RowVector<m, T> &row(int r) {
        return *(reinterpret_cast<RowVector<m, T> *>(this) + r);
    }

    /// return reference to row in a matrix
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    void set_row(int r, const RowVector<m, S> &v) {
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
    void set_diagonal(const Vector<n, S> &v) {
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

    template <typename S>
    bool operator==(const Matrix<n, m, S> &rhs) const {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) {
                if (e(i, j) != rhs.e(i, j))
                    return false;
            }
        return true;
    }

    template <typename S>
    bool operator!=(const Matrix<n, m, S> &rhs) const {
        return !(*this == rhs);
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
    // #pragma hila loop_function
    template <typename S, typename MT,
              std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    Mtype &operator-=(const Matrix_t<n, m, S, MT> &rhs) {
        for (int i = 0; i < n * m; i++) {
            c[i] -= rhs.c[i];
        }
        return *this;
    }

    /// add assign a scalar to square matrix
    // #pragma hila loop_function
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
    // #pragma hila loop_function
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
    // #pragma hila loop_function
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
    inline const RowVector<n, T> &transpose() const {
        return *reinterpret_cast<const RowVector<n, T> *>(this);
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
        Matrix<n, m, hila::scalar_type<T>> res;
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
    inline Matrix<n, m, hila::scalar_type<T>> real() const {
        Matrix<n, m, hila::scalar_type<T>> res;
        for (int i = 0; i < m * n; i++) {
            res.c[i] = ::real(c[i]);
        }
        return res;
    }

    /// return imaginary part
    inline Matrix<n, m, hila::scalar_type<T>> imag() const {
        Matrix<n, m, hila::scalar_type<T>> res;
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
    hila::scalar_type<T> squarenorm() const {
        hila::scalar_type<T> result(0);
        for (int i = 0; i < n * m; i++) {
            result += ::squarenorm(c[i]);
        }
        return result;
    }

    /// calculate vector norm - sqrt of squarenorm
    template <typename S = T,
              std::enable_if_t<hila::is_floating_point<hila::scalar_type<S>>::value, int> = 0>
    hila::scalar_type<T> norm() const {
        return sqrt(squarenorm());
    }

    template <typename S = T,
              std::enable_if_t<!hila::is_floating_point<hila::scalar_type<S>>::value, int> = 0>
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

        static_assert(hila::is_floating_point<hila::scalar_type<T>>::value,
                      "Matrix/Vector random() requires non-integral type elements");

        for (int i = 0; i < n * m; i++) {
            hila::random(c[i]);
        }
        return *this;
    }

    /// Generate gaussian random elements
    Mtype &gaussian_random(double width = 1.0) out_only {

        static_assert(hila::is_floating_point<hila::scalar_type<T>>::value,
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
            if constexpr ((n * m) % 2 > 0) {
                c[n * m - 1] = hila::gaussrand() * width;
            }
        }
        return *this;
    }

    /// Reordering utilities: swap rows or columns by the permutation vector
    /// Permutation vector must be valid permutation of cols/rows
    Mtype reorder_columns(const Vector<m, int> &permutation) const {
        Mtype res;
        for (int i = 0; i < m; i++)
            res.set_column(i, this->column(permutation[i]));
        return res;
    }

    Mtype reorder_rows(const Vector<n, int> &permutation) const {
        Mtype res;
        for (int i = 0; i < n; i++)
            res.set_row(i, this->row(permutation[i]));
        return res;
    }

    /// implement also bare reorder for vectors
    template <int N>
    Mtype reorder(const Vector<N, int> &permutation) const {
        static_assert(
            n == 1 || m == 1,
            "reorder() only for vectors, use reorder_rows() or reorder_columns() for matrices");
        static_assert(N == Mtype::size(), "Incorrect size of permutation vector");

        Mtype res;
        for (int i = 0; i < N; i++) {
            res[i] = (*this)[permutation[i]];
        }
        return res;
    }

    /// Sort a real-valued Vector
    /// Two interfaces: first returns permutation vector, which can be used to reorder other
    /// vectors/matrices second does only sort

#pragma hila novector
    template <int N>
    Mtype sort(Vector<N, int> &permutation, hila::sort order = hila::sort::ascending) const {

        static_assert(n == 1 || m == 1, "Sorting possible only for vectors");
        static_assert(hila::is_arithmetic<T>::value,
                      "Sorting possible only for arithmetic vector elements");
        static_assert(N == Mtype::size(), "Incorrect size of permutation vector");

        for (int i = 0; i < N; i++)
            permutation[i] = i;
        if (hila::sort::nonsorted == order) {
            return *this;
        }

        if (hila::sort::ascending == order) {
            for (int i = 1; i < N; i++) {
                for (int j = i; j > 0 && c[permutation[j - 1]] > c[permutation[j]]; j--)
                    hila::swap(permutation[j], permutation[j - 1]);
            }
        } else {
            for (int i = 1; i < N; i++) {
                for (int j = i; j > 0 && c[permutation[j - 1]] < c[permutation[j]]; j--)
                    hila::swap(permutation[j], permutation[j - 1]);
            }
        }

        return this->reorder(permutation);
    }

#pragma hila novector
    Mtype sort(hila::sort order = hila::sort::ascending) const {
        static_assert(n == 1 || m == 1, "Sorting possible only for vectors");

        Vector<Mtype::size(), int> permutation;
        return sort(permutation, order);
    }


    /// Multiply (nxm)-matrix from left by a matrix which is 1 except for 4 elements
    /// on rows/columns p,q.

    template <typename R, typename Mt>
    void mult_by_2x2_left(int p, int q, const Matrix_t<2, 2, R, Mt> &M) {
        static_assert(hila::contains_complex<T>::value || !hila::contains_complex<R>::value,
                      "Cannot multiply real matrix with complex 2x2 block matrix");

        Vector<2, T> a;
        for (int i = 0; i < m; i++) {
            a.e(0) = this->e(p, i);
            a.e(1) = this->e(q, i);

            a = M * a;

            this->e(p, i) = a.e(0);
            this->e(q, i) = a.e(1);
        }
    }

    template <typename R, typename Mt>
    void mult_by_2x2_right(int p, int q, const Matrix_t<2, 2, R, Mt> &M) {
        static_assert(hila::contains_complex<T>::value || !hila::contains_complex<R>::value,
                      "Cannot multiply real matrix with complex 2x2 block matrix");

        RowVector<2, T> a;
        for (int i = 0; i < n; i++) {
            a.e(0) = this->e(i, p);
            a.e(1) = this->e(i, q);

            a = a * M;

            this->e(i, p) = a.e(0);
            this->e(i, q) = a.e(1);
        }
    }

    /// Calculate eigenvalues and -vectors of an hermitean (or symmetric) matrix.
    /// Returns the number of Jacobi iterations, or -1 if did not converge.
    /// Arguments will contain eigenvalues and eigenvectors in columns of the "eigenvectors" matrix
    /// Computation is done in double precision despite the input matrix types

#pragma hila novector
    template <typename Et, typename Mt, typename MT,
              typename Dtype = typename std::conditional<hila::contains_complex<T>::value,
                                                         Complex<double>, double>::type>
    int eigen_jacobi(out_only Vector<n, Et> &eigenvaluevec,
                     out_only Matrix_t<n, n, Mt, MT> &eigenvectors,
                     enum hila::sort sorted = hila::sort::nonsorted) const {

        static_assert(!hila::contains_complex<T>::value || hila::contains_complex<Mt>::value,
                      "Eigenvector matrix must be complex with complex original matrix");

        static_assert(n == m, "Eigensystem can be solved only for square matrices");
        int rot;
        SquareMatrix<n, Dtype> M, V;
        Vector<n, double> eigenvalues;

        // Do it in double prec; copy fields
        V = 1;
        M = (*this);

        // don't need the imag. parts of diag (are zero)
        eigenvalues = M.diagonal().real();

        for (rot = 0; rot < 300; rot++) {

            /* find the largest element */
            int p, q;
            double abs_mpq = -1.0; // this quarantees that even trivial matrix works
            for (int i = 0; i < n - 1; i++) {
                for (int j = i + 1; j < n; j++) {
                    double t = ::squarenorm(M.e(i, j));
                    if (abs_mpq < t) {
                        abs_mpq = t;
                        p = i;
                        q = j;
                    }
                }
            }

            abs_mpq = ::sqrt(abs_mpq);

            // if off-diag elements are tiny return

            if (abs_mpq <= 1e-18 * (std::abs(eigenvalues[p]) + std::abs(eigenvalues[q]))) {
                if (sorted == hila::sort::nonsorted) {

                    // return values and vectors as is
                    eigenvaluevec = eigenvalues;
                    eigenvectors = V;

                } else {
                    // bubble sort eigenvalues to decreasing order
                    Vector<n, int> perm;
                    eigenvaluevec = eigenvalues.sort(perm, sorted);
                    eigenvectors = V.reorder_columns(perm);
                }
                return (rot);
            }

            Dtype mpq = M.e(p, q);

            /* now do the p,q-rotation:
             * M <- P^+ M P, where
             *     | c   s    |  p
             * P = |-s*  c*   |  q  = Ppq
             *     |        1 |  r
             *
             * with P^+ = P^-1 ->  cc* + ss* = 1, det P = 1
             * Thus, P is SU(2) or O(2) matrix
             *
             *     | c*  -s   | | mpp mpq mpr| | c   s   |
             * M = | s*   c   |	| mqp mqq mqr| |-s*  c*  |
             *     |         1|	| mrp mrq mrr| |        1|
             *
             *   = Pip* Mij Pjq,  mqp = mpq*
             *
             * Set now Mpq = (c*mpp - s mqp)s + (c*mpq -smqq)c* = 0
             *             = c*^2 mpq - s^2 mpq* + c*s (mpp - mqq)
             *             = |mpq| [ c*^2 e - s^2 e* + c*s (mpp-mqq)/|mpq| ]
             * where e = exp(i arg(mpq)) = mpq/|mpq|, e* = 1/e
             *
             * Now the "rotation angle" (c~cos\phi, s~sin\phi)
             * a = "cot2\phi" = c*^2 e - s^2 e* / 2c*s = (mqq - mpp) / 2|mpq|
             * Def t = s/c*e, so that the above is
             * t^2 + 2ta - 1 = 0 ->
             * t = -a +- sqrt(a^2+1).  Choose one with smaller |t|!
             * This is prone to cancellation, ~ a - a, so write it as
             * t = sgn(a)/(|a| + sqrt(a^2 + 1).
             * Now can choose real c*
             * c* = 1/sqrt(t^2+1)
             * s = t c* e   (and c*c + s*s = 1)
             */

            double a = (eigenvalues[q] - eigenvalues[p]) / (2 * abs_mpq);
            double t = 1.0 / (std::abs(a) + std::sqrt(a * a + 1.0));
            if (a < 0.0)
                t = -t;
            double c = 1.0 / std::sqrt(t * t + 1.0);

            Dtype s = mpq * (t * c / abs_mpq);
            Matrix<2, 2, Dtype> P;
            P.e(0, 0) = P.e(1, 1) = 1.0 / std::sqrt(t * t + 1.0);
            P.e(0, 1) = s;
            P.e(1, 0) = -::conj(s);

            M.mult_by_2x2_left(p, q, P.dagger());
            M.mult_by_2x2_right(p, q, P);

            eigenvalues[p] = ::real(M.e(p, p));
            eigenvalues[q] = ::real(M.e(q, q));

            // p,q -elements of m should be 0 - set explictly to avoid rounding erros
            M.e(p, q) = 0;
            M.e(q, p) = 0;

            /* Now M done, take care of the ev's too ..
             * V' = V P = |vpp vpq vpr| | c  s   | = V_ik P_kj
             *            |vqp vqq vqr| |-s* c   |
             * 	          |vrp vrq vrr| |       1|
             * vip <- vip c - viq s*
             * viq <- vip s + viq c
             * vir <- vir
             */

            V.mult_by_2x2_right(p, q, P);
        }

        return (-1);
    }

    /// Convert to string for printing
    ///
    std::string str(int prec = 8, char separator = ' ') const {
        std::stringstream text;
        text.precision(prec);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                text << e(i, j);
                if (i < n - 1 || j < m - 1)
                    text << separator;
            }
        }
        return text.str();
    }
};

/**
 * @brief \f$ n \times m \f$ Matrix class which defines matrix operations.
 * @details The Matrix class is a derived class of Matrix_t, which is the general definition of a
 * matrix class. See Matrix_t details section for reasoning.
 *
 * All mathematical operations for Matrix are inherited from Matrix_t and are visible on this page.
 *
 * To see all methods of initializing a matrix see constructor method #Matrix::Matrix
 *
 * @tparam n row length
 * @tparam m column length
 * @tparam T Data type Matrix
 */
template <int n, int m, typename T>
class Matrix : public Matrix_t<n, m, T, Matrix<n, m, T>> {

  public:
    /// std incantation for field types
    using base_type = hila::scalar_type<T>;
    using argument_type = T;

    /**
     * @brief Construct a new Matrix object
     * @details The following ways of initializing a matrix are:
     *
     * __default constructor__:
     * \code {.cpp}
     * .
     * .
     * .
     * Matrix<n,m,MyType> M;
     * \endcode
     *
     *
     */
    Matrix() = default;
    ~Matrix() = default;
    Matrix(const Matrix &v) = default;


    /// constructor from scalar -- keep it explicit!  Not good for auto use
    template <typename S, int nn = n, int mm = m,
              std::enable_if_t<(hila::is_assignable<T &, S>::value && nn == mm), int> = 0>
    explicit Matrix(const S rhs) out_only {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                if (i == j)
                    this->e(i, j) = rhs;
                else
                    this->e(i, j) = 0;
            }
    }

    /// Construct from a different type matrix
    template <typename S, typename MT,
              std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    Matrix(const Matrix_t<n, m, S, MT> &rhs) out_only {
        for (int i = 0; i < n * m; i++) {
            this->c[i] = rhs.c[i];
        }
    }

    /// construct from 0
    Matrix(const std::nullptr_t &z) out_only {
        for (int i = 0; i < n * m; i++) {
            this->c[i] = 0;
        }
    }

    /// Construct matrix automatically from right-size initializer list
    /// This does not seem to be dangerous, so keep non-explicit

    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    Matrix(std::initializer_list<S> rhs) {
        assert(rhs.size() == n * m &&
               "Matrix/Vector initializer list size must match variable size");
        int i = 0;
        for (auto it = rhs.begin(); it != rhs.end(); it++, i++) {
            this->c[i] = *it;
        }
    }


    // use the Base::Base -trick to inherit assignments
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
//  Matrix<Complex<type_sum(scalar_type(Mt),scalar_type(S))>>
//  - otherwise return Matrix<type_sum>

template <typename Mt, typename S, typename Enable = void>
struct matrix_scalar_op_s {
    using type = Matrix<Mt::rows(), Mt::columns(),
                        Complex<hila::type_plus<hila::scalar_type<Mt>, hila::scalar_type<S>>>>;
};

template <typename Mt, typename S>
struct matrix_scalar_op_s<
    Mt, S,
    typename std::enable_if_t<std::is_convertible<hila::type_plus<hila::number_type<Mt>, S>,
                                                  hila::number_type<Mt>>::value>> {
    // using type = Mt;
    using type = typename std::conditional<
        hila::is_floating_point<hila::scalar_type<Mt>>::value, Mt,
        Matrix<Mt::rows(), Mt::columns(),
               hila::type_plus<hila::scalar_type<Mt>, hila::scalar_type<S>>>>::type;
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
    typename Rtype = Matrix<Mtype::rows() - 1, Mtype::columns() - 1, hila::number_type<Mtype>>>
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
          typename T = hila::number_type<Mtype>, int n = Mtype::rows(), int m = Mtype::columns(),
          std::enable_if_t<(n > 3), int> = 0>
T det(const Mtype &mat) {
    static_assert(n == m, "determinants defined only for square matrices");
    T result(0);
    hila::scalar_type<T> parity = 1, opposite = -1;
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
          typename R = hila::ntype_op<hila::number_type<Mt1>, hila::number_type<Mt2>>,
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
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++) {
            strm << A.e(i, j);
            if (i < n - 1 || j < m - 1)
                strm << ' ';
        }
    return strm;
}

/// Convert to string for "pretty" printing
///

namespace hila {


template <int n, int m, typename T, typename MT>
std::string to_string(const Matrix_t<n, m, T, MT> &A, int prec = 8, char separator = ' ') {
    std::stringstream strm;
    strm.precision(prec);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++) {
            strm << hila::to_string(A.e(i, j), prec, separator);
            if (i < n - 1 || j < m - 1)
                strm << separator;
        }
    return strm.str();
}

template <int n, int m, typename T, typename MT>
std::string prettyprint(const Matrix_t<n, m, T, MT> &A, int prec = 8) {
    std::stringstream strm;
    strm.precision(prec);

    if constexpr (n == 1) {
        // print a vector, horizontally
        strm << '[';
        for (int i = 0; i < n * m; i++)
            strm << ' ' << prettyprint(A.e(i), prec);
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
                item << prettyprint(A.e(i, j), prec);
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
        for (int i = 0; i < n - 1; i++) {
            strm << lines[i] << "]\n";
        }
        strm << lines[n - 1] << "]";
    }
    return strm.str();
}

} // namespace hila

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


/// find determinant using LU decomposition. Algorithm: numerical Recipes, 2nd ed.
/// p. 47 ff
template <int n, int m, typename T, typename MT, typename radix = hila::scalar_type<T>,
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
    hila::scalar_type<T> one = 1.0;

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
