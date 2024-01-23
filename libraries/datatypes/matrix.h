/**
 * @file matrix.h
 * @brief Definition of Matrix types
 * @details This file contains base matrix type Matrix_t which defines all general matrix type
 * operations Matrix types are Matrix, #Vector, #RowVector, #SquareMatrix of which Matrix is defined
 * as a class and the rest are special cases of the Matrix class.
 *
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <type_traits>
#include <sstream>
#include "plumbing/defs.h"
#include "datatypes/cmplx.h"

// forward definitions of needed classes
template <const int n, const int m, typename T, typename Mtype>
class Matrix_t;

template <const int n, const int m, typename T = double>
class Array;

template <int n, int m, typename T>
class Matrix;

template <int n, typename T>
class DiagonalMatrix;

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


/// @brief  type to store the return combo of svd:
///   {U, D, V} where U and V are nxn unitary / orthogonal,
/// and D is real diagonal singular value matrices.
/// @tparam M  - type of input matrix
template <typename M>
struct svd_result {
    static_assert(M::is_matrix() && M::rows() == M::columns(), "SVD only for square matrix");
    M U;
    DiagonalMatrix<M::size(), hila::arithmetic_type<M>> singularvalues;
    M V;
};

/// @brief  type to store the return value of eigen_hermitean():
///   {E, U} where E is nxn DiagonalMatrix containing eigenvalues and
/// U nxn unitary matrix, with eigenvector columns
/// @tparam M  - type of input matrix
template <typename M>
struct eigen_result {
    static_assert(M::is_matrix() && M::rows() == M::columns(),
                  "Eigenvalues only for square matrix");
    DiagonalMatrix<M::size(), hila::arithmetic_type<M>> eigenvalues;
    M eigenvectors;
};

/**
 * @brief The main \f$ n \times m \f$ matrix type template Matrix_t. This is a base class type for
 * "useful" types which are derived from this.
 *
 * @details Uses curiously recurring template pattern (CRTP), where the last template parameter is
 * the template itself
 *
 * Example: the Matrix type below is defined as
 * @code{.cpp}
 * template <int n, int m, typename T>
 * class Matrix : public Matrix_t<n, m, T, Matrix<n, m, T>> { .. }
 * @endcode
 *
 * This pattern is used because stupid c++ makes it complicated to write generic code, in this case
 * derived functions to return derived type
 *
 * @tparam n Number of rows
 * @tparam m Number of columns
 * @tparam T Matrix element type
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
    using base_type = hila::arithmetic_type<T>;
    using argument_type = T;

    // help for templates, can use T::is_matrix()
    // Not very useful outside template parameters
    static constexpr bool is_matrix() {
        return true;
    }

    /**
     * @brief Returns true if Matrix is a vector
     *
     * @return true
     * @return false
     */
    static constexpr bool is_vector() {
        return (n == 1 || m == 1);
    }

    /**
     * @brief Returns true if matrix is a square matrix
     *
     * @return true
     * @return false
     */
    static constexpr bool is_square() {
        return (n == m);
    }

    /// Define default constructors to ensure std::is_trivial
    Matrix_t() = default;
    ~Matrix_t() = default;
    Matrix_t(const Matrix_t &v) = default;

    // constructor from scalar -- keep it explicit!  Not good for auto use
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

    // Construct from a different type matrix
    // template <typename S, typename MT,
    //           std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    // Matrix_t(const Matrix_t<n, m, S, MT> &rhs) out_only {
    //     for (int i = 0; i < n * m; i++) {
    //         c[i] = rhs.c[i];
    //     }
    // }

    // construct from 0
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
    // #pragma hila loop_function
    template <typename Tm = Mtype,
              std::enable_if_t<!std::is_same<Tm, Matrix<n, m, T>>::value, int> = 0>
    inline operator Mtype &() {
        return *reinterpret_cast<Mtype *>(this);
    }

    // #pragma hila loop_function
    template <typename Tm = Mtype,
              std::enable_if_t<!std::is_same<Tm, Matrix<n, m, T>>::value, int> = 0>
    inline operator const Mtype &() const {
        return *reinterpret_cast<const Mtype *>(this);
    }

    // automatically cast to generic matrix
    // #pragma hila loop_function
    inline operator Matrix<n, m, T> &() {
        return *reinterpret_cast<Matrix<n, m, T> *>(this);
    }
    // #pragma hila loop_function
    inline operator const Matrix<n, m, T> &() const {
        return *reinterpret_cast<const Matrix<n, m, T> *>(this);
    }

    /// Define constant methods rows(), columns() - may be useful in template code

    /**
     * @brief Returns row length
     *
     * @return constexpr int
     */
    static constexpr int rows() {
        return n;
    }
    /**
     * @brief Returns column length
     *
     * @return constexpr int
     */
    static constexpr int columns() {
        return m;
    }


    /**
     * @brief Returns size of #Vector or square Matrix
     *
     * @tparam q row size n
     * @tparam p column size m
     * @return constexpr int
     */
    // size for row vector
    template <int q = n, int p = m, std::enable_if_t<q == 1, int> = 0>
    static constexpr int size() {
        return p;
    }
    // size for column vector
    template <int q = n, int p = m, std::enable_if_t<p == 1, int> = 0>
    static constexpr int size() {
        return q;
    }
    // size for square matrix
    template <int q = n, int p = m, std::enable_if_t<q == p, int> = 0>
    static constexpr int size() {
        return q;
    }

    /**
     * @brief Standard array indexing operation for matrices
     * @details Accessing singular elements is insufficient, but Matrix elements are often quite
     * small.
     *
     * Exammple for matrix:
     * \code
     *  Matrix<n,m,MyType> M;
     *  MyType a = M.e(i,j); \\ i <= n, j <= m
     * \endcode
     *
     *
     * @param i Row index
     * @param j Column index
     * @return T matrix element type
     */
    inline T e(const int i, const int j) const {
        // return elem[i][j];
        return c[i * m + j];
    }

    /// @internal const_function implementation. See const_function for details
    inline T &e(const int i, const int j) const_function {
        // return elem[i][j];
        return c[i * m + j];
    }
    // declare single e here too in case we have a vector
    // (n || m == 1)
    /**
     * @brief Standard array indexing operation for vectors
     * @details  Accessing singular elements is insufficient, but Vector elements are often quite
     * small.
     * \code {.cpp}
     *  Vector<n,MyType> V;
     * MyType a = V.e(i) \\ i <= n
     * \endcode
     * @tparam q Number of rows
     * @tparam p Number of columns
     * @param i Index
     * @return T
     */
    template <int q = n, int p = m, std::enable_if_t<(q == 1 || p == 1), int> = 0>
    inline T e(const int i) const {
        return c[i];
    }

    /// @internal const_function implementation. See const_function for details
    template <int q = n, int p = m, std::enable_if_t<(q == 1 || p == 1), int> = 0>
    inline T &e(const int i) const_function {
        return c[i];
    }

    /**
     * @brief Indexing operation [] defined only for vectors.
     *
     * @details Example:
     *
     * \code {.cpp}
     * Vector<n,MyType> V;
     * MyType a = V[i] \\ i <= n
     * \endcode
     *
     * @tparam q row size n
     * @tparam p column size m
     * @param i row or vector index depending on which is being indexed
     * @return T
     */
    template <int q = n, int p = m, std::enable_if_t<(q == 1 || p == 1), int> = 0>
    inline T operator[](const int i) const {
        return c[i];
    }
    /// @internal const_function implementation. See const_function for details
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


    /**
     * @brief Return reference to row in a matrix
     * @details Since the Matrix data is ordered in a row major order accessing a row returns a
     * reference to the row.
     *
     * \code{.cpp}
     * Matrix<n,n,T> M;
     * M.random();
     * RowVector<n,T> V = M.row(i);
     * \endcode
     *
     * @param r index of row to be referenced
     * @return const RowVector<m, T>&
     */
    const RowVector<m, T> &row(int r) const {
        return *(reinterpret_cast<const RowVector<m, T> *>(this) + r);
    }

    /**
     * @brief Set row of Matrix with #RowVector if types are assignable
     * @details
     * \code{.cpp}
     * RowVector<n,T> V;
     * V.random();
     * Matrix<n,n,T> M;
     * M.set_row(i,V);
     * \endcode
     * @tparam S RowVector type
     * @param r Index of row to be set
     * @param v RowVector to be set
     */
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    void set_row(int r, const RowVector<m, S> &v) {
        for (int i = 0; i < m; i++)
            e(r, i) = v[i];
    }

    /**
     * @brief Returns column vector as value at index c
     * \code{.cpp}
     * Matrix<n,n,T> M;
     * M.random();
     * Vector<n,T> V = M.column(i);
     * \endcode
     * @param c index of column vector to be returned
     * @return const Vector<n, T>
     */
    Vector<n, T> column(int c) const {
        Vector<n, T> v;
        for (int i = 0; i < n; i++)
            v[i] = e(i, c);
        return v;
    }

    /// get column of a matrix
    // hila_matrix_column_t<n, T, Mtype> column(int c) {
    //     return hila_matrix_column_t<n, T, Mtype>(*this, c);
    // }

    /**
     * @brief Set column of Matrix with #Vector if types are assignable
     * @details
     * \code{.cpp}
     * Vector<n,T> V;
     * V.random();
     * Matrix<n,n,T> M;
     * M.set_column(i,V);
     * \endcode
     * @tparam S Vector type
     * @param c Index of column to be set
     * @param v #Vector to be set
     */
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    void set_column(int c, const Vector<n, S> &v) {
        for (int i = 0; i < n; i++)
            e(i, c) = v[i];
    }

    /**
     * @brief Return diagonal of square matrix
     * @details If called for non square matrix the program will throw an error.
     *
     * \code{.cpp}
     * Matrix<n,n,T> M = 1;
     * DiagonalMatrix<n,T> D = M.diagonal();
     * \endcode
     * @return Vector<n, T> returned vector.
     */
    DiagonalMatrix<n, T> diagonal() {
        static_assert(n == m, "diagonal() method defined only for square matrices");
        DiagonalMatrix<n, T> res;
        for (int i = 0; i < n; i++)
            res.e(i) = (*this).e(i, i);
        return res;
    }

    /**
     * @brief Set diagonal of square matrix to #Vector which is passed to the method
     * @details If called for non square matrix the program will throw an error.
     *
     * Example:
     *
     * \code {.cpp}
     * SquareMatrix<n,MyType> S = 0; \\ Zero matrix
     * Vector<n,MyType> V = 1; \\ Vector assigned to 1 at all elements
     * S.set_diagonal(V); \\ Results in Identity matrix of size n
     * \endcode
     *
     * @tparam S type vector to assign values to
     * @param v Vector to assign to diagonal
     */
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    void set_diagonal(const Vector<n, S> &v) {
        static_assert(n == m, "set_diagonal() method defined only for square matrices");
        for (int i = 0; i < n; i++)
            (*this).e(i, i) = v.e(i);
    }


    /**
     * @brief Cast Matrix to Array
     * @details used for array operations
     *
     * @return Array<n, m, T>&
     */
    const Array<n, m, T> &asArray() const {
        return *reinterpret_cast<const Array<n, m, T> *>(this);
    }
    // Same as above but with const_function, see const_function for details
    Array<n, m, T> &asArray() const_function {
        return *reinterpret_cast<Array<n, m, T> *>(this);
    }

    /**
     * @brief Cast Vector to DiagonalMatrix
     *
     * @return DiagonalMatrix<n,T>
     */
    template <int mm = m, std::enable_if_t<mm == 1, int> = 0>
    const DiagonalMatrix<n, T> &asDiagonalMatrix() const {
        return *reinterpret_cast<const DiagonalMatrix<n, T> *>(this);
    }
    // Same as above but with const_function, see const_function for details
    template <int mm = m, std::enable_if_t<mm == 1, int> = 0>
    DiagonalMatrix<n, T> &asDiagonalMatrix() const_function {
        return *reinterpret_cast<DiagonalMatrix<n, T> *>(this);
    }


    // casting from one Matrix (number) type to another: do not do this automatically.
    // but require an explicit cast operator.  This makes it easier to write code.
    // or should it be automatic?  keep/remove explicit?
    // TODO: CHECK AVX CONVERSIONS

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

    /**
     * @brief Unary operator
     * @details
     *
     * Negation operation
     *
     * \code {.cpp}
     * M == -M;
     * \endcode
     *
     */
    inline Mtype operator-() const {
        Mtype res;
        for (int i = 0; i < n * m; i++) {
            res.c[i] = -c[i];
        }
        return res;
    }

    /**
     * @brief Addition operator
     * @memberof Matrix_t
     * @details
     *
     * Identity operation
     *
     * \code {.cpp}
     * M == +M;
     * \endcode
     *
     * @tparam Mtype1 Matrix type for a
     * @tparam Mtype2 Matrix type for b
     * @param a Left matrix
     * @param b Right matrix
     * @return Rtype
     */
    inline const Mtype &operator+() const {
        return *this;
    }

    /**
     * @brief Boolean operator == to determine if two matrices are exactly the same
     * @details Tolerance for equivalence is zero, meaning that the matrices must be exactly the
     * same.
     *
     * @tparam S Type for Matrix which is being compared to
     * @param rhs right hand side Matrix which we are comparing
     * @return true
     * @return false
     */
    template <typename S>
    bool operator==(const Matrix<n, m, S> &rhs) const {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) {
                if (e(i, j) != rhs.e(i, j))
                    return false;
            }
        return true;
    }

    /**
     * @brief Boolean operator != to check if matrices are exactly different
     * @details if matrices are exactly the same then this will return false
     *
     * @tparam S Type for MAtrix which is being compared to
     * @param rhs right hand side Matrix which we are comparing
     * @return true
     * @return false
     */
    template <typename S>
    bool operator!=(const Matrix<n, m, S> &rhs) const {
        return !(*this == rhs);
    }

    /**
     * @brief Copy matrix assignment
     * @details
     * \code {.cpp}
     * Matrix<n,m,MyType> M_0;
     * .
     * . M_0 has values assigned to it
     * .
     * Matrix<n,m,MyType> M; \\ undefined matrix
     * M = M_0; \\ Assignment from M_0
     * \endcode
     * @param rhs Matrix to assign
     * @return
     */
    // #pragma hila loop_function
    template <typename S, typename MT,
              std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    inline Mtype &operator=(const Matrix_t<n, m, S, MT> &rhs) out_only {
        for (int i = 0; i < n * m; i++) {
            c[i] = rhs.c[i];
        }
        return *this;
    }

    // assign from 0
    // #pragma hila loop_function
    /**
     * @brief Zero assignment
     *
     * \code {.cpp}
     * Matrix<n,m,MyType> M;
     * M = 0; Zero matrix;
     * \endcode
     *
     * @param z 0
     * @return Mtype&
     */
    inline Mtype &operator=(const std::nullptr_t &z) out_only {
        for (int i = 0; i < n * m; i++) {
            c[i] = 0;
        }
        return *this;
    }

    // Assign from "scalar" for square matrix
    // #pragma hila loop_function
    /**
     * @brief Assignment from scalar
     * @details Assigns the scalar to the diagonal elements as \f$ M = I\cdot a\f$
     *
     * \code {.cpp}
     * MyType a = hila::random();
     * Matrix<n,m,MyType> M;
     * M = a; M = I*a
     * \endcode

     * @tparam S Scalar type to assign
     * @param rhs Scalar to assign
     * @return Mtype&
     */
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

    // Assign from diagonal matrix
    // #pragma hila loop_function

    /**
     * @brief Assignment from diagonal matrix
     *
     * @tparam S Element type of Diagonal matrix
     * @param rhs Diagonal matrix to assign
     * @return Mtype&
     */
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    inline Mtype &operator=(const DiagonalMatrix<n, S> &rhs) out_only {
        static_assert(n == m,
                      "Assigning DiagonalMatrix to Matrix possible only for square matrices");

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                if (i == j)
                    e(i, j) = rhs.e(i);
                else
                    e(i, j) = 0;
            }
        return *this;
    }

    /**
     * @brief Initializer list assignment
     * @details
     * \code{.cpp}
     * Matrix<2,2,int> M ;
     * M = {1, 0
     *      0, 1};
     * \endcode
     * @tparam S Element type of initializer list
     * @param rhs Initializer list to assign
     * @return Mtype&
     */
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    inline Mtype &operator=(std::initializer_list<S> rhs) out_only {
        assert(rhs.size() == n * m && "Initializer list has a wrong size in assignment");
        int i = 0;
        for (auto it = rhs.begin(); it != rhs.end(); it++, i++) {
            c[i] = *it;
        }
        return *this;
    }


    /**
     * @brief Add assign operator Matrix to Matrix
     * @details
     *
     * \code {.cpp}
     * Matrix<n,m,MyType> M,N;
     * M = 1;
     * N = 1;
     * M += N; \\M = 2*I
     * \endcode
     * @tparam S Element type of rhs
     * @tparam MT Matrix type of rhs
     * @param rhs Matrix to add
     * @return Mtype&
     */
    // #pragma hila loop_function
    template <typename S, typename MT,
              std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    Mtype &operator+=(const Matrix_t<n, m, S, MT> &rhs) {
        for (int i = 0; i < n * m; i++) {
            c[i] += rhs.c[i];
        }
        return *this;
    }

    /**
     * @brief Subtract assign operator Matrix to MAtrix
     * @details
     *
     * \code {.cpp}
     * Matrix<n,m,MyType> M,N;
     * M = 3;
     * N = 1;
     * M -= N; \\M = 2*I
     * \endcode
     *
     * @param rhs Matrix to subtract with
     * @return Mtype&
     */
    // #pragma hila loop_function
    template <typename S, typename MT,
              std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    Mtype &operator-=(const Matrix_t<n, m, S, MT> &rhs) {
        for (int i = 0; i < n * m; i++) {
            c[i] -= rhs.c[i];
        }
        return *this;
    }

    // add assign a scalar to square matrix
    // #pragma hila loop_function
    /**
     * @brief Add assign operator scalar to Matrix
     * @details
     * Adds scalar \f$ a \f$ to __square__ matrix as \f$ M + a\cdot\mathbb{1} \f$
     *
     * \code {.cpp}
     * Matrix<n,m,MyType> M = 1;
     * M += 1 ; \\ M = 2*I
     * \endcode
     *
     * @tparam S Element type of scalar rhs
     * @param rhs Sclar to add
     * @return Mtype&
     */
    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_plus<T, S>>::value, int> = 0>
    Mtype &operator+=(const S &rhs) {

        static_assert(n == m, "rows != columns : scalar addition possible for square matrix only!");

        for (int i = 0; i < n; i++) {
            e(i, i) += rhs;
        }
        return *this;
    }

    // subtract assign type T and convertible
    // #pragma hila loop_function
    /**
     * @brief Subtract assign operator scalar to Matrix
     * @details Subtract scalar \f$ a \f$ to __square__ matrix as \f$ M - a\cdot\mathbb{1} \f$
     *
     * \code {.cpp}
     * Matrix<n,m,MyType> M = 3;
     * M -= 1 ; \\ M = 2*I
     * \endcode
     * @tparam S Element type of scalar rhs
     * @param rhs scalar to add
     * @return Mtype&
     */
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

    /**
     * @brief Multiply assign scalar or matrix
     * @details Multiplication works as defined for matrix multiplication
     *
     * Matrix multiply assign only defined for square matrices, since the matrix dimensions would
     * change otherwise.
     *
     * \code {.cpp}
     * Matrix<n,m,MyType> M,N;
     * .
     * . Fill matrices M and N
     * .
     * M *= N; \\ M = M*N
     * \endcode
     *
     * @param rhs Matrix to multiply with
     * @return template <int p, typename S, typename MT,
     * std::enable_if_t<hila::is_assignable<T &, hila::type_mul<T, S>>::value, int> = 0>&
     */
    // #pragma hila loop_function
    template <int p, typename S, typename MT,
              std::enable_if_t<hila::is_assignable<T &, hila::type_mul<T, S>>::value, int> = 0>
    Mtype &operator*=(const Matrix_t<m, p, S, MT> &rhs) {
        static_assert(m == p, "can't assign result of *= to lhs Matrix, because doing so "
                              "would change it's dimensions");
        *this = *this * rhs;
        return *this;
    }

    // multiply assign with scalar
    // #pragma hila loop_function
    /**
     * @brief Multiply assign operator scalar to Matrix
     * @details Multiply Matrix uniformly with scalar
     *
     * \code {.cpp}
     * Matrix<n,m,MyType> M;
     * .
     * . Fill whole matrix with 1
     * .
     * M *= 2 ; \\ M is filled with 2
     * \endcode
     * @tparam S Scalar type of rhs
     * @param rhs Scalar to multiply
     * @return Mtype&
     */
    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_mul<T, S>>::value, int> = 0>
    Mtype &operator*=(const S rhs) {
        for (int i = 0; i < n * m; i++) {
            c[i] *= rhs;
        }
        return *this;
    }

    /**
     * @brief Divide assign oprator scalar with matrix
     * @details Divide works as defined for scalar matrix division.
     *
     *
     * \code {.cpp}
     * Matrix<n,m,MyType> M;
     * .
     * . Fill whole matrix with 2
     * .
     * M /= 2 ; \\ M is filled with 1
     * \endcode
     * @tparam S Scalar type of rhs
     * @param rhs Scalar to divide with
     * @return
     */
    // #pragma hila loop_function
    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_div<T, S>>::value, int> = 0>
    Mtype &operator/=(const S rhs) {
        for (int i = 0; i < n * m; i++) {
            c[i] /= rhs;
        }
        return *this;
    }

    /**
     * @brief Addition assign operator for DiagonalMatrix to Matrix
     * @details Works as long as Matrix to assign to is square
     * @tparam S Element type of rhs
     * @param rhs DiagonalMatrix to add
     * @return Mtype&
     */
    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_plus<T, S>>::value, int> = 0>
    Mtype &operator+=(const DiagonalMatrix<n, S> &rhs) {
        static_assert(n == m, "Assigning DiagonalMatrix possible only for square matrix");

        for (int i = 0; i < n; i++)
            e(i, i) += rhs.e(i);
        return *this;
    }

    /**
     * @brief Subtract assign operator for DiagonalMatrix to Matrix
     * @details Works as long as Matrix to assign to is square
     * @tparam S Element type of rhs
     * @param rhs DiagonalMatrix to add
     * @return Mtype&
     */
    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_plus<T, S>>::value, int> = 0>
    Mtype &operator-=(const DiagonalMatrix<n, S> &rhs) {
        static_assert(n == m, "Assigning DiagonalMatrix possible only for square matrix");

        for (int i = 0; i < n; i++)
            e(i, i) -= rhs.e(i);
        return *this;
    }

    /**
     * @brief Multiply assign operator for DiagonalMatrix to Matrix
     * @details Simply defined as matrix multiplication, but since rhs is guaranteed to be diagonal
     * the method is optimized to skip most of the elements.
     *
     * @tparam S Element type of rhs
     * @param rhs DiagonalMatrix to multiply
     * @return Mtype&
     */
    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_mul<T, S>>::value, int> = 0>
    Mtype &operator*=(const DiagonalMatrix<m, S> &rhs) {

        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                e(i, j) *= rhs.e(j);

        return *this;
    }

    /**
     * @brief Divide assign operator for DiagonalMatrix to Matrix
     * @details Well defined since rhs is guaranteed to be Diagonal.
     *
     * Let M be the Matrix which we divide and D be the DiagonalMatrix which we divide with then the
     * operation is defined as \f$ M \cdot D^{-1}  \f$.
     * @tparam S Element type of rhs
     * @param rhs DiagonalMatrix to divide
     * @return Mtype&
     */
    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_div<T, S>>::value, int> = 0>
    Mtype &operator/=(const DiagonalMatrix<m, S> &rhs) {

        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                e(i, j) /= rhs.e(j);

        return *this;
    }


    /**
     * @brief Matrix fill
     * @details Fills the matrix with element if it is assignable to matrix type T
     *
     * Works as follows:
     *
     * \code {.cpp}
     * Matrix<n,m,MyType> M;
     * M.fill(2) \\ Matrix is filled with 2
     * \endcode
     *
     * @tparam S Scalar type to of rhs
     * @param rhs Element to fill matrix with
     * @return const Mtype&
     */
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    const Mtype &fill(const S rhs) out_only {
        for (int i = 0; i < n * m; i++)
            c[i] = rhs;
        return *this;
    }

    /**
     * @brief Transpose of matrix
     * @details Return type for square matrix is same input, for non square return type is
     * Matrix<n,m,MyType>
     * @tparam mm
     * @tparam std::conditional<n, Mtype, Matrix<m, n, T>>::type
     * @return const Rtype
     */
    template <int mm = m,
              typename Rtype = typename std::conditional<n == m, Mtype, Matrix<m, n, T>>::type,
              std::enable_if_t<(mm != 1), int> = 0>
    inline Rtype transpose() const {
        Rtype res;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) {
                res.e(j, i) = e(i, j);
            }
        return res;
    }

    /**
     * @brief Transpose of vector
     * @details Returns reference
     *
     * @tparam mm
     * @return const RowVector<n, T>&
     */
    template <int mm = m, std::enable_if_t<mm == 1, int> = 0>
    inline const RowVector<n, T> &transpose() const {
        return *reinterpret_cast<const RowVector<n, T> *>(this);
    }

    /**
     * @brief Transpose of RowVector
     * @details Returns reference
     *
     * @tparam mm
     * @return const RowVector<n, T>&
     */
    template <int nn = n, std::enable_if_t<nn == 1 && m != 1, int> = 0>
    inline const Vector<n, T> &transpose() const {
        return *reinterpret_cast<const Vector<n, T> *>(this);
    }


    /**
     * @brief Returns complex conjugate of Matrix
     *
     * @return const Mtype
     */
    inline Mtype conj() const {
        Mtype res;
        for (int i = 0; i < n * m; i++) {
            res.c[i] = ::conj(c[i]);
        }
        return res;
    }

    /**
     * @brief Hermitian conjugate of matrix
     * @details for square matrix return type is same, for non square it is Matrix<m,n,MyType>
     *
     * @tparam std::conditional<n, Mtype, Matrix<m, n, T>>::type
     * @return const Rtype
     */
    template <typename Rtype = typename std::conditional<n == m, Mtype, Matrix<m, n, T>>::type>
    inline Rtype dagger() const {
        Rtype res;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) {
                res.e(j, i) = ::conj(e(i, j));
            }
        return res;
    }

    /**
     * @brief Adjoint of matrix
     * @details Alias to dagger
     *
     * @tparam std::conditional<n, Mtype, Matrix<m, n, T>>::type
     * @return Rtype
     */
    template <typename Rtype = typename std::conditional<n == m, Mtype, Matrix<m, n, T>>::type>
    inline Rtype adjoint() const {
        return dagger();
    }

    /**
     * @brief Returns absolute value of Matrix
     * @details For Matrix<n,m,Complex<T>> case type is changed to Matrix<n,m,T> as expected since
     * absolute value of a complex number is real
     *
     * @tparam M
     * @return Mtype
     */
    auto abs() const {
        Matrix<n, m, hila::arithmetic_type<T>> res;
        for (int i = 0; i < n * m; i++) {
            res.c[i] = ::abs(c[i]);
        }
        return res;
    }


    // It seems that using special "Dagger" type makes the code slower!
    // Disable it now
    // inline const DaggerMatrix<m,n,T> & dagger() const {
    //     return *reinterpret_cast<const DaggerMatrix<m,n,T>*>(this);
    // }

    /**
     * @brief Returns real part of Matrix or #Vector
     *
     * @return Matrix<n, m, hila::arithmetic_type<T>>
     */
    inline Matrix<n, m, hila::arithmetic_type<T>> real() const {
        Matrix<n, m, hila::arithmetic_type<T>> res;
        for (int i = 0; i < m * n; i++) {
            res.c[i] = ::real(c[i]);
        }
        return res;
    }

    /**
     * @brief Returns imaginary part of Matrix or #Vector
     *
     * @return Matrix<n, m, hila::arithmetic_type<T>>
     */
    inline Matrix<n, m, hila::arithmetic_type<T>> imag() const {
        Matrix<n, m, hila::arithmetic_type<T>> res;
        for (int i = 0; i < m * n; i++) {
            res.c[i] = ::imag(c[i]);
        }
        return res;
    }

    /**
     * @brief Computes Trace for Matrix
     * @details Not define for non square matrices. Static assert will stop code from compiling if
     * executed for non square matrices.
     *
     * @return T
     */
    T trace() const {
        static_assert(n == m, "trace not defined for non square matrices");
        T result(0);
        for (int i = 0; i < n; i++) {
            result += e(i, i);
        }
        return result;
    }

    /**
     * @brief Multiply with given matrix and compute trace of result
     * @details Slightly cheaper operation than \code (M*N).trace() \endcode
     *
     * @tparam p
     * @tparam q
     * @tparam S
     * @tparam MT
     * @param rm
     * @return hila::type_mul<T, S>
     */
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

    /**
     * @brief Calculate square norm - sum of squared elements
     *
     * @return hila::arithmetic_type<T>
     */
    hila::arithmetic_type<T> squarenorm() const {
        hila::arithmetic_type<T> result(0);
        for (int i = 0; i < n * m; i++) {
            result += ::squarenorm(c[i]);
        }
        return result;
    }

    /**
     * @brief Calculate vector norm - sqrt of squarenorm
     *
     * @tparam S
     * @return hila::arithmetic_type<T>
     */
    template <typename S = T,
              std::enable_if_t<hila::is_floating_point<hila::arithmetic_type<S>>::value, int> = 0>
    hila::arithmetic_type<T> norm() const {
        return sqrt(squarenorm());
    }

    template <typename S = T,
              std::enable_if_t<!hila::is_floating_point<hila::arithmetic_type<S>>::value, int> = 0>
    double norm() const {
        return sqrt(static_cast<double>(squarenorm()));
    }

    /**
     * @brief Find max of Matrix only for arithmetic types
     */
    template <typename S = T, std::enable_if_t<hila::is_arithmetic<S>::value, int> = 0>
    T max() const {
        T res = c[0];
        for (int i = 1; i < n * m; i++) {
            if (res < c[i])
                res = c[i];
        }
        return res;
    }

    /**
     * @brief Find min of Matrix only for arithmetic types
     */
    template <typename S = T, std::enable_if_t<hila::is_arithmetic<S>::value, int> = 0>
    T min() const {
        T res = c[0];
        for (int i = 1; i < n * m; i++) {
            if (res > c[i])
                res = c[i];
        }
        return res;
    }


    // dot product - (*this).dagger() * rhs
    // could be done as well by writing the operation as above!
    /**
     * @brief Dot product
     * @details Only works between two #Vector objects
     *
     * \code {.cpp}
     * Vector<m,MyType> V,W;
     * .
     * .
     * .
     * V.dot(W); // valid operation
     * \endcode
     *
     * \code {.cpp}
     * RowVector<m,MyType> V,W;
     * .
     * .
     * .
     * V.dot(W); // not valid operation
     * \endcode
     *
     *
     * @tparam p Row length for rhs
     * @tparam q Column length for rhs
     * @tparam S Type for rhs
     * @tparam R Gives resulting type of lhs and rhs multiplication
     * @param rhs Vector to compute dot product with
     * @return R Value of dot product
     */
    template <int p, int q, typename S, typename R = hila::type_mul<T, S>>
    inline R dot(const Matrix<p, q, S> &rhs) const {
        static_assert(m == 1 && q == 1 && p == n,
                      "dot() product only for vectors of the same length");

        R r = 0;
        for (int i = 0; i < n; i++) {
            r += ::conj(c[i]) * rhs.e(i);
        }
        return r;
    }

    // outer product - (*this) * rhs.dagger(), sizes (n,1) and (p,1)
    // gives n * p matrix
    // could be done as well by the above operation!

    /**
     * @brief Outer product
     * @details Only works between two #Vector objects
     *
     * \code{.cpp}
     * Vector<n,MyType> V,W;
     * .
     * .
     * .
     * V.outer_product(W); \\ Valid operation
     * \endcode
     *
     * \code {.cpp}
     * RowVector<m,MyType> V,W;
     * .
     * .
     * .
     * V.outer_product(W); // not valid operation
     * \endcode
     *
     * @tparam p Row length for rhs
     * @tparam q Column length for rhs
     * @tparam S Element type for rhs
     * @tparam R Type between lhs and rhs multiplication
     * @param rhs Vector to compute outer product with
     * @return Matrix<n, p, R>
     */
    template <int p, int q, typename S, typename R = hila::type_mul<T, S>>
    inline Matrix<n, p, R> outer_product(const Matrix<p, q, S> &rhs) const {
        static_assert(m == 1 && q == 1, "outer_product() only for vectors");

        Matrix<n, p, R> res;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < p; j++) {
                res.e(i, j) = c[i] * ::conj(rhs.e(j));
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

    /**
     * @brief Fills Matrix with random elements from uniform distribution
     * @details
     * \code{.cpp}
     * Matrix<n,m,T> M;
     * M.random();
     * \endcode
     * @return Mtype&
     */
    Mtype &random() out_only {

        static_assert(hila::is_floating_point<hila::arithmetic_type<T>>::value,
                      "Matrix/Vector random() requires non-integral type elements");

        for (int i = 0; i < n * m; i++) {
            hila::random(c[i]);
        }
        return *this;
    }

    /**
     * @brief Fills Matrix with gaussian random elements from gaussian distribution
     * @details
     * \code {.cpp}
     * Matrix<n,m,T> M;
     * M.gaussian_random();
     * \endcode
     *
     * @param width
     * @return Mtype&
     */
    Mtype &gaussian_random(double width = 1.0) out_only {

        static_assert(hila::is_floating_point<hila::arithmetic_type<T>>::value,
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

    /**
     * @brief permute columns of Matrix
     * @details Reordering is done based off of permutation vector
     *
     * @param permutation Vector of integers to permute columns
     * @return Mtype
     */
    Mtype permute_columns(const Vector<m, int> &permutation) const {
        Mtype res;
        for (int i = 0; i < m; i++)
            res.set_column(i, this->column(permutation[i]));
        return res;
    }

    /**
     * @brief permute rows of Matrix
     * @details Reordering is done based off of permutation vector
     *
     * @param permutation Vector of integers to permute rows
     * @return Mtype
     */
    Mtype permute_rows(const Vector<n, int> &permutation) const {
        Mtype res;
        for (int i = 0; i < n; i++)
            res.set_row(i, this->row(permutation[i]));
        return res;
    }


    /**
     * @brief Permute elements of vector
     * @details Reordering is done based off of permutation vector
     *
     * @param permutation Vector of integers to permute rows
     * @return Mtype
     */
    template <int N>
    Mtype permute(const Vector<N, int> &permutation) const {
        static_assert(
            n == 1 || m == 1,
            "permute() only for vectors, use permute_rows() or permute_columns() for matrices");
        static_assert(N == Mtype::size(), "Incorrect size of permutation vector");

        Mtype res;
        for (int i = 0; i < N; i++) {
            res[i] = (*this)[permutation[i]];
        }
        return res;
    }

    /**
     * @brief Sort method for #Vector which returns permutation order
     * @details Order can be past as argument as either ascending (default) or descending
     *
     * __Ascending__
     * @code {.cpp}
     * Vector<n,MyType> V;
     * Vector<n,int> perm;
     * V.random();
     * V.sort(perm);
     * V.permute(perm);
     * @endcode
     *
     * __Descending__
     * @code {.cpp}
     * Vector<n,MyType> V;
     * Vector<n,int> perm;
     * V.random();
     * V.sort(perm,hila::sort::descending);
     * V.permute(perm);
     * @endcode
     * @tparam N
     * @param permutation
     * @param order
     * @return Mtype
     */
#pragma hila novector
    template <int N>
    Mtype sort(Vector<N, int> &permutation, hila::sort order = hila::sort::ascending) const {

        static_assert(n == 1 || m == 1, "Sorting possible only for vectors");
        static_assert(hila::is_arithmetic<T>::value,
                      "Sorting possible only for arithmetic vector elements");
        static_assert(N == Mtype::size(), "Incorrect size of permutation vector");

        for (int i = 0; i < N; i++)
            permutation[i] = i;
        if (hila::sort::unsorted == order) {
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

        return this->permute(permutation);
    }

    /**
     * @brief  Sort method for #Vector
     * @details Order can be past as argument as either ascending (default) or descending
     *
     * __Ascending__
     * @code {.cpp}
     * Vector<n,MyType> V;
     * V.random();
     * V.sort(); // V is sorted in ascending order
     * @endcode
     *
     * __Descending__
     * @code {.cpp}
     * V.sort(hila::sort::descending);
     * @endcode
     *
     * @param order Order to sort in
     * @return Mtype
     */
#pragma hila novector
    Mtype sort(hila::sort order = hila::sort::ascending) const {
        static_assert(n == 1 || m == 1, "Sorting possible only for vectors");

        Vector<Mtype::size(), int> permutation;
        return sort(permutation, order);
    }


    /**
     * @brief Multiply \f$ n \times m \f$-matrix from the left by  \f$ n \times m \f$ matrix defined
     * by  \f$ 2 \times 2 \f$ sub matrix
     * @details The multiplication is defined as follows, let \f$M\f$ as the \f$ 2 \times 2 \f$
     * input matrix and \f$B\f$ be `(this)` matrix, being the matrix stored in the object this
     * method is called for. Let \f$A = I\f$ be a \f$ n \times m \f$ unit matrix. We then set the
     * values of A to be: \f{align}{ A_{p,p} = M_{0,0}, \hspace{5px} A_{p,q} = M_{0,1}, \hspace{5px}
     * A_{q,p} = M_{1,0}, \hspace{5px} A_{q,q} = M_{1,1}. \f}
     *
     * Then the resulting matrix will be:
     *
     * \f{align}{ B = A \cdot B  \f}
     *
     * @tparam R Element type of M
     * @tparam Mt Matrix type of M
     * @param p First row and column
     * @param q Second row and column
     * @param M \f$ 2 \times 2\f$ Matrix to multiply with
     */
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

    /**
     * @brief  Multiply \f$ n \times m \f$-matrix from the right by  \f$ n \times m \f$ matrix
     * defined by  \f$ 2 \times 2 \f$ sub matrix
     * @details See Matrix::mult_by_2x2_left, only difference being that the multiplication is from
     * the right.
     *
     * @tparam R Element type of M
     * @tparam Mt Matrix type of M
     * @param p First row and column
     * @param q Second row and column
     * @param M \f$ 2 \times 2\f$ Matrix to multiply with
     */
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


    //   Following methods are defined in matrix_linalg.h
    //

    /**
     * @brief following calculate the determinant of a square matrix
     * det() is the generic interface, using laplace for small matrices and LU for large
     */

#pragma hila novector
    T det_lu() const;
#pragma hila novector
    T det_laplace() const;
#pragma hila novector
    T det() const;

    /**
     * @brief  Calculate eigenvalues and -vectors of an hermitean (or symmetric) matrix.
     * @details Returns the number of Jacobi iterations, or -1 if did not converge.
     * Arguments will contain eigenvalues and eigenvectors in columns of the "eigenvectors" matrix.
     * Computation is done in double precision despite the input matrix types.
     * @param eigenvaluevec Vector of computed eigenvalue
     * @param eigenvectors Eigenvectors as columns of \f$ n \times n \f$ Matrix
     *
     */

#pragma hila novector
    template <typename Et, typename Mt, typename MT>
    int eigen_hermitean(out_only DiagonalMatrix<n, Et> &eigenvalues,
                        out_only Matrix_t<n, n, Mt, MT> &eigenvectors,
                        enum hila::sort sorted = hila::sort::unsorted) const;

    eigen_result<Mtype> eigen_hermitean(enum hila::sort sorted = hila::sort::unsorted) const;

    // Temporary interface to eigenvalues, to be removed
    template <typename Et, typename Mt, typename MT>
    int eigen_jacobi(out_only Vector<n, Et> &eigenvaluevec,
                     out_only Matrix_t<n, n, Mt, MT> &eigenvectors,
                     enum hila::sort sorted = hila::sort::unsorted) const {

        DiagonalMatrix<n, Et> d;
        int i = eigen_hermitean(d, eigenvectors, sorted);
        eigenvaluevec = d.asArray().asVector();
        return i;
    }


#pragma hila novector
    template <typename Et, typename Mt, typename MT>
    int svd(out_only Matrix_t<n, n, Mt, MT> &_U, out_only DiagonalMatrix<n, Et> &_D,
            out_only Matrix_t<n, n, Mt, MT> &_V,
            enum hila::sort sorted = hila::sort::unsorted) const;

    svd_result<Mtype> svd(enum hila::sort sorted = hila::sort::unsorted) const;


#pragma hila novector
    template <typename Et, typename Mt, typename MT>
    int svd_pivot(out_only Matrix_t<n, n, Mt, MT> &_U, out_only DiagonalMatrix<n, Et> &_D,
                  out_only Matrix_t<n, n, Mt, MT> &_V,
                  enum hila::sort sorted = hila::sort::unsorted) const;

    svd_result<Mtype> svd_pivot(enum hila::sort sorted = hila::sort::unsorted) const;


    //////// matrix_linalg.h

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
 * __NOTE__: Some method documentation is not being displayed, for example the assignment operator
 * documentation is not being inherited. See Matrix_t::operator= for use.
 *
 * To see all methods of initializing a matrix see constructor method #Matrix::Matrix
 *
 * __NOTE__: In the documentation examples n,m are integers and MyType is a HILA [standard
 * type](@ref standard) or Complex.
 *
 * @tparam n Number of rows
 * @tparam m Number of columns
 * @tparam T Data type Matrix
 */
template <int n, int m, typename T>
class Matrix : public Matrix_t<n, m, T, Matrix<n, m, T>> {

  public:
    /// std incantation for field types
    using base_type = hila::arithmetic_type<T>;
    using argument_type = T;

    /**
     * @brief Construct a new Matrix object
     * @details
     * __NOTE__: n,m are integers and MyType is a HILA [standard type](@ref standard) or Complex.
     *
     * __Default constructor__:
     *
     * Allocates undefined \f$ n\times m\f$ Array.
     *
     * \code{.cpp}
     * Matrix<n,m,MyType> M;
     * \endcode
     */
    Matrix() = default;
    ~Matrix() = default;
    /// Default copy constructor
    Matrix(const Matrix &v) = default;


    // constructor from scalar -- keep it explicit!  Not good for auto use
    /**
     * @brief Scalar constructor
     * @details
     *
     * Construct with given scalar at diagonal elements \f$ M = \mathbf{I}\cdot x\f$. Matrix must be
     * square \f$n == m\f$
     *
     * \code{.cpp}
     * MyType x = hila::random();
     * Matrix<n,m,MyType> M = x;
     * \endcode

     * @tparam S Type for scalar
     * @tparam nn Number of rows
     * @tparam mm Numebr of columns
     * @param rhs Scalar element to assign
     */
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

    /**
     * @brief Copy constructor
     *
     * Construction from previously defined matrix if types are compatible. For example the code
     * below only works if assignment from MyOtherType to MyType is defined.
     *
     * \code{.cpp}
     * Matrix<n,m,MyOtherType> M_0;
     * //
     * // M_0 is filled with content
     * //
     * Matrix<n,m,MyType> M = M_0;
     * \endcode

     * @tparam S Element type for copied matrix
     * @tparam MT Matrix type for copied matrix
     * @param rhs Matrix to copy
     */
    template <typename S, typename MT,
              std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    Matrix(const Matrix_t<n, m, S, MT> &rhs) out_only {
        for (int i = 0; i < n * m; i++) {
            this->c[i] = rhs.c[i];
        }
    }

    /**
     * @brief std::nullptr_t constructo
     * @details Constructing from 0 sets the whole Matrix to zero
     * @param z 0
     */
    Matrix(const std::nullptr_t &z) out_only {
        for (int i = 0; i < n * m; i++) {
            this->c[i] = 0;
        }
    }

    // Construct matrix automatically from right-size initializer list
    // This does not seem to be dangerous, so keep non-explicit
    /**
     * @brief Initializer list constructor
     *
     * Construction from c++ initializer list.
     *
     * \code{.cpp}
     * Matrix<2,2,int> M = {1, 0
     *                      0, 1};
     * \endcode
     * @tparam S Element type for initializer list
     * @param rhs
     */
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

    using Matrix_t<n, m, T, Matrix<n, m, T>>::operator-;
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
    using type =
        Matrix<Mt::rows(), Mt::columns(),
               Complex<hila::type_plus<hila::arithmetic_type<Mt>, hila::arithmetic_type<S>>>>;
};

template <typename Mt, typename S>
struct matrix_scalar_op_s<
    Mt, S,
    typename std::enable_if_t<std::is_convertible<hila::type_plus<hila::number_type<Mt>, S>,
                                                  hila::number_type<Mt>>::value>> {
    // using type = Mt;
    using type = typename std::conditional<
        hila::is_floating_point<hila::arithmetic_type<Mt>>::value, Mt,
        Matrix<Mt::rows(), Mt::columns(),
               hila::type_plus<hila::arithmetic_type<Mt>, hila::arithmetic_type<S>>>>::type;
};

template <typename Mt, typename S>
using mat_scalar_type = typename matrix_scalar_op_s<Mt, S>::type;


} // namespace hila

//////////////////////////////////////////////////////////////////////////

/**
 * @brief Return transposed Matrix or #Vector
 * @details Wrapper around Matrix::transpose
 *
 * Can be called as:
 * \code {.cpp}
 * Matrix<n,m,MyType> M,M_transposed;
 * M_transposed = transpose(M);
 * \endcode
 *
 *
 * @tparam Mtype Matrix type
 * @param arg Object to be transposed
 * @return auto Mtype transposed matrix
 */
template <typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0>
inline auto transpose(const Mtype &arg) {
    return arg.transpose();
}

/**
 * @brief Return conjugate Matrix or #Vector
 * @details Wrapper around Matrix::conj
 *
 * Can be called as:
 * \code {.cpp}
 * Matrix<n,m,MyType> M,M_conjugate;
 * M_conjugate = conj(M);
 * \endcode
 *
 *
 * @tparam Mtype Matrix type
 * @param arg Object to be conjugated
 * @return auto Mtype conjugated matrix
 */
template <typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0>
inline auto conj(const Mtype &arg) {
    return arg.conj();
}

/**
 * @brief Return adjoint Matrix
 * @details Wrapper around Matrix::adjoint
 *
 * Can be called as:
 * \code {.cpp}
 * Matrix<n,m,MyType> M,M_adjoint;
 * M_conjugate = adjoint(M);
 * \endcode
 *
 *
 * @tparam Mtype Matrix type
 * @param arg Object to compute adjoint of
 * @return auto Mtype adjoint matrix
 */
template <typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0>
inline auto adjoint(const Mtype &arg) {
    return arg.adjoint();
}

/**
 * @brief Return dagger of Matrix
 * @details Wrapper around Matrix::adjoint
 *
 * Same as adjoint
 *
 *
 * @tparam Mtype Matrix type
 * @param arg Object to compute dagger of
 * @return auto Mtype dagger matrix
 */
template <typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0>
inline auto dagger(const Mtype &arg) {
    return arg.adjoint();
}

/**
 * @brief Return absolute value Matrix or #Vector
 * @details Wrapper around Matrix::abs
 *
 * Can be called as:
 * \code {.cpp}
 * Matrix<n,m,MyType> M,M_abs;
 * M_abs = abs(M);
 * \endcode
 *
 *
 * @tparam Mtype Matrix type
 * @param arg Object to compute absolute value of
 * @return auto Mtype absolute value Matrix or Vector
 */
template <typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0>
inline auto abs(const Mtype &arg) {
    return arg.abs();
}

/**
 * @brief Return trace of square Matrix
 * @details Wrapper around Matrix::trace
 *
 * Can be called as:
 * \code {.cpp}
 * Matrix<n,m,MyType> M;
 * MyType T;
 * T = trace(M);
 * \endcode
 *
 *
 * @tparam Mtype Matrix type
 * @param arg Object to compute trace of
 * @return auto Trace of Matrix
 */
template <typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0>
inline auto trace(const Mtype &arg) {
    return arg.trace();
}

/**
 * @brief Return real of Matrix or #Vector
 * @details Wrapper around Matrix::real
 *
 * Can be called as:
 * \code {.cpp}
 * Matrix<n,m,MyType> M,M_real;
 * M_real = real(M);
 * \endcode
 *
 *
 * @tparam Mtype Matrix type
 * @param arg Object to compute real part of
 * @return auto Mtype real part of arg
 */
template <typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0>
inline auto real(const Mtype &arg) {
    return arg.real();
}

/**
 * @brief Return imaginary of Matrix or #Vector
 * @details Wrapper around Matrix::imag
 *
 * Can be called as:
 * \code {.cpp}
 * Matrix<n,m,MyType> M,M_imag;
 * M_imag = imag(M);
 * \endcode
 *
 *
 * @tparam Mtype Matrix type
 * @param arg Object to compute imag part of
 * @return auto Mtype imag part of arg
 */
template <typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0>
inline auto imag(const Mtype &arg) {
    return arg.imag();
}


///////////////////////////////////////////////////////////////////////////
// Now matrix additions: matrix + matrix


/**
 * @brief Addition operator Matrix + Matrix
 * @memberof Matrix_t
 * @details Addition operator between matrices is defined in the usual way (element wise).
 *
 * __NOTE__: Matrices must share same dimensionality.
 *
 * \code {.cpp}
 * Matrix<n,m,MyType> M, N, S;
 * M.fill(1);
 * N.fill(1);
 * S = M + N; // Resulting matrix is uniformly 2
 * \endcode
 * @tparam Mtype1 Matrix type of a
 * @tparam Mtype2 Matrix type of b
 * @param a Matrix to add
 * @param b Matrix to add
 * @return Rtype Return matrix of compatible type between Mtype1 and Mtype2
 */
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

/**
 * @brief Real micro-optimization Matrix + Matrix - no extra creation of variable and copy.
 * @internal
 * @tparam Mtype1
 * @tparam Mtype2
 * @tparam Mtype2,
 * std::enable_if_t<Mtype1::is_matrix() && Mtype2::is_matrix(), int,
 * typename Rtype, Mtype2>
 * @param a
 * @param b
 * @return Rtype Return matrix of compatible type between Mtype1 and Mtype2
 */
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

/**
 * @brief Subtraction operator Matrix - Matrix
 * @memberof Matrix_t
 * @details
 * Subtraction operator between matrices is defined in the usual way (element wise).
 *
 * __NOTE__: Matrices must share same dimensionality.
 *
 * \code {.cpp}
 * Matrix<n,m,MyType> M, N, S;
 * M.fill(1);
 * N.fill(1);
 * S = M - N; // Resulting matrix is uniformly 0
 * \endcode
 *
 * @tparam Mtype1 Matrix type of a
 * @tparam Mtype2 Matrix type of b
 * @param a Matrix to subtract from
 * @param b Matrix to subtract
 * @return Rtype Return matrix of compatible type between Mtype1 and Mtype2
 */
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

/**
 * @brief Addition operator Matrix + scalar
 * @memberof Matrix_t
 * @details
 * Addition operator between matrix and scalar is defined in the usual way, where the scalar is
 * added to the diagonal elements.
 * *
 * \f$ M + 2 = M + 2\cdot\mathbb{1} \f$
 *
 * \code {.cpp}
 * Matrix<n,m,MyType> M,S;
 * M.fill(0);
 * S = M + 1; // Resulting matrix is identity matrix
 * \endcode
 *
 * @tparam Mtype Matrix type of b
 * @tparam S Type of scalar
 * @param a Matrix to add to
 * @param b Scalar to add
 * @return Rtype
 */
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

/**
 * @brief Addition operator scalar + Matrix
 * @memberof Matrix_t
 * @details
 * Addition operator between Scalar and Matrix is defined in the usual way, where the scalar is
 * treated as diagonal matrix which is then added to.
 * \f$ 2 + M = 2\cdot\mathbb{1} + M \f$
 *
 * \code {.cpp}
 * Matrix<n,m,MyType> M,R;
 * M = 0; // M = 0*I
 * R = 1 + M; // Resulting matrix is identity matrix
 * \endcode
 * @tparam Mtype Matrix type of a
 * @tparam S scalar type of b
 * @param b Matrix to add
 * @param a Scalar to add to
 * @return Rtype
 */
template <typename Mtype, typename S,
          std::enable_if_t<Mtype::is_matrix() && hila::is_complex_or_arithmetic<S>::value, int> = 0,
          typename Rtype = hila::mat_scalar_type<Mtype, S>>
inline Rtype operator+(const S &b, const Mtype &a) {

    static_assert(Mtype::rows() == Mtype::columns(),
                  "Matrix + scalar possible only for square matrix");
    return a + b;
}

// matrix - scalar
/**
 * @brief Subtraction operator Matrix - scalar
 * @memberof Matrix_t
 * @details
 *
 * Subtraction operator between matrix and scalar is defined in the usual way, where the scalar
 * is subtracted from the diagonal elements.
 *
 * \f$ M - 2 = M - 2\cdot\mathbb{1} \f$
 *
 * \code {.cpp}
 * Matrix<n,m,MyType> M,R;
 * M = 2; // M = 2*I
 * R = M - 1; // Resulting matrix is identity matrix
 * \endcode
 * @tparam Mtype Matrix type of a
 * @tparam S Scalar type of b
 * @param a Matrix to subtract from
 * @param b Scalar to subtract
 * @return Rtype
 */
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

/**
 * @brief Subtraction operator Scalar - Matrix
 * @memberof Matrix_t
 * @details
 * Subtraction operator between Scalar and Matrix is defined in the usual way, where the scalar is
 * treated as diagonal matrix which is then subtracted from.
 *
 * \f$ 2 - M = 2\cdot\mathbb{1} - M \f$
 *
 * \code {.cpp}
 * Matrix<n,m,MyType> M,R;
 * M = 1; // M = 1*I
 * R = 2 - M; // Resulting matrix is identity matrix
 * \endcode
 * @tparam Mtype Matrix type a
 * @tparam S Scalar type of b
 * @param b Scalar to subtract from
 * @param a Matrix to subtract
 * @return Rtype
 */
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
// matrix * matrix is the most important bit

/**
 * @brief Multiplication operator Square Matrix * Square Matrix
 * @internal
 * @memberof Matrix_t
 * @tparam Mt Matrix type
 * @param a First matrix to multiply
 * @param b Second matrix to multiply
 * @return Mt
 */
template <typename Mt, std::enable_if_t<Mt::is_matrix() && Mt::rows() == Mt::columns(), int> = 0>
inline Mt operator*(const Mt &a, const Mt &b) {

    constexpr int n = Mt::rows();

    Mt res;

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            res.e(i, j) = a.e(i, 0) * b.e(0, j);
            for (int k = 1; k < n; k++) {
                res.e(i, j) += a.e(i, k) * b.e(k, j);
            }
        }
    return res;
}

/**
 * @brief Multiplication operator
 * @memberof Matrix_t
 * @details Multiplication type depends on the original types of the multiplied matrices. Defined
 * for the following operations.
 *
 * Matrix * Matrix
 *
 * [Standard matrix multiplication](https://en.wikipedia.org/wiki/Matrix_multiplication) where LHS
 * columns must match RHS rows
 *
 * \code {.cpp}
 * Matrix<2, 3, double> M;
 * Matrix<3, 2, double> N;
 * M.random();
 * N.random();
 * auto S = N * M; // Results in a 3 x 3 Matrix since N.rows() = 3 and M.columns = 3
 * \endcode
 *
 * #Vector * RowVector is same as outer product which is equivalent to a matrix
 * multiplication
 *
 * \code {.cpp}
 * auto S = W*V // Result in n x n Matrix of type MyType
 * \endcode
 *
 * @tparam Mt1 Matrix type for a
 * @tparam Mt2 Matrix type for b
 * @tparam n Number of rows
 * @tparam m Number of columns
 * @param a Left Matrix or Vector
 * @param b Right Matrix or RowVector
 * @return Matrix<n, m, R>
 */
template <typename Mt1, typename Mt2,
          std::enable_if_t<Mt1::is_matrix() && Mt2::is_matrix() && !std::is_same<Mt1, Mt2>::value,
                           int> = 0,
          typename R = hila::ntype_op<hila::number_type<Mt1>, hila::number_type<Mt2>>,
          int n = Mt1::rows(), int m = Mt2::columns()>
inline Matrix<n, m, R> operator*(const Mt1 &a, const Mt2 &b) {

    constexpr int p = Mt1::columns();
    static_assert(p == Mt2::rows(), "Matrix size: LHS columns != RHS rows");

    Matrix<n, m, R> res;

    if constexpr (n > 1 && m > 1) {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) {
                res.e(i, j) = a.e(i, 0) * b.e(0, j);
                for (int k = 1; k < p; k++) {
                    res.e(i, j) += a.e(i, k) * b.e(k, j);
                }
            }
    } else if constexpr (m == 1) {
        // matrix x vector
        for (int i = 0; i < n; i++) {
            res.e(i) = a.e(i, 0) * b.e(0);
            for (int k = 1; k < p; k++) {
                res.e(i) += a.e(i, k) * b.e(k);
            }
        }
    } else if constexpr (n == 1) {
        // vector x matrix
        for (int j = 0; j < m; j++) {
            res.e(j) = a.e(0) * b.e(0, j);
            for (int k = 1; k < p; k++) {
                res.e(j) += a.e(k) * b.e(k, j);
            }
        }
    }

    return res;
}

/**
 * @brief Multiplication operator RowVector * Vector
 * @memberof Matrix_t
 * @details
 * Defined as standard [dot product](https://en.wikipedia.org/wiki/Dot_product) between
 * vectors as long as vectors are of same length
 *
 * \code {.cpp}
 * RowVector<n,MyType> V;
 * Vector<n,MyType> W;
 * //
 * // Fill V and W
 * //
 * auto S = V*W; // Results in MyType scalar
 * \endcode
 *
 * @tparam m Dimensions of RowVector
 * @tparam n Dimension of Vector
 * @tparam T1 Element type of RowVector
 * @tparam T2 Element type of Vector
 * @tparam Rt Return type of T1 \f$ \cdot \f$ T2 product
 * @param A RowVector to multiply
 * @param B Vector to multiply
 * @return Rt
 */
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

/**
 * @brief Multiplication operator Matrix * scalar
 * @memberof Matrix_t
 * @details Multiplication operator between Matrix and Scalar are defined in the usual way (element
 * wise)
 *
 * \code {.cpp}
 * Matrix<n,n,MyType> M;
 * M.fill(1);
 * auto S = M*2; // Resulting Matrix is uniformly 2
 * \endcode
 *
 * @tparam Mt Matrix type of mat
 * @tparam s Scalar type rhs
 * @param mat Matrix to multiply
 * @param rhs Scalar to multiply
 * @return Rt
 */
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

/**
 * @brief Multiplication operator Scalar * Matrix
 * @memberof Matrix_t
 * @details Multiplication operator between Matrix and Scalar are defined in the usual way (element
 * wise)
 *
 * \code {.cpp}
 * Matrix<n,n,MyType> M;
 * M.fill(1);
 * auto S = 2*M; // Resulting Matrix is uniformly 2
 * \endcode
 *
 * @tparam Mt Matrix type of mat
 * @tparam s Scalar type rhs
 * @param mat Matrix to multiply
 * @param rhs Scalar to multiply
 * @return Rt
 */
template <typename Mt, typename S,
          std::enable_if_t<(Mt::is_matrix() && hila::is_complex_or_arithmetic<S>::value), int> = 0,
          typename Rt = hila::mat_scalar_type<Mt, S>>
inline Rt operator*(const S &rhs, const Mt &mat) {
    return mat * rhs; // assume commutes
}

// matrix / scalar

/**
 * @brief Division operator
 * @memberof Matrix_t
 * @details Defined for following operations
 *
 * __Matrix / Scalar:__
 *
 * Division operator between Matrix and Scalar are defined in the usual way (element wise)
 *
 * \code {.cpp}
 * Matrix<n,m,MyType> M;
 * M.fill(2);
 * auto S = M/2; // Resulting matrix is uniformly 1
 * \endcode
 *
 * @tparam Mt Matrix type
 * @tparam S Scalar type
 * @param mat Matrix to divide scalar with
 * @param rhs Scalar to divide with
 * @return Rt Resulting Matrix
 */
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


/**
 * @brief Returns trace of product of two matrices
 * @details Wrapper around Matrix::mul_trace in the form
 * \code {.cpp}
 * mul_trace(a,b) // -> a.mul_trace(b)
 * \endcode
 *
 * @tparam Mtype1 Matrix type for a
 * @tparam Mtype2 Matrix type for b
 * @param a Left Matrix
 * @param b Right Matrix
 * @return auto Resulting trace
 */
template <typename Mtype1, typename Mtype2,
          std::enable_if_t<Mtype1::is_matrix() && Mtype2::is_matrix(), int> = 0>
inline auto mul_trace(const Mtype1 &a, const Mtype2 &b) {

    static_assert(Mtype1::columns() == Mtype2::rows() && Mtype1::rows() == Mtype2::columns(),
                  "mul_trace(a,b): matrix sizes are not compatible");
    return a.mul_trace(b);
}

//////////////////////////////////////////////////////////////////////////////////


// Stream operator
/**
 * @brief Stream operator
 * @memberof Matrix
 * @details Naive Stream operator for formatting Matrix for printing. Simply puts elements one after
 * the other in row major order
 *
 * @tparam n
 * @tparam m
 * @tparam T
 * @tparam MT
 * @param strm
 * @param A
 * @return std::ostream&
 */
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

/**
 * @brief Converts Matrix_t object to string
 *
 * @tparam n Number of rows
 * @tparam m Number of columns
 * @tparam T Matrix element type
 * @tparam MT Matrix type
 * @param A Matrix to convert to string
 * @param prec Precision of T
 * @param separator Separator between elements
 * @return std::string
 */
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

/**
 * @brief Formats Matrix_t object in a human readable way
 * @details Example 2 x 3 matrix is the following
 * \code {.cpp}
 *  Matrix<2, 3, double> W;
 *  W.random();
 *  hila::out0 << hila::prettyprint(W ,4) << '\n';
 *  // Result:
 *  // [ 0.8555 0.006359 0.3193 ]
 *  // [  0.237   0.8871 0.8545 ]
 * \endcode
 *
 *
 * @tparam n Number of rows
 * @tparam m Number of columns
 * @tparam T Matrix element type
 * @tparam MT Matrix type
 * @param A Matrix to be formatted
 * @param prec Precision of Matrix element
 * @return std::string
 */
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

/**
 * @brief Returns square norm of Matrix
 * @details Wrapper around Matrix::squarenorm - sum of squared elements
 *
 * Can be called as:
 * \code {.cpp}
 * Matrix<n,m,MyType> M;
 * auto a = squarenorm(M);
 * \endcode
 *
 * @tparam Mt Matrix type
 * @param rhs Matrix to compute square norm of
 * @return auto
 */
template <typename Mt, std::enable_if_t<Mt::is_matrix(), int> = 0>
inline auto squarenorm(const Mt &rhs) {
    return rhs.squarenorm();
}

// Vector norm - sqrt of squarenorm()

/**
 * @brief Returns vector norm of Matrix
 * @details Wrapper around Matrix::norm - sqrt of sum of squared elements
 * @tparam Mt Matrix type
 * @param rhs Matrix to compute norm of
 * @return auto
 */
template <typename Mt, std::enable_if_t<Mt::is_matrix(), int> = 0>
inline auto norm(const Mt &rhs) {
    return rhs.norm();
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

//  Calculate exp of a square matrix
//  Go to  order ORDER in the exponential of the matrices
//  matrices, since unitarity seems important.
//  Evaluation is done as:
// 	exp(H) = 1 + H + H^2/2! + H^3/3! ..-
//	           = 1 + H*( 1 + (H/2)*( 1 + (H/3)*( ... )))
//  Do it backwards in order to reduce accumulation of errors

/**
 * @brief Calculate exp of a square matrix
 * @memberof Matrix
 * @details Computes up to certain order given as argument
 *
 * __Evaluation is done as__:
 * \f{align}{ \exp(H) &= 1 + H + \frac{H^2}{2!} + \frac{H^2}{2!} + \frac{H^3}{3!} \\
 *                    &= 1 + H\cdot(1 + (\frac{H}{2})\cdot (1 + (\frac{H}{3})\cdot(...))) \f}
 * Done backwards in order to reduce accumulation of errors
 * @tparam n Number of rowsMa
 * @tparam T Matrix element type
 * @tparam MT Matrix type
 * @param mat Matrix to compute exponential for
 * @param order order to compute exponential to
 * @return Matrix_t<n, m, T, MT>
 */
template <int n, int m, typename T, typename MT>
inline Matrix_t<n, m, T, MT> exp(const Matrix_t<n, m, T, MT> &mat, const int order = 20) {
    static_assert(n == m, "exp() only for square matrices");

    Matrix_t<n, m, T, MT> r;
    hila::arithmetic_type<T> one = 1.0;

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

#include "datatypes/diagonal_matrix.h"

#include "datatypes/matrix_linalg.h"

// #include "datatypes/dagger.h"

#endif
