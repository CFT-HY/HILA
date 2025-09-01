#ifndef DIAGONAL_MATRIX_H_
#define DIAGONAL_MATRIX_H_

#include "matrix.h"

/// Define type DiagonalMatrix<n,T>


/**
 * @brief Class for diagonal matrix
 *
 * More optimal storage and algebra than normal square matrix
 *
 * @tparam n dimensionality
 * @tparam T type
 */
template <int n, typename T>
class DiagonalMatrix {

  public:
    static_assert(hila::is_complex_or_arithmetic<T>::value,
                  "DiagonalMatrix requires Complex or arithmetic type");

    // std incantation for field types
    using base_type = hila::arithmetic_type<T>;
    using argument_type = T;

    T c[n];

    /// Define default constructors to ensure std::is_trivial
    DiagonalMatrix() = default;
    ~DiagonalMatrix() = default;
    DiagonalMatrix(const DiagonalMatrix &v) = default;

    // constructor from scalar -- keep it explicit!  Not good for auto use
    template <typename S, std::enable_if_t<(hila::is_assignable<T &, S>::value), int> = 0>
    explicit inline DiagonalMatrix(const S rhs) {
        for (int i = 0; i < n; i++)
            c[i] = rhs;
    }


    // Construct matrix automatically from right-size initializer list
    // This does not seem to be dangerous, so keep non-explicit
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    inline DiagonalMatrix(std::initializer_list<S> rhs) {
        assert(rhs.size() == n && "Matrix/Vector initializer list size must match variable size");
        int i = 0;
        for (auto it = rhs.begin(); it != rhs.end(); it++, i++) {
            c[i] = *it;
        }
    }

    /**
     * @brief Returns matrix size - all give same result
     */
    static constexpr int rows() {
        return n;
    }
    static constexpr int columns() {
        return n;
    }
    static constexpr int size() {
        return n;
    }

    /**
     * @brief Element access - e(i) gives diagonal element i
     *
     */
    
    #pragma hila loop_function
    inline T e(const int i) const {
        return c[i];
    }

    // Same as above but with const_function, see const_function for details
    inline T &e(const int i) const_function {
        return c[i];
    }

    // For completeness, add e(i,j) - only for const
    T e(const int i, const int j) const {
        T ret(0);
        if (i == j)
            ret = c[i];
        return ret;
    }


    /**
     * @brief Return row from Diagonal matrix
     *
     */
    RowVector<n, T> row(int i) const {
        RowVector<n, T> res = 0;
        res.e(i) = c[i];
        return res;
    }

    /**
     * @brief Returns column vector i
     *
     */
    Vector<n, T> column(int i) const {
        Vector<n, T> v = 0;
        v.e(i) = c[i];
        return v;
    }


    /**
     * @brief Unary - operator
     *
     */
    inline DiagonalMatrix<n, T> operator-() const {
        DiagonalMatrix<n, T> res;
        for (int i = 0; i < n; i++) {
            res.e(i) = -c[i];
        }
        return res;
    }

    /**
     * @brief Unary + operator
     */
    inline const auto &operator+() const {
        return *this;
    }


    /**
     * @brief Boolean operator != to check if matrices are exactly different
     */
    template <typename S>
    bool operator!=(const S &rhs) const {
        return !(*this == rhs);
    }


    /**
     * Assignment operators: assign from another DiagonalMatrix, scalar or initializer list
     */

    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    inline DiagonalMatrix &operator=(const DiagonalMatrix<n, S> &rhs) out_only & {

        for (int i = 0; i < n; i++) {
            c[i] = rhs.e(i);
        }
        return *this;
    }

    // Assign from "scalar"

    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    inline DiagonalMatrix &operator=(const S &rhs) out_only & {

        for (int i = 0; i < n; i++) {
            c[i] = rhs;
        }
        return *this;
    }

    // Assign from initializer list

    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    inline DiagonalMatrix &operator=(std::initializer_list<S> rhs) out_only & {
        assert(rhs.size() == n && "Initializer list has a wrong size in assignment");
        int i = 0;
        for (auto it = rhs.begin(); it != rhs.end(); it++, i++) {
            c[i] = *it;
        }
        return *this;
    }


    // +=

    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    DiagonalMatrix &operator+=(const DiagonalMatrix<n, S> &rhs) & {
        for (int i = 0; i < n; i++) {
            c[i] += rhs.e(i);
        }
        return *this;
    }


    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    DiagonalMatrix &operator+=(const S &rhs) & {
        for (int i = 0; i < n; i++) {
            c[i] += rhs;
        }
        return *this;
    }

    // -=

    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    DiagonalMatrix &operator-=(const DiagonalMatrix<n, S> &rhs) & {
        for (int i = 0; i < n; i++) {
            c[i] -= rhs.e(i);
        }
        return *this;
    }


    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    DiagonalMatrix &operator-=(const S &rhs) & {
        for (int i = 0; i < n; i++) {
            c[i] -= rhs;
        }
        return *this;
    }

    // *=

    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_mul<T, S>>::value, int> = 0>
    DiagonalMatrix &operator*=(const DiagonalMatrix<n, S> &rhs) & {
        for (int i = 0; i < n; i++) {
            c[i] *= rhs.e(i);
        }
        return *this;
    }


    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_mul<T, S>>::value, int> = 0>
    DiagonalMatrix &operator*=(const S &rhs) & {
        for (int i = 0; i < n; i++) {
            c[i] *= rhs;
        }
        return *this;
    }

    // /=

    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_div<T, S>>::value, int> = 0>
    DiagonalMatrix &operator/=(const DiagonalMatrix<n, S> &rhs) & {
        for (int i = 0; i < n; i++) {
            c[i] /= rhs.e(i);
        }
        return *this;
    }


    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_div<T, S>>::value, int> = 0>
    DiagonalMatrix &operator/=(const S &rhs) & {
        for (int i = 0; i < n; i++) {
            c[i] /= rhs;
        }
        return *this;
    }


    /**
     * @brief fill with constant value - same as assignment
     */
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    DiagonalMatrix &fill(const S rhs) out_only {
        for (int i = 0; i < n; i++)
            c[i] = rhs;
        return *this;
    }

    /**
     * @brief transpose - leaves diagonal matrix as is
     */
    DiagonalMatrix &transpose() const {
        return *this;
    }

    /**
     * @brief dagger - conjugate elements
     */
    DiagonalMatrix dagger() const {
        DiagonalMatrix<n, T> ret;
        for (int i = 0; i < n; i++)
            ret.e(i) = ::conj(c[i]);
        return ret;
    }

    DiagonalMatrix adjoint() const {
        return this->dagger();
    }

    /**
     * @brief conj - conjugate elements
     */
    DiagonalMatrix conj() const {
        return this->dagger();
    }


    /**
     * @brief absolute value of all elements
     *
     * Returns a real DiagonalMatrix even if the original is complex
     */
    auto abs() const {
        DiagonalMatrix<n, hila::arithmetic_type<T>> res;
        for (int i = 0; i < n; i++) {
            res.c[i] = ::abs(c[i]);
        }
        return res;
    }

    /**
     * @brief real and imaginary parts of diagonal matrix
     *
     * Returns a real DiagonalMatrix even if the original is complex
     */
    auto real() const {
        DiagonalMatrix<n, hila::arithmetic_type<T>> res;
        for (int i = 0; i < n; i++) {
            res.c[i] = ::real(c[i]);
        }
        return res;
    }

    auto imag() const {
        DiagonalMatrix<n, hila::arithmetic_type<T>> res;
        for (int i = 0; i < n; i++) {
            res.c[i] = ::imag(c[i]);
        }
        return res;
    }

    /**
     * @brief Find max or min value - only for arithmetic types
     */

    template <typename S = T, std::enable_if_t<hila::is_arithmetic<S>::value, int> = 0>
    T max() const {
        T res = c[0];
        for (int i = 1; i < n; i++) {
            if (res < c[i])
                res = c[i];
        }
        return res;
    }

    template <typename S = T, std::enable_if_t<hila::is_arithmetic<S>::value, int> = 0>
    T min() const {
        T res = c[0];
        for (int i = 1; i < n; i++) {
            if (res > c[i])
                res = c[i];
        }
        return res;
    }

    T trace() const {
        T res = c[0];
        for (int i = 1; i < n; i++)
            res += c[i];
        return res;
    }

    T det() const {
        T res = c[0];
        for (int i = 1; i < n; i++)
            res *= c[i];
        return res;
    }

    auto squarenorm() const {
        hila::arithmetic_type<T> res(0);
        for (int i = 0; i < n; i++)
            res += ::squarenorm(c[i]);
        return res;
    }

    hila::arithmetic_type<T> norm() const {
        return sqrt(squarenorm());
    }

    /**
     * @brief Fills Matrix with random elements
     * @details Works only for non-integer valued elements
     */
    DiagonalMatrix &random() out_only {

        static_assert(hila::is_floating_point<hila::arithmetic_type<T>>::value,
                      "DiagonalMatrix random() requires non-integral type elements");
        for (int i = 0; i < n; i++) {
            hila::random(c[i]);
        }
        return *this;
    }

    DiagonalMatrix &gaussian_random(double width = 1.0) out_only {

        static_assert(hila::is_floating_point<hila::arithmetic_type<T>>::value,
                      "DiagonaMatrix gaussian_random() requires non-integral type elements");
        for (int i = 0; i < n; i++) {
            hila::gaussian_random(c[i], width);
        }
        return *this;
    }


    /**
     * @brief convert to string for printing
     */
    std::string str(int prec = 8, char separator = ' ') const {
        return this->asArray().str(prec, separator);
    }

    /**
     * @brief convert to generic matrix
     */
    Matrix<n, n, T> toMatrix() const {
        Matrix<n, n, T> res;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j)
                    res.e(i, j) = c[i];
                else
                    res.e(i, j) = 0;
            }
        }
        return res;
    }


    // temp cast to Array, for some arithmetic ops

    Array<n, 1, T> &asArray() const_function {
        return *(reinterpret_cast<Array<n, 1, T> *>(this));
    }

    const Array<n, 1, T> &asArray() const {
        return *(reinterpret_cast<const Array<n, 1, T> *>(this));
    }

    Vector<n, T> &asVector() const_function {
        return *(reinterpret_cast<Vector<n, T> *>(this));
    }

    const Vector<n, T> &asVector() const {
        return *(reinterpret_cast<const Vector<n, T> *>(this));
    }

    /// implement sort as casting to matrix
#pragma hila novector
    template <int N>
    DiagonalMatrix<n, T> sort(Vector<N, int> &permutation,
                              hila::sort order = hila::sort::ascending) const {
        return this->asArray().sort(permutation, order).asDiagonalMatrix();
    }

#pragma hila novector
    DiagonalMatrix<n, T> sort(hila::sort order = hila::sort::ascending) const {
        return this->asArray().sort(order).asDiagonalMatrix();
    }
};

///////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Boolean operator == to determine if two diagonal matrices are exactly the same
 */
template <typename A, typename B, int n, int m>
inline bool operator==(const DiagonalMatrix<n, A> &lhs, const DiagonalMatrix<m, B> &rhs) {
    if constexpr (m != n)
        return false;

    for (int i = 0; i < n; i++) {
        if (lhs.e(i) != rhs.e(i))
            return false;
    }
    return true;
}

/**
 * @brief Boolean operator == to compare with square matrix
 * Matrices are equal if diagonals are equal and off-diag is zero
 */
template <typename A, typename S, typename Mtype, int n, int m1, int m2>
inline bool operator==(const DiagonalMatrix<n, A> &lhs, const Matrix_t<m1, m2, S, Mtype> &rhs) {
    if constexpr (m1 != n || m2 != n)
        return false;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                if (lhs.e(i) != rhs.e(i, i))
                    return false;
            } else {
                if (rhs.e(i, j) != 0)
                    return false;
            }
        }
    }
    return true;
}

template <typename A, typename S, typename Mtype, int n, int m1, int m2>
inline bool operator==(const Matrix_t<m1, m2, S, Mtype> &lhs, const DiagonalMatrix<n, A> &rhs) {
    return rhs == lhs;
}

// Compiler generates operator!= automatically

template <int n, typename T>
inline const auto &transpose(const DiagonalMatrix<n, T> &arg) {
    return arg;
}

template <int n, typename T>
inline auto dagger(const DiagonalMatrix<n, T> &arg) {
    return arg.dagger();
}

template <int n, typename T>
inline auto conj(const DiagonalMatrix<n, T> &arg) {
    return arg.conj();
}

template <int n, typename T>
inline auto adjoint(const DiagonalMatrix<n, T> &arg) {
    return arg.adjoint();
}

template <int n, typename T>
inline auto abs(const DiagonalMatrix<n, T> &arg) {
    return arg.abs();
}

template <int n, typename T>
inline auto real(const DiagonalMatrix<n, T> &arg) {
    return arg.real();
}

template <int n, typename T>
inline auto imag(const DiagonalMatrix<n, T> &arg) {
    return arg.imag();
}

template <int n, typename T>
inline auto trace(const DiagonalMatrix<n, T> &arg) {
    return arg.trace();
}

template <int n, typename T>
inline auto squarenorm(const DiagonalMatrix<n, T> &arg) {
    return arg.squarenorm();
}

template <int n, typename T>
inline auto norm(const DiagonalMatrix<n, T> &arg) {
    return arg.norm();
}

template <int n, typename T>
inline auto det(const DiagonalMatrix<n, T> &arg) {
    return arg.det();
}


namespace hila {

////////////////////////////////////////////////////////////////////////////////
// DiagonalMatrix + scalar result type:
// hila::diagonalmatrix_scalar_type<Mt,S>
//  - if result is convertible to Mt, return Mt
//  - if Mt is not complex and S is, return
//  DiagonalMatrix<Complex<type_sum(scalar_type(Mt),scalar_type(S))>>
//  - otherwise return DiagonalMatrix<type_sum>

template <typename Mt, typename S, typename Enable = void>
struct diagonalmatrix_scalar_op_s {
    using type = DiagonalMatrix<
        Mt::rows(), Complex<hila::type_plus<hila::arithmetic_type<Mt>, hila::arithmetic_type<S>>>>;
};

template <typename Mt, typename S>
struct diagonalmatrix_scalar_op_s<
    Mt, S,
    typename std::enable_if_t<std::is_convertible<hila::type_plus<hila::number_type<Mt>, S>,
                                                  hila::number_type<Mt>>::value>> {
    // using type = Mt;
    using type = typename std::conditional<
        hila::is_floating_point<hila::arithmetic_type<Mt>>::value, Mt,
        DiagonalMatrix<Mt::rows(),
                       hila::type_plus<hila::arithmetic_type<Mt>, hila::arithmetic_type<S>>>>::type;
};

template <typename Mt, typename S>
using diagonalmatrix_scalar_type = typename diagonalmatrix_scalar_op_s<Mt, S>::type;

} // namespace hila


/**
 * operators: can add scalar, diagonal matrix or square matrix.
 */

// diagonal + scalar
template <int n, typename T, typename S,
          std::enable_if_t<hila::is_complex_or_arithmetic<S>::value, int> = 0,
          typename Rtype = hila::diagonalmatrix_scalar_type<DiagonalMatrix<n, T>, S>>
inline Rtype operator+(const DiagonalMatrix<n, T> &a, const S &b) {
    Rtype res;
    for (int i = 0; i < n; i++)
        res.e(i) = a.e(i) + b;
    return res;
}

// diagonal + scalar
template <int n, typename T, typename S,
          std::enable_if_t<hila::is_complex_or_arithmetic<S>::value, int> = 0,
          typename Rtype = hila::diagonalmatrix_scalar_type<DiagonalMatrix<n, T>, S>>
inline Rtype operator+(const S &b, const DiagonalMatrix<n, T> &a) {
    return a + b;
}

// diagonal - scalar
template <int n, typename T, typename S,
          std::enable_if_t<hila::is_complex_or_arithmetic<S>::value, int> = 0,
          typename Rtype = hila::diagonalmatrix_scalar_type<DiagonalMatrix<n, T>, S>>
inline Rtype operator-(const DiagonalMatrix<n, T> &a, const S &b) {
    Rtype res;
    for (int i = 0; i < n; i++)
        res.e(i) = a.e(i) - b;
    return res;
}

// scalar - diagonal
template <int n, typename T, typename S,
          std::enable_if_t<hila::is_complex_or_arithmetic<S>::value, int> = 0,
          typename Rtype = hila::diagonalmatrix_scalar_type<DiagonalMatrix<n, T>, S>>
inline Rtype operator-(const S &b, const DiagonalMatrix<n, T> &a) {
    Rtype res;
    for (int i = 0; i < n; i++)
        res.e(i) = b - a.e(i);
    return res;
}

// diagonal * scalar
template <int n, typename T, typename S,
          std::enable_if_t<hila::is_complex_or_arithmetic<S>::value, int> = 0,
          typename Rtype = hila::diagonalmatrix_scalar_type<DiagonalMatrix<n, T>, S>>
inline Rtype operator*(const DiagonalMatrix<n, T> &a, const S &b) {
    Rtype res;
    for (int i = 0; i < n; i++)
        res.e(i) = a.e(i) * b;
    return res;
}

// scalar * diagonal
template <int n, typename T, typename S,
          std::enable_if_t<hila::is_complex_or_arithmetic<S>::value, int> = 0,
          typename Rtype = hila::diagonalmatrix_scalar_type<DiagonalMatrix<n, T>, S>>
inline Rtype operator*(const S &b, const DiagonalMatrix<n, T> &a) {
    return a * b;
}

// diagonal / scalar
template <int n, typename T, typename S,
          std::enable_if_t<hila::is_complex_or_arithmetic<S>::value, int> = 0,
          typename Rtype = hila::diagonalmatrix_scalar_type<DiagonalMatrix<n, T>, S>>
inline Rtype operator/(const DiagonalMatrix<n, T> &a, const S &b) {
    Rtype res;
    for (int i = 0; i < n; i++)
        res.e(i) = a.e(i) / b;
    return res;
}

// scalar / diagonal
template <int n, typename T, typename S,
          std::enable_if_t<hila::is_complex_or_arithmetic<S>::value, int> = 0,
          typename Rtype = hila::diagonalmatrix_scalar_type<DiagonalMatrix<n, T>, S>>
inline Rtype operator/(const S &b, const DiagonalMatrix<n, T> &a) {
    Rtype res;
    for (int i = 0; i < n; i++)
        res.e(i) = b / a.e(i);
    return res;
}

/////
// diagonal X diagonal
template <int n, typename A, typename B, typename R = hila::type_plus<A, B>>
inline auto operator+(const DiagonalMatrix<n, A> &a, const DiagonalMatrix<n, B> &b) {
    DiagonalMatrix<n, R> res;
    for (int i = 0; i < n; i++)
        res.e(i) = a.e(i) + b.e(i);
    return res;
}

template <int n, typename A, typename B, typename R = hila::type_minus<A, B>>
inline auto operator-(const DiagonalMatrix<n, A> &a, const DiagonalMatrix<n, B> &b) {
    DiagonalMatrix<n, R> res;
    for (int i = 0; i < n; i++)
        res.e(i) = a.e(i) - b.e(i);
    return res;
}

template <int n, typename A, typename B, typename R = hila::type_mul<A, B>>
inline auto operator*(const DiagonalMatrix<n, A> &a, const DiagonalMatrix<n, B> &b) {
    DiagonalMatrix<n, R> res;
    for (int i = 0; i < n; i++)
        res.e(i) = a.e(i) * b.e(i);
    return res;
}

template <int n, typename A, typename B, typename R = hila::type_div<A, B>>
inline auto operator/(const DiagonalMatrix<n, A> &a, const DiagonalMatrix<n, B> &b) {
    DiagonalMatrix<n, R> res;
    for (int i = 0; i < n; i++)
        res.e(i) = a.e(i) / b.e(i);
    return res;
}

//// Finally, diagonal X Matrix ops - gives Matrix

template <int n, typename T, typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0,
          typename Rtype = hila::mat_x_mat_type<Matrix<n, n, T>, Mtype>>
inline Rtype operator+(const DiagonalMatrix<n, T> &a, const Mtype &b) {

    constexpr int mr = Mtype::rows();
    constexpr int mc = Mtype::columns();

    static_assert(mc == n && mr == n, "Matrix sizes do not match");

    Rtype r;
    r = b;
    for (int i = 0; i < n; i++)
        r.e(i, i) += a.e(i);
    return r;
}

template <int n, typename T, typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0,
          typename Rtype = hila::mat_x_mat_type<Matrix<n, n, T>, Mtype>>
inline Rtype operator+(const Mtype &b, const DiagonalMatrix<n, T> &a) {
    return a + b;
}

template <int n, typename T, typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0,
          typename Rtype = hila::mat_x_mat_type<Matrix<n, n, T>, Mtype>>
inline Rtype operator-(const DiagonalMatrix<n, T> &a, const Mtype &b) {

    constexpr int mr = Mtype::rows();
    constexpr int mc = Mtype::columns();

    static_assert(mc == n && mr == n, "Matrix sizes do not match");

    Rtype r = -b;
    for (int i = 0; i < n; i++)
        r.e(i, i) += a.e(i);
    return r;
}

template <int n, typename T, typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0,
          typename Rtype = hila::mat_x_mat_type<Matrix<n, n, T>, Mtype>>
inline Rtype operator-(const Mtype &b, const DiagonalMatrix<n, T> &a) {
    constexpr int mr = Mtype::rows();
    constexpr int mc = Mtype::columns();

    static_assert(mc == n && mr == n, "Matrix sizes do not match");

    Rtype r = b;
    for (int i = 0; i < n; i++)
        r.e(i, i) -= a.e(i);
    return r;
}


// multiply by matrix

template <int n, typename T, typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0,
          typename Rtype = hila::mat_x_mat_type<Matrix<n, n, T>, Mtype>>
inline Rtype operator*(const DiagonalMatrix<n, T> &a, const Mtype &b) {

    constexpr int mr = Mtype::rows();
    constexpr int mc = Mtype::columns();

    static_assert(mr == n, "Matrix sizes do not match");

    Rtype r;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < mc; j++)
            r.e(i, j) = a.e(i) * b.e(i, j);
    return r;
}

template <int n, typename T, typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0,
          typename Rtype = hila::mat_x_mat_type<Matrix<n, n, T>, Mtype>>
inline Rtype operator*(const Mtype &b, const DiagonalMatrix<n, T> &a) {

    constexpr int mr = Mtype::rows();
    constexpr int mc = Mtype::columns();

    static_assert(mc == n, "Matrix sizes do not match");

    Rtype r;
    for (int i = 0; i < mr; i++)
        for (int j = 0; j < n; j++)
            r.e(i, j) = b.e(i, j) * a.e(j);
    return r;
}

// division
template <int n, typename T, typename Mtype, std::enable_if_t<Mtype::is_matrix(), int> = 0,
          typename Rtype = hila::mat_x_mat_type<Matrix<n, n, T>, Mtype>>
inline Rtype operator/(const Mtype &b, const DiagonalMatrix<n, T> &a) {

    constexpr int mr = Mtype::rows();
    constexpr int mc = Mtype::columns();

    static_assert(mc == n, "Matrix sizes do not match");

    Rtype r;
    for (int i = 0; i < mr; i++)
        for (int j = 0; j < n; j++)
            r.e(i, j) = b.e(i, j) / a.e(j);
    return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Standard arithmetic functions - do element by element
////////////////////////////////////////////////////////////////////////////////

template <int n, typename T>
inline DiagonalMatrix<n, T> sqrt(DiagonalMatrix<n, T> a) {
    for (int i = 0; i < n; i++)
        a.c[i] = sqrt(a.c[i]);
    return a;
}

template <int n, typename T>
inline DiagonalMatrix<n, T> cbrt(DiagonalMatrix<n, T> a) {
    for (int i = 0; i < n; i++)
        a.c[i] = cbrt(a.c[i]);
    return a;
}

template <int n, typename T>
inline DiagonalMatrix<n, T> exp(DiagonalMatrix<n, T> a) {
    for (int i = 0; i < n; i++)
        a.c[i] = exp(a.c[i]);
    return a;
}

template <int n, typename T>
inline DiagonalMatrix<n, T> log(DiagonalMatrix<n, T> a) {
    for (int i = 0; i < n; i++)
        a.c[i] = log(a.c[i]);
    return a;
}

template <int n, typename T>
inline DiagonalMatrix<n, T> sin(DiagonalMatrix<n, T> a) {
    for (int i = 0; i < n; i++)
        a.c[i] = sin(a.c[i]);
    return a;
}

template <int n, typename T>
inline DiagonalMatrix<n, T> cos(DiagonalMatrix<n, T> a) {
    for (int i = 0; i < n; i++)
        a.c[i] = cos(a.c[i]);
    return a;
}

template <int n, typename T>
inline DiagonalMatrix<n, T> tan(DiagonalMatrix<n, T> a) {
    for (int i = 0; i < n; i++)
        a.c[i] = tan(a.c[i]);
    return a;
}

template <int n, typename T>
inline DiagonalMatrix<n, T> asin(DiagonalMatrix<n, T> a) {
    for (int i = 0; i < n; i++)
        a.c[i] = asin(a.c[i]);
    return a;
}

template <int n, typename T>
inline DiagonalMatrix<n, T> acos(DiagonalMatrix<n, T> a) {
    for (int i = 0; i < n; i++)
        a.c[i] = acos(a.c[i]);
    return a;
}

template <int n, typename T>
inline DiagonalMatrix<n, T> atan(DiagonalMatrix<n, T> a) {
    for (int i = 0; i < n; i++)
        a.c[i] = atan(a.c[i]);
    return a;
}

template <int n, typename T>
inline DiagonalMatrix<n, T> sinh(DiagonalMatrix<n, T> a) {
    for (int i = 0; i < n; i++)
        a.c[i] = sinh(a.c[i]);
    return a;
}

template <int n, typename T>
inline DiagonalMatrix<n, T> cosh(DiagonalMatrix<n, T> a) {
    for (int i = 0; i < n; i++)
        a.c[i] = cosh(a.c[i]);
    return a;
}

template <int n, typename T>
inline DiagonalMatrix<n, T> tanh(DiagonalMatrix<n, T> a) {
    for (int i = 0; i < n; i++)
        a.c[i] = tanh(a.c[i]);
    return a;
}

template <int n, typename T>
inline DiagonalMatrix<n, T> asinh(DiagonalMatrix<n, T> a) {
    for (int i = 0; i < n; i++)
        a.c[i] = asinh(a.c[i]);
    return a;
}

template <int n, typename T>
inline DiagonalMatrix<n, T> acosh(DiagonalMatrix<n, T> a) {
    for (int i = 0; i < n; i++)
        a.c[i] = acosh(a.c[i]);
    return a;
}

template <int n, typename T>
inline DiagonalMatrix<n, T> atanh(DiagonalMatrix<n, T> a) {
    for (int i = 0; i < n; i++)
        a.c[i] = atanh(a.c[i]);
    return a;
}


// return pow of diagonalMatrix as original type if power is scalar or diagonal is complex
template <
    int n, typename T, typename S,
    std::enable_if_t<hila::is_arithmetic<S>::value || hila::contains_complex<T>::value, int> = 0>
inline DiagonalMatrix<n, T> pow(DiagonalMatrix<n, T> a, S p) {
    for (int i = 0; i < n; i++)
        a.c[i] = pow(a.c[i], p);
    return a;
}

// if power is complex but matrix is scalar need to upgrade return type
template <int n, typename T, typename S,
          std::enable_if_t<!hila::contains_complex<T>::value, int> = 0>
inline auto pow(const DiagonalMatrix<n, T> &a, const Complex<S> &p) {
    DiagonalMatrix<n, hila::type_mul<T, Complex<S>>> res;
    for (int i = 0; i < n; i++)
        res.e(i) = pow(a.e(i), p);
    return res;
}

/// Stream operator
template <int n, typename T>
std::ostream &operator<<(std::ostream &strm, const DiagonalMatrix<n, T> &A) {
    return operator<<(strm, A.asArray());
}

namespace hila {

/**
 * @brief Cast to different basic number:
 * 
 * @details
 * hila::cast_to<Ntype>(a);
 * 
 * where a is DiagonalMatrix<n,T> does 
 * DiagonalMatrix<n,T> -> DiagonalMatrix<n,Ntype>
 * 
 */

template <typename Ntype, typename T, int n,
          std::enable_if_t<hila::is_arithmetic_or_extended<T>::value, int> = 0>
DiagonalMatrix<n, Ntype> cast_to(const DiagonalMatrix<n, T> &mat) {
    DiagonalMatrix<n, Ntype> res;
    for (int i = 0; i < n; i++)
        res.c[i] = cast_to<Ntype>(mat.c[i]);
    return res;
}

template <typename Ntype, typename T, int n, std::enable_if_t<hila::is_complex<T>::value, int> = 0>
DiagonalMatrix<n, Complex<Ntype>> cast_to(const DiagonalMatrix<n, T> &mat) {
    DiagonalMatrix<n, Complex<Ntype>> res;
    for (int i = 0; i < n; i++)
        res.c[i] = cast_to<Ntype>(mat.c[i]);
    return res;
}

template <int n, typename T>
std::string to_string(const DiagonalMatrix<n, T> &A, int prec = 8, char separator = ' ') {
    return to_string(A.asArray(), prec, separator);
}

template <int n, typename T>
std::string prettyprint(const DiagonalMatrix<n, T> &A, int prec = 8) {
    return prettyprint(A.toMatrix(), prec);
}

} // namespace hila


#endif