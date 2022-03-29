#ifndef DAGGER_H_
#define DAGGER_H_

#include <type_traits>
#include <sstream>
#include "plumbing/defs.h"
#include "datatypes/cmplx.h"
#include "datatypes/matrix.h"

////////////////////////////////////////////////////////////////
/// The dagger of the nxm matrix - special type for arithmetics (is mxn matrix)
///   Matrix<n,m,T>.dagger()  ->   DaggerMatrix<m,n,T>
////////////////////////////////////////////////////////////////
template <const int n, const int m, typename T>
class DaggerMatrix {
  public:
    static_assert(hila::is_complex_or_arithmetic<T>::value,
                  "DaggerMatrix requires Complex or arithmetic type");

    /// std incantation for field types
    using base_type = hila::number_type<T>;
    using argument_type = T;

    /// Same as Matrix, data as a one dimensional array
    T c[n * m];

    /// define default constructors to ensure std::is_trivial
    DaggerMatrix() =
        delete; // this should mean variables of the type cannot be declared
    ~DaggerMatrix() = default;
    DaggerMatrix(const DaggerMatrix<n, m, T> &v) = delete;

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
    /// For Matrix<n,m>:  e(i,j) = c[i*m + j]
    /// So dagger is c[j*n + i] conjugate
    inline T e(const int i, const int j) const {
        return ::conj(c[j * n + i]);
    }
    /// no ref accessor, cannot assign to a dagger matrix.

    /// declare single e here too in case we have a vector
    /// (one size == 1)
    template <int q = n, int p = m, std::enable_if_t<(q == 1 || p == 1), int> = 0>
    inline T e(const int i) const {
        return ::conj(c[i]);
    }

    /// And also [] for vectors (not matrices!)
    /// (one size == 1)
    template <int q = n, int p = m, std::enable_if_t<(q == 1 || p == 1), int> = 0>
    inline T operator[](const int i) const {
        return ::conj(c[i]);
    }

    Matrix<1, m, T> row(int r) const {
        Matrix<1, m, T> v;
        for (int i = 0; i < m; i++) v[i] = e(r, i);
        return v;
    }

    Matrix<n, 1, T> column(int c) const {
        Matrix<n, 1, T> v;
        for (int i = 0; i < n; i++) v[i] = e(i, c);
        return v;
    }

    /// interpret as Array -  for array ops
    const Array<n, m, T> asArray() const {
        Array<n, m, T> v;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) v.e(i, j) = e(i, j);
        return v;
    }

    /// automatic cast to matrix - hopefully takes care of many ops
    operator Matrix<n, m, T>() const {
        Matrix<n, m, T> v;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) v.e(i, j) = e(i, j);
        return v;
    }

    /// casting from one Matrix (number) type to another: do not do this automatically.
    /// but require an explicit cast operator.  This makes it easier to write code.
    /// or should it be automatic?  keep/remove explicit?
    /// TODO: CHECK AVX CONVERSIONS

    template <typename S, std::enable_if_t<hila::is_assignable<S &, T>::value, int> = 0>
    explicit operator Matrix<n, m, S>() const {
        Matrix<n, m, S> v;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) v.e(i, j) = e(i, j);
        return v;
    }

    /// unary -
    inline Matrix<n, m, T> operator-() const {
        Matrix<n, m, T> v;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) v.e(i, j) = -e(i, j);
        return v;
    }

    /// unary +
    inline DaggerMatrix<n, m, T> operator+() const {
        return *this;
    }

    // no assignments!

    /// numpy style fill
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    DaggerMatrix<n, m, T> &fill(const S rhs) out_only {
        for (int i = 0; i < n * m; i++) c[i] = rhs;
        return *this;
    }

    // return copy of transpose of this matrix
    //
    inline Matrix<m, n, T> transpose() const {
        return dagger().conj();
    }

    /// and complex conjugate
    inline Matrix<n, m, T> conj() const {
        return dagger().transpose();
    }

    /// return copy of adjoint/dagger of this Matrix
    inline const Matrix<m, n, T> &dagger() const {
        return *reinterpret_cast<const Matrix<m, n, T> *>(this);
    }

    /// alias dagger() to adjoint()
    inline const Matrix<m, n, T> &adjoint() const {
        return dagger();
    }

    /// return real part
    inline Matrix<n, m, hila::number_type<T>> real() const {
        Matrix<n, m, hila::number_type<T>> v;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) v.e(i, j) = real(e(i, j));
        return v;
    }

    /// return imaginary part
    inline Matrix<n, m, hila::number_type<T>> imag() const {
        Matrix<n, m, hila::number_type<T>> v;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) v.e(i, j) = imag(e(i, j));
        return v;
    }

    /// matrix trace (return type T)
    T trace() const {
        static_assert(n == m, "trace not defined for non square matrices!");
        T result = 0;
        for (int i = 0; i < n; i++) {
            result += e(i, i);
        }
        return result;
    }

    /// calculate square norm - sum of squared elements
    hila::number_type<T> squarenorm() const {
        hila::number_type<T> result = 0;
        for (int i = 0; i < n * m; i++) {
            result += ::squarenorm(c[i]);
        }
        return result;
    }

    /// calculate vector norm - sqrt of squarenorm
    hila::number_type<T> norm() const {
        return sqrt(squarenorm());
    }

    /// dot product - vector . vector
    template <int p, std::enable_if_t<(p == 1 && n == 1), int> = 0>
    inline T dot(const Matrix<m, p, T> &rhs) const {
        T r = 0;
        for (int i = 0; i < m; i++) {
            r += ::conj(c[i]) * rhs.c[i];
        }
        return r;
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

    /// Convert to string for printing
    std::string str() const {
        return ((Matrix<n, m, T>)*this).str();
    }
};

/// do transpose and adjoint functions here
template <const int n, const int m, typename T>
inline Matrix<n, m, T> transpose(const DaggerMatrix<m, n, T> &arg) {
    return arg.transpose();
}
/// conjugate
template <const int n, const int m, typename T>
inline Matrix<n, m, T> conj(const DaggerMatrix<n, m, T> &arg) {
    return arg.conj();
}
/// adjoint (=conjugate)
template <const int n, const int m, typename T>
inline Matrix<n, m, T> adjoint(const DaggerMatrix<m, n, T> &arg) {
    return arg.adjoint();
}
/// dagger (=conjugate)
template <const int n, const int m, typename T>
inline Matrix<n, m, T> dagger(const DaggerMatrix<m, n, T> &arg) {
    return arg.adjoint();
}
/// trace
template <const int n, const int m, typename T>
inline T trace(const DaggerMatrix<n, m, T> &arg) {
    return arg.trace();
}
/// real part
template <const int n, const int m, typename T>
inline Matrix<n, m, hila::number_type<T>> real(const DaggerMatrix<n, m, T> &arg) {
    return arg.real();
}
/// imaginary part
template <const int n, const int m, typename T>
inline Matrix<n, m, hila::number_type<T>> imag(const DaggerMatrix<n, m, T> &arg) {
    return arg.imag();
}

/// determinant -> use LU factorization later
template <int n, int m, typename T>
T det(const DaggerMatrix<n, m, T> &mat) {
    static_assert(n == m, "determinants defined only for square matrices");
    return det((Matrix<n, m, T>)mat);
}

// /// Now matrix additions: matrix + matrix
// template <int n, int m, typename T>
// inline Matrix<n, m, T> operator+(Matrix<n, m, T> a, const Matrix<n, m, T> &b) {
//     a += b;
//     return a;
// }

// /// Matrix subtract
// template <int n, int m, typename T>
// inline Matrix<n, m, T> operator-(Matrix<n, m, T> a, const Matrix<n, m, T> &b) {
//     a -= b;
//     return a;
// }

// /// Matrix + scalar
// template <
//     int n, int m, typename T, typename S,
//     std::enable_if_t<std::is_convertible<hila::type_plus<T, S>, T>::value, int> = 0>
// inline Matrix<n, m, T> operator+(Matrix<n, m, T> a, const S b) {
//     a += b;
//     return a;
// }

// /// scalar + matrix
// template <
//     int n, int m, typename T, typename S,
//     std::enable_if_t<std::is_convertible<hila::type_plus<T, S>, T>::value, int> = 0>
// inline Matrix<n, m, T> operator+(const S b, Matrix<n, m, T> a) {
//     a += b;
//     return a;
// }

// /// matrix - scalar
// template <
//     int n, int m, typename T, typename S,
//     std::enable_if_t<std::is_convertible<hila::type_minus<T, S>, T>::value, int> = 0>
// Matrix<n, m, T> operator-(Matrix<n, m, T> a, const S b) {
//     a -= b;
//     return a;
// }

// /// scalar - matrix
// template <
//     int n, int m, typename T, typename S,
//     std::enable_if_t<std::is_convertible<hila::type_minus<S, T>, T>::value, int> = 0>
// inline Matrix<n, m, T> operator-(const S b, Matrix<n, m, T> a) {
//     static_assert(
//         n == m,
//         "rows != columns : scalar subtraction possible for square matrices only!");
//     for (int i = 0; i < n; i++) a.e(i, i) = b - a.e(i, i);
//     return a;
// }

////////////////////////////////////////
/// matrix * matrix is the most important bit
template <
    int n, int m, int p, typename T, template <int, int, typename> typename M1,
    template <int, int, typename> typename M2,
    std::enable_if_t<(std::is_same<Matrix<n, m, T>, M1<n, m, T>>::value ||
                      std::is_same<DaggerMatrix<n, m, T>, M1<n, m, T>>::value) &&
                         (std::is_same<Matrix<m, p, T>, M2<m, p, T>>::value ||
                          std::is_same<DaggerMatrix<m, p, T>, M2<m, p, T>>::value),
                     int> = 0>
inline Matrix<n, p, T> operator*(const M1<n, m, T> &A, const M2<m, p, T> &B) {
    Matrix<n, p, T> res;

    if constexpr (n > 1 && p > 1) {
        // normal matrix*matrix
        for (int i = 0; i < n; i++)
            for (int j = 0; j < p; j++) {
                res.e(i, j) = 0;
                for (int k = 0; k < m; k++) {
                    res.e(i, j) += A.e(i, k) * B.e(k, j);
                }
            }
    } else if constexpr (p == 1) {
        // matrix * vector
        for (int i = 0; i < n; i++) {
            res.e(i) = 0;
            for (int k = 0; k < m; k++) {
                res.e(i) += A.e(i, k) * B.e(k);
            }
        }
    } else if constexpr (n == 1) {
        // horiz. vector * matrix
        for (int j = 0; j < p; j++) {
            res.e(j) = 0;
            for (int k = 0; k < m; k++) {
                res.e(j) += A.e(k) * B.e(k, j);
            }
        }
    }
    return res;
}

/// and treat separately horiz. vector * vector
template <
    int m, typename T, template <int, int, typename> typename M1,
    template <int, int, typename> typename M2,
    std::enable_if_t<(std::is_same<Matrix<1, m, T>, M1<1, m, T>>::value ||
                      std::is_same<DaggerMatrix<1, m, T>, M1<1, m, T>>::value) &&
                         (std::is_same<Matrix<m, 1, T>, M2<m, 1, T>>::value ||
                          std::is_same<DaggerMatrix<m, 1, T>, M2<m, 1, T>>::value),
                     int> = 0>
T operator*(const Matrix<1, m, T> &A, const Matrix<m, 1, T> &B) {
    T res = 0;
    for (int i = 0; i < m; i++) {
        res += A.e(i) * B.e(i);
    }
    return res;
}

/// matrix * scalar
template <
    int n, int m, typename T, typename S,
    std::enable_if_t<std::is_convertible<hila::type_mul<T, S>, T>::value, int> = 0>
Matrix<n, m, T> operator*(DaggerMatrix<n, m, T> mat, const S rhs) {
    Matrix<n, m, T> v;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++) v.e(i, j) = mat.e(i, j) * rhs;
    return v;
}

/// scalar * matrix
template <
    int n, int m, typename T, typename S,
    std::enable_if_t<std::is_convertible<hila::type_mul<S, T>, T>::value, int> = 0>
Matrix<n, m, T> operator*(const S rhs, DaggerMatrix<n, m, T> mat) {
    return mat * rhs; // assume commutativity, should be ok
}

/// matrix / scalar
template <
    int n, int m, typename T, typename S,
    std::enable_if_t<std::is_convertible<hila::type_div<T, S>, T>::value, int> = 0>
Matrix<n, m, T> operator/(DaggerMatrix<n, m, T> mat, const S rhs) {
    Matrix<n, m, T> v;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++) v.e(i, j) = mat.e(i, j) / rhs;
    return v;
}

/// Stream operator
template <int n, int m, typename T>
std::ostream &operator<<(std::ostream &strm, const DaggerMatrix<n, m, T> &A) {
    return operator<<(strm, (Matrix<n, m, T>)A);
}

/// Norm squared function
template <int n, int m, typename T>
inline hila::number_type<T> squarenorm(const DaggerMatrix<n, m, T> &rhs) {
    return rhs.squarenorm();
}

/// Vector norm - sqrt of squarenorm()
template <int n, int m, typename T>
inline hila::number_type<T> norm(const DaggerMatrix<n, m, T> &rhs) {
    return rhs.norm();
}

/// find determinant using LU decomposition. Algorithm: numerical Recipes, 2nd ed.
/// p. 47 ff
template <int n, typename T, typename radix = hila::number_type<T>,
          std::enable_if_t<hila::is_complex<T>::value, int> = 0>
Complex<radix> det_lu(const DaggerMatrix<n, n, T> &mat) {
    return det_lu((Matrix<n, n, T>)mat);
}

/// find determinant using LU decomposition. Algorithm: numerical Recipes, 2nd ed.
/// p. 47 ff
template <int n, typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
T det_lu(const DaggerMatrix<n, n, T> &mat) {
    return det_lu((Matrix<n, n, T>)mat);
}

// Cast operators to different number or Complex type
// cast_to<double>(a);
// cast_to<Complex<float>>(b);
// Cast from number->number, number->Complex, Complex->Complex OK,
//     Complex->number not.

template <typename Ntype, typename T, int n, int m,
          std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
Matrix<n, m, Ntype> cast_to(const DaggerMatrix<n, m, T> &mat) {
    Matrix<n, m, Ntype> res;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++) res.e(i, j) = mat.e(i, j);
    return res;
}

template <typename Ntype, typename T, int n, int m,
          std::enable_if_t<hila::is_complex<T>::value, int> = 0>
Matrix<n, m, Ntype> cast_to(const DaggerMatrix<n, m, T> &mat) {
    Matrix<n, m, Ntype> res;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++) res.e(i, j) = cast_to<Ntype>(mat.e(i, j));
    return res;
}

#endif
