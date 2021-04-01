#ifndef MATRIX_H_
#define MATRIX_H_

#include <type_traits>
#include <sstream>
#include "plumbing/defs.h"
#include "operations.h"
#include "datatypes/cmplx.h"

// Do this macro here to ease "switching" between using
// mul_add operation and the normal sum
// THE mul_add METHOD SEEMS TO BE SLOWER?
#define MUL_SUM(a, b, c) c += a * b
//#define MUL_SUM(a, b, c) c = mul_add(a, b, c)

template <const int n, const int m, typename T = double> class Array;

////////////////////////////////////////////////////////////////
/// The main nxm matrix type template
////////////////////////////////////////////////////////////////
template <const int n, const int m, typename T> class Matrix {
  public:
    static_assert(is_complex_or_arithmetic<T>::value,
                  "Matrix requires Complex or arithmetic type");

    /// std incantation for field types
    using base_type = number_type<T>;
    using argument_type = T;

    /// The data as a one dimensional array
    T c[n * m];

    /// define default constructors to ensure std::is_trivial
    Matrix() = default;
    ~Matrix() = default;
    Matrix(const Matrix<n, m, T> &v) = default;

    /// constructor from scalar -- keep it explicit!  Not good for auto use
    template <typename S, int nn = n, int mm = m,
              std::enable_if_t<(is_assignable<T &, S>::value && nn == mm), int> = 0>
    explicit inline Matrix(const S rhs) {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                if (i == j)
                    e(i, j) = rhs;
                else
                    e(i, j) = 0;
            }
    }

    /// construct from 0
    inline Matrix(const std::nullptr_t &z) {
        for (int i = 0; i < n * m; i++) {
            c[i] = 0;
        }
    }

    /// Construct matrix automatically from right-size initializer list
    /// This does not seem to be dangerous, so keep non-explicit

    template <typename S, std::enable_if_t<is_assignable<T &, S>::value, int> = 0>
    inline Matrix(std::initializer_list<S> rhs) {
        assert(rhs.size() == n * m &&
               "Matrix/Vector initializer list size must match variable size");
        int i = 0;
        for (auto it = rhs.begin(); it != rhs.end(); it++, i++) {
            c[i] = *it;
        }
    }

    /// Define constant methods rows(), columns() - may be useful in template code
    static constexpr int rows() { return n; }
    static constexpr int columns() { return m; }

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
    inline T e(const int i, const int j) const { return c[i * m + j]; }
    /// standard access ops m.e(i,j) - assume T is small, as it should
    inline T &e(const int i, const int j) { return c[i * m + j]; }

    /// declare single e here too in case we have a vector
    /// (one size == 1)
    template <int q = n, int p = m, std::enable_if_t<(q == 1 || p == 1), int> = 0>
    inline T e(const int i) const {
        return c[i];
    }

    template <int q = n, int p = m, std::enable_if_t<(q == 1 || p == 1), int> = 0>
    inline T &e(const int i) {
        return c[i];
    }

    /// And also [] for vectors (not matrices!)
    /// (one size == 1)
    template <int q = n, int p = m, std::enable_if_t<(q == 1 || p == 1), int> = 0>
    inline T operator[](const int i) const {
        return c[i];
    }

    template <int q = n, int p = m, std::enable_if_t<(q == 1 || p == 1), int> = 0>
    inline T &operator[](const int i) {
        return c[i];
    }

    /// interpret Matrix as Array -  for array ops
    Array<n, m, T> &asArray() { return *reinterpret_cast<Array<n, m, T> *>(this); }
    const Array<n, m, T> asArray() const {
        return *reinterpret_cast<const Array<n, m, T> *>(this);
    }

    /// casting from one Matrix (number) type to another: do not do this automatically.
    /// but require an explicit cast operator.  This makes it easier to write code.
    /// or should it be automatic?  keep/remove explicit?
    /// TODO: CHECK AVX CONVERSIONS

    template <typename S, std::enable_if_t<std::is_convertible<T, S>::value, int> = 0>
    operator Matrix<n, m, S>() {
        Matrix<n, m, S> res;
        for (int i = 0; i < n * m; i++)
            res.c[i] = c[i];
        return res;
    }

    /// unary -
    inline Matrix<n, m, T> operator-() const {
        Matrix<n, m, T> res;
        for (int i = 0; i < n * m; i++) {
            res.c[i] = -c[i];
        }
        return res;
    }

    /// unary +
    inline Matrix<n, m, T> operator+() const { return *this; }

    /// assign from 0
    inline Matrix<n, m, T> &operator=(const std::nullptr_t &z) {
        for (int i = 0; i < n * m; i++) {
            c[i] = 0;
        }
        return *this;
    }

    /// Assign from different type matrix
    template <typename S, std::enable_if_t<is_assignable<T &, S>::value, int> = 0>
    inline Matrix<n, m, T> &operator=(const Matrix<n, m, S> &rhs) {
        for (int i = 0; i < n * m; i++) {
            c[i] = rhs.c[i];
        }
        return *this;
    }

    /// Assign from "scalar" for square matrix
    template <typename S, int nn = n, int mm = m,
              std::enable_if_t<(is_assignable<T &, S>::value && nn == mm), int> = 0>
    inline Matrix<n, m, T> &operator=(const S rhs) {
        // static_assert( n==m, "rows != columns : assigning a scalar only possible for
        // a square Matrix");

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
    template <typename S, std::enable_if_t<is_assignable<T &, S>::value, int> = 0>
    inline Matrix<n, m, T> &operator=(std::initializer_list<S> rhs) {
        assert(rhs.size() == n * m &&
               "Initializer list has a wrong size in assignment");
        int i = 0;
        for (auto it = rhs.begin(); it != rhs.end(); it++, i++) {
            c[i] = *it;
        }
        return *this;
    }

    /// add assign a Matrix
#pragma hila loop_function
    template <typename S, std::enable_if_t<std::is_assignable<T &, S>::value, int> = 0>
    Matrix<n, m, T> &operator+=(const Matrix<n, m, S> &rhs) {
        for (int i = 0; i < n * m; i++) {
            c[i] += rhs.c[i];
        }
        return *this;
    }

    /// subtract assign a Matrix
    template <typename S, std::enable_if_t<std::is_assignable<T &, S>::value, int> = 0>
    Matrix<n, m, T> &operator-=(const Matrix<n, m, S> &rhs) {
        for (int i = 0; i < n * m; i++) {
            c[i] -= rhs.c[i];
        }
        return *this;
    }

    /// add assign type T and convertible
    template <typename S, std::enable_if_t<
                              std::is_assignable<T &, type_plus<T, S>>::value, int> = 0>
    Matrix<n, m, T> &operator+=(const S rhs) {
        static_assert(
            n == m,
            "rows != columns : scalar addition possible for square matrices only!");
        for (int i = 0; i < n; i++) {
            e(i, i) += rhs;
        }
        return *this;
    }

    /// subtract assign type T and convertible
    template <
        typename S,
        std::enable_if_t<std::is_assignable<T &, type_minus<T, S>>::value, int> = 0>
    Matrix<n, m, T> &operator-=(const S rhs) {
        static_assert(
            n == m,
            "rows != columns : scalar subtraction possible for square matrices only!");
        for (int i = 0; i < n; i++) {
            e(i, i) -= rhs;
        }
        return *this;
    }

    /// multiply assign with matrix
    template <int p, typename S,
              std::enable_if_t<std::is_assignable<T &, type_mul<T, S>>::value, int> = 0>
    Matrix<n, m, T> &operator*=(const Matrix<m, p, S> &rhs) {
        static_assert(m == p,
                      "can't assign result of *= to lhs Matrix, because doing so "
                      "would change it's dimensions");
        *this = *this * rhs;
        return *this;
    }

    /// multiply assign with scalar
    template <typename S,
              std::enable_if_t<std::is_assignable<T &, type_mul<T, S>>::value, int> = 0>
    Matrix<n, m, T> &operator*=(const S rhs) {
        for (int i = 0; i < n * m; i++) {
            c[i] *= rhs;
        }
        return *this;
    }

    /// divide assign with scalar
    template <typename S,
              std::enable_if_t<std::is_assignable<T &, type_div<T, S>>::value, int> = 0>
    Matrix<n, m, T> &operator/=(const S rhs) {
        for (int i = 0; i < n * m; i++) {
            c[i] /= rhs;
        }
        return *this;
    }

    /// numpy style matrix fill
    template <typename S, std::enable_if_t<is_assignable<T &, S>::value, int> = 0>
    Matrix<n, m, T> &fill(const S rhs) output_only {
        for (int i = 0; i < n * m; i++)
            c[i] = rhs;
        return *this;
    }

    // return copy of transpose of this matrix
    inline Matrix<m, n, T> transpose() const {
        Matrix<m, n, T> res;
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++) {
                res.e(i, j) = e(j, i);
            }
        return res;
    }

    /// and complex conjugate
    inline Matrix<n, m, T> conj() const {
        Matrix<n, m, T> res;
        for (int i = 0; i < n * m; i++) {
            res.c[i] = ::conj(c[i]);
        }
        return res;
    }

    /// return copy of adjoint/dagger of this Matrix
    inline Matrix<m, n, T> adjoint() const {
        Matrix<m, n, T> res;
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++) {
                res.e(i, j) = ::conj(e(j, i));
            }
        return res;
    }

    /// alias dagger() to adjoint()
    inline Matrix<m, n, T> dagger() const { return adjoint(); }

    /// return real part
    inline Matrix<n, m, number_type<T>> real() const {
        Matrix<n, m, number_type<T>> res;
        for (int i = 0; i < m * n; i++) {
            res.c[i] = real(c[i]);
        }
        return res;
    }

    /// return imaginary part
    inline Matrix<n, m, number_type<T>> imag() const {
        Matrix<n, m, number_type<T>> res;
        for (int i = 0; i < m * n; i++) {
            res.c[i] = imag(c[i]);
        }
        return res;
    }

    /// matrix trace (return type T)
    T trace() const {
        static_assert(n == m, "trace not defined for non square matrices!");
        T result = 9;
        for (int i = 0; i < n; i++) {
            result += e(i, i);
        }
        return result;
    }

    /// calculate square norm - sum of squared elements
    number_type<T> norm_sq() const {
        number_type<T> result = 0;
        for (int i = 0; i < n * m; i++) {
            result += norm_squared(c[i]);
        }
        return result;
    }

    /// dot product - vector . vector
    template <int p, std::enable_if_t<(p == 1 && m == 1), int> = 0>
    inline T dot(const Matrix<n, p, T> &rhs) const {
        T r = 0;
        for (int i = 0; i < n; i++) {
            r += ::conj(c[i]) * rhs.c[i];
        }
        return r;
    }

    /// dot with matrix - matrix
    template <int p, std::enable_if_t<(p > 1 || m > 1), int> = 0>
    inline Matrix<m, p, T> dot(const Matrix<n, p, T> &rhs) const {
        Matrix<m, p, T> res;
        for (int i = 0; i < m; i++)
            for (int j = 0; j < p; j++) {
                res.e(i, j) = 0;
                for (int k = 0; k < n; j++)
                    res.e(i, j) += ::conj(e(k, i)) * rhs.e(k, j);
            }
    }

    /// Generate random elements
    Matrix<n, m, T> &random() output_only {
        for (int i = 0; i < n * m; i++) {
            ::random(c[i]);
        }
        return *this;
    }

    /// Generate gaussian random elements
    inline Matrix<n, m, T> &gaussian() output_only {
        for (int i = 0; i < n * m; i++) {
            ::gaussian_random(c[i]);
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

/// do transpose and adjoint functions here
template <const int n, const int m, typename T>
inline Matrix<n, m, T> transpose(const Matrix<m, n, T> &arg) {
    return arg.transpose();
}
/// conjugate
template <const int n, const int m, typename T>
inline Matrix<n, m, T> conj(const Matrix<n, m, T> &arg) {
    return arg.conj();
}
/// adjoint (=conjugate)
template <const int n, const int m, typename T>
inline Matrix<n, m, T> adjoint(const Matrix<m, n, T> &arg) {
    return arg.adjoint();
}
/// dagger (=conjugate)
template <const int n, const int m, typename T>
inline Matrix<n, m, T> dagger(const Matrix<m, n, T> &arg) {
    return arg.adjoint();
}
/// trace
template <const int n, const int m, typename T>
inline T trace(const Matrix<n, m, T> &arg) {
    return arg.trace();
}
/// real part
template <const int n, const int m, typename T>
inline Matrix<n, m, number_type<T>> real(const Matrix<n, m, T> &arg) {
    return arg.real();
}
/// imaginary part
template <const int n, const int m, typename T>
inline Matrix<n, m, number_type<T>> imag(const Matrix<n, m, T> &arg) {
    return arg.imag();
}

/// templates needed for naive calculation of determinants
template <int n, int m, typename T>
Matrix<n - 1, m - 1, T> Minor(const Matrix<n, m, T> &bigger, int i, int j) {
    Matrix<n - 1, m - 1, T> result;
    int ri = 0, bi = 0;
    for (int p = 0; p < n; p++)
        for (int l = 0; l < m; l++, bi++) {
            if (p != i && l != j) {
                result.c[ri] = bigger.c[bi];
                ri++;
            }
        }
    return result;
}

/// determinant -> use LU factorization later
template <int n, int m, typename T> T det(const Matrix<n, m, T> &mat) {
    static_assert(n == m, "determinants defined only for square matrices");
    T result = 0;
    number_type<T> parity = 1, opposite = -1;
    for (int i = 0; i < n; i++) {
        Matrix<n - 1, m - 1, T> minor = Minor(mat, 0, i);
        result += parity * det(minor) * mat.e(0, i);
        parity *= opposite;
    }
    return result;
}

/// 2x2 matrix det
template <typename T> T det(const Matrix<2, 2, T> &mat) {
    return mat.e(0, 0) * mat.e(1, 1) - mat.e(1, 0) * mat.e(0, 1);
}

/// 1x1 matrix det (trivial)
template <typename T> T det(const Matrix<1, 1, T> &mat) { return mat.e(0, 0); }

/// Now matrix additions: matrix + matrix
template <int n, int m, typename T>
inline Matrix<n, m, T> operator+(Matrix<n, m, T> a, const Matrix<n, m, T> &b) {
    a += b;
    return a;
}

/// Matrix subtract
template <int n, int m, typename T>
inline Matrix<n, m, T> operator-(Matrix<n, m, T> a, const Matrix<n, m, T> &b) {
    a -= b;
    return a;
}

/// Matrix + scalar
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_plus<T, S>, T>::value, int> = 0>
inline Matrix<n, m, T> operator+(Matrix<n, m, T> a, const S b) {
    a += b;
    return a;
}

/// scalar + matrix
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_plus<T, S>, T>::value, int> = 0>
inline Matrix<n, m, T> operator+(const S b, Matrix<n, m, T> a) {
    a += b;
    return a;
}

/// matrix - scalar
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_minus<T, S>, T>::value, int> = 0>
Matrix<n, m, T> operator-(Matrix<n, m, T> a, const S b) {
    a -= b;
    return a;
}

/// scalar - matrix
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_minus<S, T>, T>::value, int> = 0>
inline Matrix<n, m, T> operator-(const S b, Matrix<n, m, T> a) {
    static_assert(
        n == m,
        "rows != columns : scalar subtraction possible for square matrices only!");
    for (int i = 0; i < n; i++)
        a.e(i, i) = b - a.e(i, i);
    return a;
}

////////////////////////////////////////
/// matrix * matrix is the most important bit
template <int n, int m, int p, typename T>
inline Matrix<n, p, T> operator*(const Matrix<n, m, T> &A, const Matrix<m, p, T> &B) {
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
template <int m, typename T>
T operator*(const Matrix<1, m, T> &A, const Matrix<m, 1, T> &B) {
    T res = 0;
    for (int i = 0; i < m; i++) {
        res += A.e(i) * B.e(i);
    }
    return res;
}

/// matrix * scalar
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_mul<T, S>, T>::value, int> = 0>
Matrix<n, m, T> operator*(Matrix<n, m, T> mat, const S rhs) {
    mat *= rhs;
    return mat;
}

/// scalar * matrix
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_mul<S, T>, T>::value, int> = 0>
Matrix<n, m, T> operator*(const S rhs, Matrix<n, m, T> mat) {
    mat *= rhs; // assume commutativity, should be ok
    return mat;
}

/// matrix / scalar
template <int n, int m, typename T, typename S,
          std::enable_if_t<std::is_convertible<type_div<T, S>, T>::value, int> = 0>
Matrix<n, m, T> operator/(Matrix<n, m, T> mat, const S rhs) {
    mat /= rhs;
    return mat;
}

/// Stream operator
template <int n, int m, typename T>
std::ostream &operator<<(std::ostream &strm, const Matrix<n, m, T> &A) {
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
template <int n, int m, typename T>
inline number_type<T> norm_squared(const Matrix<n, m, T> &rhs) {
    return rhs.norm_sq();
}

/// Function that calls random()-method
template <int n, int m, typename T>
inline void random(output_only Matrix<n, m, T> &mat) {
    mat.random();
}

/// Function that calls the gaussian()-method
template <int n, int m, typename T>
inline void gaussian_random(output_only Matrix<n, m, T> &mat) {
    mat.gaussian();
}

/// find determinant using LU decomposition. Algorithm: numerical Recipes, 2nd ed. p. 47
/// ff
template <int n, typename T, typename radix = number_type<T>,
          std::enable_if_t<is_arithmetic<radix>::value, int> = 0,
          std::enable_if_t<std::is_same<T, Complex<radix>>::value, int> = 0>
Complex<radix> det_lu(const Matrix<n, n, T> &mat) {

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

/// find determinant using LU decomposition. Algorithm: numerical Recipes, 2nd ed. p. 47
/// ff
template <int n, typename T, std::enable_if_t<is_arithmetic<T>::value, int> = 0>
T det_lu(const Matrix<n, n, T> &mat) {
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

///////////////////////////////////////////////////////////////
/// Define Vector as a special case of Matrix
///////////////////////////////////////////////////////////////
template <int n, typename T> using Vector = Matrix<n, 1, T>;

///////////////////////////////////////////////////////////////
/// Define HorizontalVector as a special case of Matrix
///////////////////////////////////////////////////////////////
template <int n, typename T> using HorizontalVector = Matrix<1, n, T>;

///////////////////////////////////////////////////////////////
/// Define SquareMatrix as a special case of Matrix
///////////////////////////////////////////////////////////////
template <int n, typename T> using SquareMatrix = Matrix<n, n, T>;

#include "datatypes/array.h"

#endif
