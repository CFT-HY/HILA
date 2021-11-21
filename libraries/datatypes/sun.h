#ifndef SUN_M
#define SUN_M

#include "plumbing/defs.h"
#include "datatypes/cmplx.h"
#include "datatypes/matrix.h"
#include "plumbing/mersenne.h" //has to be included
#include <cmath>

// Macros for sped-up operations

#define CMULJJ(a, b, c)                                                                \
    do {                                                                               \
        (c).re = (a).re * (b).re - (a).im * (b).im;                                    \
        (c).im = -(a).re * (b).im - (a).im * (b).re;                                   \
    } while (0)
#define CMULJJ_ADD(a, b, c)                                                            \
    do {                                                                               \
        (c).re += (a).re * (b).re - (a).im * (b).im;                                   \
        (c).im += -(a).re * (b).im - (a).im * (b).re;                                  \
    } while (0)
#define CMULJJ_SUB(a, b, c)                                                            \
    do {                                                                               \
        (c).re -= (a).re * (b).re - (a).im * (b).im;                                   \
        (c).im += (a).re * (b).im + (a).im * (b).re;                                   \
    } while (0)

#define su3_vector_crossprod_conj(av, bv, res)                                         \
    do {                                                                               \
                                                                                       \
        CMULJJ((av).c[1], (bv).c[2], (res).c[0]);                                      \
        CMULJJ_SUB((av).c[2], (bv).c[1], (res).c[0]);                                  \
                                                                                       \
        CMULJJ((av).c[2], (bv).c[0], (res).c[1]);                                      \
        CMULJJ_SUB((av).c[0], (bv).c[2], (res).c[1]);                                  \
                                                                                       \
        CMULJJ((av).c[0], (bv).c[1], (res).c[2]);                                      \
        CMULJJ_SUB((av).c[1], (bv).c[0], (res).c[2]);                                  \
    } while (0)

// SU2 matrix multiplication routines ------------------------

#define nn_a(x, y) (x.d * y.a + x.a * y.d - x.b * y.c + x.c * y.b)
#define nn_b(x, y) (x.d * y.b + x.b * y.d - x.c * y.a + x.a * y.c)
#define nn_c(x, y) (x.d * y.c + x.c * y.d - x.a * y.b + x.b * y.a)
#define nn_d(x, y) (x.d * y.d - x.a * y.a - x.b * y.b - x.c * y.c)

#define na_a(x, y) (-x.d * y.a + x.a * y.d + x.b * y.c - x.c * y.b)
#define na_b(x, y) (-x.d * y.b + x.b * y.d + x.c * y.a - x.a * y.c)
#define na_c(x, y) (-x.d * y.c + x.c * y.d + x.a * y.b - x.b * y.a)
#define na_d(x, y) (x.d * y.d + x.a * y.a + x.b * y.b + x.c * y.c)

#define an_a(x, y) (x.d * y.a - x.a * y.d + x.b * y.c - x.c * y.b)
#define an_b(x, y) (x.d * y.b - x.b * y.d + x.c * y.a - x.a * y.c)
#define an_c(x, y) (x.d * y.c - x.c * y.d + x.a * y.b - x.b * y.a)
#define an_d(x, y) (x.d * y.d + x.a * y.a + x.b * y.b + x.c * y.c)

#define aa_a(x, y) (-x.d * y.a - x.a * y.d - x.b * y.c + x.c * y.b)
#define aa_b(x, y) (-x.d * y.b - x.b * y.d - x.c * y.a + x.a * y.c)
#define aa_c(x, y) (-x.d * y.c - x.c * y.d - x.a * y.b + x.b * y.a)
#define aa_d(x, y) (x.d * y.d - x.a * y.a - x.b * y.b - x.c * y.c)

//--------------------------------------------------------

///
/// SU(N) class
///
/// Specializes the matrix class with certain additional methods
/// specific to SU(N) matrices. Allows matching SU(N) vectors and
/// SU(N) vectors in multiplication.
///
template <int n, typename radix = double>
class SU : public SquareMatrix<n, Complex<radix>> {
  public:
    using base_type = hila::number_type<radix>;
    using argument_type = radix;

    static constexpr int size = n;

    /// Explicitly include default constructor. Necessary for use in fields.
    SU() = default;

    /// Construct from scalar by setting diagonal
    template <typename scalart,
              std::enable_if_t<hila::is_complex_or_arithmetic<scalart>::value, int> = 0>
    SU(const scalart rhs) {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                if (i == j)
                    this->e(i, j) = rhs;
                else
                    this->e(i, j) = (0);
            }
    }

    /// Construct from a square matrix
    template <typename S>
    SU(const SquareMatrix<n, Complex<S>> &m) {
        for (int i = 0; i < n * n; i++) {
            this->c[i] = m.c[i];
        }
    }

    /// Casting to different Complex type
    template <typename S>
    operator SquareMatrix<n, Complex<S>>() {
        SquareMatrix<n, Complex<radix>> r;
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) {
                r.e(i, j) = this->e(i, j);
            }
        return r;
    }

    /// Reunitarize the SU(N) matrix.  This should be OK when removing
    /// small deviations, but it is not a projective method which
    /// should be used when "unbiased" reunitarization is needed
    void reunitarize() {
        make_unitary();
        fix_det();
    };

    /// Mame the matrix unitary by orthogonalizing the rows
    /// There must be a faster way to do this, but this is simple
    ///  i = 0 ... n-1
    ///     normalize row i
    ///     make rows i+1 .. (n-1) orthogonal to row i
    void make_unitary() {

        for (int r = 0; r < n; r++) {

            // normalize row r
            radix n2 = 0;
            for (int c = 0; c < n; c++) n2 += this->e(r, c).squarenorm();
            n2 = 1.0 / sqrt(n2);
            for (int c = 0; c < n; c++) this->e(r, c) *= n2;

            // Now make rows r+1 .. n-1 orthogonal to row r,
            // by doing j = j - (r^* j) r

            Complex<radix> d;
            for (int j = r + 1; j < n; j++) {
                // dot productof r^* j
                d = 0;
                for (int i = 0; i < n : i++) {
                    d += this->e(r, i).conj_mul(this->e(j, i));
                }
                // and j -= d * r
                for (int i = 0; i < n : i++) {
                    this->e(j, i) -= d * this->e(r, i);
                }
            }
        }
    }

    /// Set the determinant of the SU(N) matrix to 1
    void fix_det() {
        Complex<radix> d, factor;
        radix t;
        int i, j;

        d = det(*(this));
        t = d.arg() / static_cast<radix>(n);
        factor = Complex<radix>(cos(-t), sin(-t));
        for (j = 0; j < n; j++)
            for (i = 0; i < n; i++) {
                this->e(j, i) = this->e(j, i) * factor;
            }
    }

    /// exp() - method for square matrices using Taylor expansion
    /// Evaluate as exp(H) = 1 + H + H^2/2! + H^3/3! ...

    // TODO:do this with eigevalues!
    SU exp(const int order = 20) const {

        SU A, Mp;
        A = *this;
        Mp = *this;
        radix multiplier;
        for (int k = 2; k <= order; k++) {
            multiplier = 1.0 / k;
            Mp *= *this;
            Mp *= multiplier;
            A += Mp;
        }
        A += 1;
        return A;
    }

    /// generate random SU(N) element by expanding exp(A), where A is a traceless
    /// hermitian matrix. more iterations are needed to generate larger elements: 12
    /// works well for n < 10.
    void random(const int depth = 12) {
        Matrix<n, n, Complex<radix>> A, An, res;
        An = 1;
        res = 1;
        Complex<radix> tr(1, 0), factor(1, 0);
        for (int i = 0; i < n; i++) {
            A.e(i, i) = Complex<radix>(hila::random(), 0.0);
            for (int j = 0; j < i; j++) {
                Complex<radix> a(static_cast<radix>(hila::random() / n),
                                 static_cast<radix>(hila::random() / n));
                A.e(i, j) = a;
                A.e(j, i) = a.conj();
            }
        }
        tr = A.trace() * (static_cast<radix>(1) / static_cast<radix>(n));
        for (int i = 0; i < n; i++) {
            A.e(i, i) -= tr;
        }
        An = A;
        for (int k = 1; k <= depth; k++) {
            factor = factor * Complex<radix>(0, 1) *
                     (static_cast<radix>(1) / static_cast<radix>(k));
            res += An * factor;
            An *= A;
        }
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                (*this).e(i, j) = res.e(i, j);
            }
    }

    /// Generate a random algebra (antihermitean) matrix
    void gaussian_algebra() {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                double a = hila::gaussrand();
                double b = hila::gaussrand();
                (*this).e(i, j).re = a;
                (*this).e(j, i).re = -a;
                (*this).e(i, j).im = b;
                (*this).e(j, i).im = b;
            }
        }

        for (int i = 0; i < n; i++) {
            (*this).e(i, i).re = 0;
            (*this).e(i, i).im = 0;
        }
        for (int i = 1; i < n; i++) {
            double a = hila::gaussrand() * sqrt(2.0 / (i * (i + 1)));
            for (int j = 0; j < i; j++) (*this).e(j, j).im += a;
            (*this).e(i, i).im -= i * a;
        }
    }

    /// The norm of an algebra matrix
    radix algebra_norm() {
        radix thissum = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                thissum += (*this).e(i, j).squarenorm();
            }
            radix diag = (*this).e(i, i).im;
            thissum += 0.5 * diag * diag;
        }
        return thissum;
    }




    /// find determinant using LU decomposition. Algorithm: numerical Recipes, 2nd ed.
    /// p. 47 ff
    Complex<radix> det_lu() {

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
            for (j = 0; j < n; j++) a[i][j] = this->e(i, j);

        for (i = 0; i < n; i++) {
            big = 0;
            for (j = 0; j < n; j++) {
                if ((temp = a[i][j].squarenorm()) > big) big = temp;
            }
            big = sqrt(big);
            if (big == 0.0) return 0;
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
                cdum = one / a[j][j]; // check Complex division
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

    static SU generator(int ng) {
        // SUN generators normalized as tr(T^2) = 1/2
        static bool initialize = true;
        static SU generators[n * n - 1];
        if (initialize)
            for (int g = 0; g < n * n - 1; g++) {
                generators[g] = 0;
                if (g < n - 1) {
                    // Diagonal generators
                    double w = 0.5 / sqrt((g + 1) * (g + 2) * 0.5);
                    for (int i = 0; i < g + 1; i++) {
                        generators[g].e(i, i).re = w;
                    }
                    generators[g].e(g + 1, g + 1).re = -(g + 1) * w;
                } else {
                    // Nondiagonal ones. Just run through the indexes and
                    // count until they match...
                    int k = n - 1;
                    for (int m1 = 0; m1 < n; m1++)
                        for (int m2 = m1 + 1; m2 < n; m2++) {
                            if (g == k) {
                                generators[g].e(m1, m2).re = 0.5;
                                generators[g].e(m2, m1).re = 0.5;
                            } else if (g == k + 1) {
                                generators[g].e(m1, m2).im = 0.5;
                                generators[g].e(m2, m1).im = -0.5;
                            }
                            k += 2;
                        }
                }
                initialize = false;
            }
        return generators[ng];
    }

    static constexpr int generator_count() {
        return n * n - 1;
    }
};

template <typename radix>
class SU2;

template <typename radix>
class SU3;

template <typename radix>
class SU2adjoint;

/// SU3<radix> is equivalent to SU<3, radix>
template <typename radix>
class SU3 : public SU<3, radix> {};

/// SU2 matrices are equivalent to SU<2, radix>, but are implemented more
/// efficiently here. This implementation represents the matrices as
/// U = d + a * sigma_1 + b * sigma_2 + c * sigma_3. This is in SU(2) if
/// a^2 + b^2 + c^2 + d^2 = 1.
template <typename radix>
class SU2 {
  public:
    using base_type = hila::number_type<radix>;
    using argument_type = radix;

    SU2() : a(0), b(0), c(0), d(1) {}

    SU2(radix *vals) : a(vals[0]), b(vals[1]), c(vals[2]), d(vals[3]) {
        normalize();
    }

    SU2(const SU2<radix> &rhs) {
        b = rhs.b;
        a = rhs.a;
        d = rhs.d;
        c = rhs.c;
    };

    friend SU2<radix> operator*(const SU2<radix> &x, const SU2<radix> &y) {
        SU2<radix> r;
        r.a = nn_a(x, y);
        r.b = nn_b(x, y);
        r.c = nn_c(x, y);
        r.d = nn_d(x, y);
        return r;
    }

    friend SU2<radix> operator*(const SU2<radix> &x, const SU2adjoint<radix> &y) {
        SU2<radix> r;
        r.a = na_a(x, y.ref);
        r.b = na_b(x, y.ref);
        r.c = na_c(x, y.ref);
        r.d = na_d(x, y.ref);
        return r;
    }

    friend SU2<radix> operator*(const SU2adjoint<radix> &x, const SU2<radix> &y) {
        SU2<radix> r;
        r.a = an_a(x.ref, y);
        r.b = an_b(x.ref, y);
        r.c = an_c(x.ref, y);
        r.d = an_d(x.ref, y);
        return r;
    }

    friend SU2<radix> operator*(const SU2adjoint<radix> &x,
                                const SU2adjoint<radix> &y) {
        SU2<radix> r;
        r.a = aa_a(x.ref, y.ref);
        r.b = aa_b(x.ref, y.ref);
        r.c = aa_c(x.ref, y.ref);
        r.d = aa_d(x.ref, y.ref);
        return r;
    }

    SU2<radix> &operator=(const SU2<radix> &);
    SU2<radix> &operator=(const SU2adjoint<radix> &);
    SU2<radix> &normalize();
    SU2<radix> &reunitarize();
    SU2<radix> &random();
    SU2<radix> &inv();
    SU2adjoint<radix> &adj();
    radix sqr() const;
    radix tr() const;
    radix det() const;
    SU2<radix> operator+(const SU2<radix> &); // basic operations. These versions return
                                              // new matrix as result
    SU2<radix> operator-(const SU2<radix> &);
    SU2<radix> &operator+=(
        const SU2<radix> &); // same ops as above, except store result in lhs matrix
    SU2<radix> &operator-=(const SU2<radix> &);
    SU2<radix> &operator*=(const SU2<radix> &);
    SU2<radix> &operator+=(const SU2adjoint<radix> &); // same ops as above, except
                                                       // store result in lhs matrix
    SU2<radix> &operator-=(const SU2adjoint<radix> &);
    SU2<radix> &operator*=(const SU2adjoint<radix> &);
    SU2<radix> &operator*=(const radix &);
    SU2<radix> operator*(const radix &);

  private:
    radix a, b, c, d;
};

/// An adjoint operation for SU2 matrices, which only contains a reference
/// to the original. Expands when the actual value of the matrix is needed
/// or an operation is applied.
template <typename radix>
class SU2adjoint {
  public:
    SU2adjoint(const SU2<radix> &rhs) : ref(rhs){};
    const SU2<radix> &ref;

  private:
    SU2adjoint() {}
    SU2adjoint<radix> &operator=(const SU2adjoint &rhs) {}
};

/// Conj returns the adjoing of the matrix
template <typename radix>
inline SU2adjoint<radix> conj(SU2<radix> &ref) {
    SU2adjoint<radix> result(ref);
    return result;
}

/// Take the adjoint. Returns an SU2adjoint
template <typename radix>
SU2adjoint<radix> &SU2<radix>::adj() {
    return SU2adjoint<radix>(*this);
};

template <typename radix>
radix SU2<radix>::sqr() const {
    return a * a + b * b + c * c + d * d;
}

template <typename radix>
radix SU2<radix>::det() const {
    return a * a + b * b + c * c + d * d;
}

template <typename radix>
radix SU2<radix>::tr() const {
    return 2 * d;
}

/// Apply normalization to make sure this is an
/// SU2 matrix
template <typename radix>
SU2<radix> &SU2<radix>::normalize() {
    radix sq = sqrt(this->sqr());
    a /= sq;
    b /= sq;
    c /= sq;
    d /= sq;
    return *this;
}

/// Apply normalization to make sure this is an
/// SU2 matrix
template <typename radix>
SU2<radix> &SU2<radix>::reunitarize() {
    return this->normalize();
}

/// Create a random SU2 matrix
template <typename radix>
SU2<radix> &SU2<radix>::random() {
    radix one, two;
    one = hila::gaussrand2(&two);
    a = one;
    b = two;
    one = hila::gaussrand2(&two);
    c = one;
    d = two;
    return this->normalize();
}

/// The inverse of an SU(N) matrix is its conjugate
template <typename radix>
SU2<radix> &SU2<radix>::inv() {
    a *= static_cast<radix>(-1);
    b *= static_cast<radix>(-1);
    c *= static_cast<radix>(-1);
    return *this;
}

template <typename radix>
SU2<radix> &SU2<radix>::operator=(const SU2<radix> &rhs) {
    a = rhs.a;
    b = rhs.b;
    c = rhs.c;
    d = rhs.d;
    return *this;
};

template <typename radix>
SU2<radix> &SU2<radix>::operator=(const SU2adjoint<radix> &rhs) {
    a = -rhs.a;
    b = -rhs.b;
    c = -rhs.c;
    d = rhs.d;
    return *this;
};

template <typename radix>
SU2<radix> &SU2<radix>::operator*=(const SU2<radix> &y) {
    a = nn_a((*this), y);
    b = nn_b((*this), y);
    c = nn_c((*this), y);
    d = nn_d((*this), y);
    return *this;
}

template <typename radix>
SU2<radix> &SU2<radix>::operator*=(const SU2adjoint<radix> &y) {
    a = na_a((*this), y);
    b = na_b((*this), y);
    c = na_c((*this), y);
    d = na_d((*this), y);
    return *this;
}

template <typename radix>
SU2<radix> &SU2<radix>::operator*=(const radix &rhs) {
    a *= rhs;
    b *= rhs;
    c *= rhs;
    d *= rhs;
    return *this;
};

template <typename radix>
SU2<radix> SU2<radix>::operator*(const radix &rhs) {
    SU2<radix> r;
    r.a = a * rhs;
    r.b = b * rhs;
    r.c = c * rhs;
    r.d = d * rhs;
    return r;
};

template <typename radix>
SU2<radix> SU2<radix>::operator+(const SU2<radix> &y) {
    SU2<radix> r;
    r.a = a + y.a;
    r.b = b + y.b;
    r.c = c + y.c;
    r.d = d + y.d;
    return r;
};

template <typename radix>
SU2<radix> SU2<radix>::operator-(const SU2<radix> &y) {
    SU2<radix> r;
    r.a = a - y.a;
    r.b = b - y.b;
    r.c = c - y.c;
    r.d = d - y.d;
    return r;
};

/// Project to the antihermitean part of a matrix
template <int N, typename radix>
void project_antihermitean(Matrix<N, N, Complex<radix>> &matrix) {
    radix tr = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i; j++) {
            radix a = 0.5 * (matrix.e(i, j).re - matrix.e(j, i).re);
            radix b = 0.5 * (matrix.e(i, j).im + matrix.e(j, i).im);
            matrix.e(i, j).re = a;
            matrix.e(j, i).re = -a;
            matrix.e(i, j).im = b;
            matrix.e(j, i).im = b;
        }
        tr += matrix.e(i, i).im;
        matrix.e(i, i).re = 0;
    }
    for (int i = 0; i < N; i++) {
        matrix.e(i, i).im -= tr / N;
    }
};

/// A new vector type for color vectors. These can be multiplied with appropriate SU(N)
/// matrices.
template <int n, typename radix>
class SU_vector : public Vector<n, Complex<radix>> {
  public:
    using base_type = hila::number_type<radix>;
    using argument_type = radix;

    static constexpr int size = n;

    SU_vector() = default;

    template <typename scalart,
              std::enable_if_t<hila::is_arithmetic<scalart>::value, int> = 0>
    SU_vector(const scalart rhs) {
        for (int i = 0; i < n; i++) {
            this->c[i] = (rhs);
        }
    }

    template <typename scalart,
              std::enable_if_t<hila::is_arithmetic<scalart>::value, int> = 0>
    SU_vector(const SU_vector<n, scalart> m) {
        for (int i = 0; i < n; i++) {
            this->c[i] = m.c[i];
        }
    }

    SU_vector(Vector<n, Complex<radix>> m) {
        for (int i = 0; i < n; i++) {
            this->c[i] = m.c[i];
        }
    }

    SU_vector operator-() const {
        SU_vector r;
        for (int i = 0; i < n; i++) {
            r.c[i] = -this->c[i];
        }
        return r;
    }

    inline radix rdot(const SU_vector &rhs) const {
        radix r = 0;
        for (int i = 0; i < n; i++) {
            r += this->c[i].re * rhs.c[i].re;
            r += this->c[i].im * rhs.c[i].im;
        }
        return r;
    }

    Matrix<n, n, Complex<radix>> outer_product(const SU_vector &rhs) const {
        Matrix<n, n, Complex<radix>> r;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                r.e(i, j) += this->c[i] * rhs.c[j];
            }
        return r;
    }
};

/// Why do I need a separate implementation for SU_vector?
template <int n, typename T>
inline SU_vector<n, T> operator*(const Matrix<n, n, T> &A, const SU_vector<n, T> &B) {
    SU_vector<n, T> res;

    for (int i = 0; i < n; i++) {
        res.e(i) = 0;
        for (int j = 0; j < n; j++) {
            res.e(i) += A.e(i, j) * B.e(j);
        }
    }
    return res;
}

#endif