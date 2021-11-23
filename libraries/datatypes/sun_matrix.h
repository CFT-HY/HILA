#ifndef SUN_MATRIX_H_
#define SUN_MATRIX_H_

#include "matrix.h"

/// Define type SUmatrix<n,type>
/// Derives from square Matrix<Complex> type

template <typename G>
class Algebra;

template <int N, typename T>
class SUmatrix : public Matrix_t<N, N, Complex<T>, SUmatrix<N,T>> {

  public:
    // std incantation for field types
    using base_type = hila::number_type<T>;
    using argument_type = T;

    // get all constructors from base
    using Matrix_t<N,N,Complex<T>,SUmatrix<N,T>>::Matrix_t;
    using Matrix_t<N,N,Complex<T>,SUmatrix<N,T>>::operator=;
    using Matrix_t<N,N,Complex<T>,SUmatrix<N,T>>::operator+=;
    using Matrix_t<N,N,Complex<T>,SUmatrix<N,T>>::operator-=;
    using Matrix_t<N,N,Complex<T>,SUmatrix<N,T>>::operator*=;
    using Matrix_t<N,N,Complex<T>,SUmatrix<N,T>>::operator/=;
    
    /// Make the matrix unitary by orthogonalizing the rows
    /// There must be a faster way to do this, but this is simple
    ///  i = 0 ... n-1
    ///     normalize row i
    ///     make rows i+1 .. (n-1) orthogonal to row i
    SUmatrix &make_unitary() {

        for (int r = 0; r < N; r++) {

            // normalize row r
            T n2 = 0;
            // use here function instead of method, works for double/float too
            for (int c = 0; c < N; c++)
                n2 += ::squarenorm(this->e(r, c));
            n2 = 1.0 / sqrt(n2);
            for (int c = 0; c < N; c++)
                this->e(r, c) *= n2;

            // Now make rows r+1 .. n-1 orthogonal to row r,
            // by doing j = j - (r^* j) r

            Complex<T> d;
            for (int j = r + 1; j < N; j++) {
                // dot productof r^* j
                d = 0;
                for (int i = 0; i < N; i++) {
                    d += ::conj(this->e(r, i)) * this->e(j, i);
                }
                // and j -= d * r
                for (int i = 0; i < N; i++) {
                    this->e(j, i) -= d * this->e(r, i);
                }
            }
        }
        return *this;
    }

    /// Set the determinant of the SU(N) matrix to 1
    inline SUmatrix &fix_det() {

        Complex<T> d, factor;
        T t;

        d = det(*(this));
        t = d.arg() / N;
        factor = Complex<T>(cos(-t), sin(-t));
        this->asArray() *= factor;
        return *this;
    }

    /// Make the matrix special unitary
    inline SUmatrix &make_group_element() {
        make_unitary();
        fix_det();
        return *this;
    }

    /// Project matrix to antihermitean and traceless algebra
    /// of the group.
    Algebra<SUmatrix<N, T>> project_to_algebra() const {
        Algebra<SUmatrix<N, T>> a;
        T tr = 0;

        // get first the traceless and imag diag part
        for (int i = 0; i < N - 1; i++)
            tr += a.diag[i] = this->e(i, i).im;
        tr += this->e(N - 1, N - 1).im;
        tr /= N;
        for (int i = 0; i < N - 1; i++)
            a.diag[i] -= tr;
        // Then off-diag bits

        int k = 0;
        for (int i = 0; i < N - 1; i++)
            for (int j = i + 1; j < N; j++)
                a.offdiag[k++] = (this->e(i, j) - this->e(j, i).conj()) * 0.5;

        return a;
    }
};

///////////////////////////////////////////////////////////
/// Specialize Algebra type to SU(N)

template <int N, typename T>
class Algebra<SUmatrix<N, T>> {
  public:
    // std incantation for field types
    using base_type = hila::number_type<T>;
    using argument_type = T;

    // storage for the diagonal and off-diag
    // components of the antihermitean traceless matrix
    static constexpr int n_offdiag = N * (N - 1) / 2;
    static constexpr int n_diag = N - 1;

    Complex<T> offdiag[n_offdiag];
    T diag[n_diag];

    /// std constructors and delete
    Algebra() = default;
    ~Algebra() = default;
    Algebra(const Algebra &v) = default;

    /// construct from 0
    inline Algebra(const std::nullptr_t &z) {
        for (int i = 0; i < n_offdiag; i++) {
            offdiag[i] = 0;
        }
        for (int i = 0; i < n_diag; i++) {
            diag[i] = 0;
        }
    }

    /// assignments
    inline Algebra &operator=(const Algebra &v) = default;

    template <typename S>
    inline Algebra &operator=(const Algebra<SUmatrix<N, S>> &v) {
        for (int i = 0; i < n_offdiag; i++) {
            offdiag[i] = v.offdiag[i];
        }
        for (int i = 0; i < n_diag; i++) {
            diag[i] = v.diag[i];
        }
        return *this;
    }

    /// from 0
    inline Algebra &operator=(const std::nullptr_t n) {
        for (int i = 0; i < n_offdiag; i++) {
            offdiag[i] = 0;
        }
        for (int i = 0; i < n_diag; i++) {
            diag[i] = 0;
        }
        return *this;
    }

    /// add assign
    template <typename S>
    inline Algebra &operator+=(const Algebra<SUmatrix<N, S>> &v) {
        for (int i = 0; i < n_offdiag; i++) {
            offdiag[i] += v.offdiag[i];
        }
        for (int i = 0; i < n_diag; i++) {
            diag[i] += v.diag[i];
        }
        return *this;
    }

    /// and sub
    template <typename S>
    inline Algebra &operator-=(const Algebra<SUmatrix<N, S>> &v) {
        for (int i = 0; i < n_offdiag; i++) {
            offdiag[i] -= v.offdiag[i];
        }
        for (int i = 0; i < n_diag; i++) {
            diag[i] -= v.diag[i];
        }
        return *this;
    }

    /// multiply by real scalar
    template <typename S, std::enable_if_t<hila::is_arithmetic<S>::value, int> = 0>
    inline Algebra &operator*=(const S &rhs) {
        for (int i = 0; i < n_offdiag; i++) {
            offdiag[i] *= rhs;
        }
        for (int i = 0; i < n_diag; i++) {
            diag[i] *= rhs;
        }
        return *this;
    }

    /// Unary - and +
    inline Algebra operator-() const {
        Algebra res;
        for (int i = 0; i < n_offdiag; i++)
            res.offdiag[i] = -offdiag[i];
        for (int i = 0; i < n_diag; i++) {
            res.diag[i] = -diag[i];
        }
        return res;
    }

    /// unary +
    inline const Algebra &operator+() const {
        return *this;
    }

    /// expand algebra to matrix rep - antihermitean
    SUmatrix<N, T> expand() const {
        SUmatrix<N, T> m;
        int k = 0;
        T trace = 0.0;
        for (int i = 0; i < N - 1; i++) {
            trace += diag[i];
            m.e(i, i) = Complex<T>((T)0, diag[i]);
            for (int j = i + 1; j < N; j++) {
                m.e(i, j) = offdiag[k];
                m.e(j, i) = -offdiag[k].conj();
                k++;
            }
        }
        m.e(N - 1, N - 1) = Complex<T>((T)0, -trace);

        return m;
    }

    // suN generators, normalized as
    //  Tr(\lambda_i\lambda_j) = 1/2 \delta_ij
    // off-diagonal are just that:
    //    \lambda^od_ij,r = 1/2 for elements ij and ji
    //    \lambda^od_ij,i = i/2 for ij; -i for ji
    //
    // diagonals: su(N) has N-1 diag generators
    // parametrize these recursively:
    //  \lambda_1 = diag(1,-1,0,0,..)/sqrt(1)/2
    //  \lambda_2 = diag(1,1,-2,0,..)/sqrt(3)/2
    //  \lambda_3 = diag(1,1,1,-3,..)/sqrt(6)/2
    //  \lambda_i = diag(1,.. ,-i,..)/sqrt(i(i+1)/2)/2
    //  ..
    //  \lambda_N-1 = diag(1,.. ,1,-(N-1))/sqrt(N(N-1)/2)/2
    //
    // Normalize so that the coefficients of \lambda_i are
    // set by gaussrand(), i.e. their <> = 1.  Then
    //     ah = i h = i \lambda_i \xi_i,
    // and 
    //     p = exp(-Tr h^2 ) = exp(-1/2 \xi_i^2)
    //
    //   <|od|^2> = 1/2 = <od.re^2> + <od.im^2> = 1/4 + 1/4
    //

    Algebra &gaussian_random() output_only {

        const T inv2 = 1.0 / 2.0;

        // off-diag elements, easy
        for (int i = 0; i < n_offdiag; i++) {
            T r;
            offdiag[i].re = inv2 * hila::gaussrand2(r);
            offdiag[i].im = inv2 * r;
        }

        // then, diagonal elements - sum over generators i

        diag[0] = 0;
        for (int i = 1; i < N; i++) {
            // sqrt(0.125) = sqrt(1/2)/2 
            T r = hila::gaussrand() * sqrt(0.125 / (i * (i + 1)));
            // add the normalized 'ones' in generators
            for (int j = 0; j < i; j++)
                diag[j] += r;
            // and set the negative element - no sum here, first contrib
            if (i < N - 1)
                diag[i] = -i * r;
        }

        return *this;
    }

    /// calculate the norm as sum of squares of elements of the coeffs of
    /// the generators:  if h = a_i \lambda_i, gives a_i^2
    /// Thus, 2 Tr h^2 = a_i^2
    /// Or, sum |h_ij|^2 = sum h_ji^* h_ij = Tr h^2 = 1/2 a_i^2

    T squarenorm() const {
        T res = 0;
        T rem = 0;
        for (int i = 0; i < n_diag; i++) {
            rem += diag[i];
            res += 2.0 * diag[i] * diag[i];
        }
        res += 2.0 * rem * rem;
        for (int i = 0; i < n_offdiag; i++) {
            res += 4.0 * offdiag[i].squarenorm();
        }
        return res;
    }
};

template <int N, typename T>
SUmatrix<N, T> exp(const Algebra<SUmatrix<N, T>> &a) {
    SUmatrix<N, T> m = a.expand();
    return exp(m);
}

#endif