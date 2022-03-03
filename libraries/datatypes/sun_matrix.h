#ifndef SUN_MATRIX_H_
#define SUN_MATRIX_H_

#include "matrix.h"

/// Define type SUmatrix<n,type>
/// Derives from square Matrix<Complex> type

template <typename G>
class Algebra;

template <int N, typename T>
class SUmatrix : public Matrix_t<N, N, Complex<T>, SUmatrix<N, T>> {

  public:
    // std incantation for field types
    using base_type = hila::number_type<T>;
    using argument_type = T;

    // get all constructors from base
    using Matrix_t<N, N, Complex<T>, SUmatrix<N, T>>::Matrix_t;
    using Matrix_t<N, N, Complex<T>, SUmatrix<N, T>>::operator=;
    using Matrix_t<N, N, Complex<T>, SUmatrix<N, T>>::operator+=;
    using Matrix_t<N, N, Complex<T>, SUmatrix<N, T>>::operator-=;
    using Matrix_t<N, N, Complex<T>, SUmatrix<N, T>>::operator*=;
    using Matrix_t<N, N, Complex<T>, SUmatrix<N, T>>::operator/=;

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
    inline SUmatrix &reunitarize() {
        make_unitary();
        fix_det();
        return *this;
    }

    /// multiply matrix rows i,j by SU2 "subgroup" from left
    void mult_by_SU2_left(int r, int q, const SUmatrix<2, T> &m) {
        // copy 2xN matrix
        Vector<2, Complex<T>> a;
        for (int i = 0; i < N; i++) {
            a.e(0) = this->e(r, i);
            a.e(1) = this->e(q, i);

            a = m * a;

            this->e(r, i) = a.e(0);
            this->e(q, i) = a.e(1);
        }
    }

    SUmatrix &random(int nhits = 16) out_only {

        // use Pauli matrix representation to generate SU(2) random matrix
        if constexpr (N == 2) {
            Vector<4, T> v;
            v.gaussian_random();
            v /= v.norm();
            this->e(0, 0) = Complex<T>(v[0], v[3]);
            this->e(1, 1) = Complex<T>(v[0], -v[3]);
            this->e(0, 1) = Complex<T>(v[2], v[1]);
            this->e(1, 0) = Complex<T>(-v[2], v[1]);

        } else {

            *this = 1;
            SUmatrix<2, T> m2;

            for (int h = 1; h <= nhits; h++) {
                for (int r = 0; r < N - 1; r++)
                    for (int q = r + 1; q < N; q++) {
                        m2.random();
                        this->mult_by_SU2_left(r, q, m2);
                    }

                // keep it SU(N)
                if (h % 16 == 0) {
                    this->reunitarize();
                }
            }
            if (nhits % 16 != 0)
                this->reunitarize();
        }

        return *this;
    }


    /// Project matrix to antihermitean and traceless algebra
    /// of the group.
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
    // Define \lambda's so that diagonals come first
    //
    // Dividing U = U_ah + U_h, U_ah = 1/2 (U - U^+) = i a_i \lambda_i + tr.im I/N
    // =>  Tr \lambda_i (U_ah) = 1/2 i a_i = 1/2 (\lambda_i)_jk (u_kj - u_jk^*)
    // =>  a_i = -i (\lambda_i)_jk (u_kj - u_jk^*)
    //
    // Thus, for diags,
    //     a_i = (u_00 + ..u_(i-1)(i-1) - i*u_ii).im 2/(sqrt(2i(i+1)))
    //
    // and off-diags:
    // symm: a_i = -i (u_kj - u_kj^* + u_jk - u_jk^*)/2 = -i i(u_kj.im + u_jk.im)
    //           = (u_kj.im + u_jk.im)
    // antisymm:  a_i = -i (i u_kj - i u_jk^* - i u_jk + i u_kj^*)/2
    //                = (u_kj.re - u_jk.re)

    Algebra<SUmatrix<N, T>> project_to_algebra() const {
        Algebra<SUmatrix<N, T>> a;

        // diagonal generators
        T sum = this->e(0, 0).im;
        for (int i = 1; i < N; i++) {
            a.e(i - 1) = (sum - i * this->e(i, i).im) / sqrt(0.5 * i * (i + 1));
            sum += this->e(i, i).im;
        }

        // Then off-diag bits
        int k = a.n_diag;
        for (int i = 0; i < N - 1; i++) {
            for (int j = i + 1; j < N; j++) {
                auto od = this->e(i, j) - this->e(j, i).conj();
                a.e(k) = od.re;
                a.e(k + 1) = od.im;
                k += 2;
            }
        }

        return a;
    }
};

///////////////////////////////////////////////////////////
/// Specialize Algebra type to SU(N)
/// Derive from (real) Vector of N*N-1 elements

template <int N, typename T>
class Algebra<SUmatrix<N, T>>
    : public Matrix_t<N * N - 1, 1, T, Algebra<SUmatrix<N, T>>> {
  public:
    // std incantation for field types
    using base_type = hila::number_type<T>;
    using argument_type = T;

    // storage for the diagonal and off-diag
    // components of the antihermitean traceless matrix
    static constexpr int n_offdiag = N * (N - 1);
    static constexpr int n_diag = N - 1;
    static constexpr int N_a = N * N - 1;

    /// std constructors and operators derived from vector
    using Matrix_t<N * N - 1, 1, T, Algebra<SUmatrix<N, T>>>::Matrix_t;
    using Matrix_t<N * N - 1, 1, T, Algebra<SUmatrix<N, T>>>::operator=;
    using Matrix_t<N * N - 1, 1, T, Algebra<SUmatrix<N, T>>>::operator+=;
    using Matrix_t<N * N - 1, 1, T, Algebra<SUmatrix<N, T>>>::operator-=;
    using Matrix_t<N * N - 1, 1, T, Algebra<SUmatrix<N, T>>>::operator*=;
    using Matrix_t<N * N - 1, 1, T, Algebra<SUmatrix<N, T>>>::operator/=;

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
    // Define \lambda's so that diagonals come first

    /// expand algebra to matrix rep - antihermitean
    SUmatrix<N, T> expand() const {
        SUmatrix<N, T> m;

        Vector<N, T> d;

        d.e(0) = (T)0;
        for (int i = 1; i < N; i++) {
            T r = this->e(i - 1) * sqrt(0.5 / (i * (i + 1)));
            // the contributions from 1's above
            for (int j = 0; j < i; j++)
                d.e(j) += r;
            // and set the negative element - no sum here, first contrib
            d.e(i) = -i * r;
        }

        for (int i = 0; i < N; i++)
            m.e(i, i) = Complex<T>(0, d.e(i));

        int k = n_diag;
        T inv2 = 1.0 / 2;
        for (int i = 0; i < N - 1; i++)
            for (int j = i + 1; j < N; j++) {
                Complex<T> v(this->c[k] * inv2, this->c[k + 1] * inv2);
                m.e(i, j) = v;
                m.e(j, i) = -v.conj();
                k += 2;
            }

        return m;
    }


    /// Produce gaussian random distributed algebra
    /// element.
    /// Set default normalisation so that the algebra matrix
    ///    ah = i h = i xi_i \lambda_i
    /// is from distribution
    ///   exp(-Tr h^2 ) = exp(- xi_i^2 / 2 )
    /// I.E. the coefficients of the generators have
    ///   < xi_i^2 > = 1
    ///
    /// Now the off-diag elements have
    ///   <|od|^2> = 1/2 = <od.re^2> + <od.im^2> = 1/4 + 1/4
    ///   1/4 (<xi_a^2> + <xi_b^2>)
    ///
    /// With this convention the inherited method from Vector
    /// is fine and the code below is not needed

    // Algebra &gaussian_random() out_only {
    //
    //     for (int i=0; i<N_a; i++) {
    //         a.e(i) = hila::gaussrand();
    //     }
    //
    //     return *this;
    // }
};

template <int N, typename T>
SUmatrix<N, T> exp(const Algebra<SUmatrix<N, T>> &a) {
    SUmatrix<N, T> m = a.expand();
    return exp(m);
}


#endif