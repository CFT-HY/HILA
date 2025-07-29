#ifndef HILA_SUN_MATRIX_H_
#define HILA_SUN_MATRIX_H_

/**
 * @file sun_matrix.h
 * @brief SU(N) Matrix definitions
 * @details This file contains definitions for SU(N) matrix class and specialized Algebra
 * class for SU(N) matrix.
 */

#include "matrix.h"
#include "su2.h"
#include "tools/floating_point_epsilon.h"

/// Define type SU<n,type>
/// Derives from square Matrix<Complex> type

template <typename G>
class Algebra;

template <typename T>
class sparse_el {
  public:
    using base_type = hila::arithmetic_type<T>;
    using argument_type = T;
    int ind1;
    int ind2;
    T val;
    sparse_el<T>() = default;
    ~sparse_el<T>() = default;
    sparse_el<T>(const sparse_el<T> &a) = default;
    sparse_el<T>(int i1, int i2, const T &v) : ind1(i1), ind2(i2), val(v) {}

    inline sparse_el<T> &operator=(const sparse_el<T> &s) & = default;
    sparse_el<T> &set(const sparse_el<T> &s) {
        ind1 = s.ind1;
        ind2 = s.ind2;
        val = s.val;
        return *this;
    }

    sparse_el<T> &set(int i1, int i2, const T &v) {
        ind1 = i1;
        ind2 = i2;
        val = v;
        return *this;
    }
};

template <int N, typename T>
class Alg_gen {
  public:
    using base_type = T;
    using argument_type = Complex<T>;
    int size;
    sparse_el<Complex<T>> c[N];
    Alg_gen() : size(0){};
    ~Alg_gen() = default;
    inline Alg_gen &operator=(const Alg_gen &s) & = default;
    Alg_gen &push(sparse_el<T> tel) {
        c[size] = tel;
        ++size;
        return *this;
    }
    Alg_gen &push(int i1, int i2, const T &v) {
        c[size].set(i1, i2, Complex<T>(v,0.0));
        ++size;
        return *this;
    }
    Alg_gen &push(int i1, int i2, const Complex<T> &v) {
        c[size].set(i1, i2, v);
        ++size;
        return *this;
    }
    sparse_el<T> pop() {
        sparse_el<T> res(c[size - 1]);
        --size;
        return res;
    }
    Matrix<N, N, Complex<T>> to_matrix() const {
        Matrix<N, N, Complex<T>> res=0;
        for (int i = 0; i < size; ++i) {
            res.e(c[i].ind1, c[i].ind2) = c[i].val;
        }
        return res;
    }
    Alg_gen &empty() {
        size = 0;
        return *this;
    }
};

/**
 * @brief Class for SU(N) matrix
 * @details Class for Special unitary group SU(N).
 *
 * SU class is a special case inherited from Matrix_t class. SU specific or overloaded methods are:
 *
 * - SU::make_unitary
 * - SU::fix_det
 * - SU::reunitarize
 * - SU::random
 * - SU::project_to_algebra
 *
 *
 * @tparam N Dimensionality of SU(N) matrix
 * @tparam T Arithmetic type of Complex<T> number which SU(N) matrix elements consist of
 */
template <int N, typename T>
class SU : public Matrix_t<N, N, Complex<T>, SU<N, T>> {

  public:
    // std incantation for field types
    using base_type = T;
    using argument_type = Complex<T>; // constructed from complex

    // get all constructors from base
    using Matrix_t<N, N, Complex<T>, SU<N, T>>::Matrix_t;
    // same with assignent, but delete rvalue assign
    using Matrix_t<N, N, Complex<T>, SU<N, T>>::operator=;
    SU &operator=(const SU &su) && = delete;


    /**
     * @brief Makes the matrix unitary by orthogonalizing the rows
     * @details This is not the most optimal method, but it is simple:
     *
     * Let `SU<N> H `  be a special unitary matrix and define the indicies \f$ i \in
     * \{0,...,N-1\} \f$. The method is as follows:
     * 1. Normalize  ` H.row(i) `
     * 2. Make rows `H.row(i+1) ` ...` H.row(n-1)` orthogonal with respect
     * to row `H.row(i) `
     *
     *
     * @return const SU&
     */
    const SU &make_unitary() {

        for (int r = 0; r < N; r++) {

            // normalize row r
            T n2 = 0;
            // use here function instead of method, works for double/float too
            for (int i = 0; i < N; i++)
                n2 += ::squarenorm(this->e(r, i));
            n2 = 1.0 / sqrt(n2);
            for (int i = 0; i < N; i++)
                this->e(r, i) *= n2;

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

    /**
     * @brief Fix determinant
     * @details Set the determinant of the SU(N) matrix to 1
     * @return const SU&
     */
    const SU &fix_det() {

        Complex<T> d = det(*(this));
        T t = d.arg() / N;
        d = Complex<T>(cos(-t), sin(-t));
        for (int i = 0; i < N * N; i++) {
            this->c[i] *= d;
        }
        return *this;
    }

    /**
     * @brief Reunitarize SU(N) matrix
     * @details Steps to reunitarize are:
     * 1. Make SU(N) matrix unitary
     * 2. Fix determinant
     *
     * @return const SU&
     */
    const SU &reunitarize() {
        make_unitary();
        fix_det();
        return *this;
    }

    /**
     * @brief Generate random SU(N) matrix
     * @details If N=2 random SU(N) matrix is generated by using Pauli matrix representation.
     *
     * If N > 2 then first we generate a random SU(2) matrix, and multiply the Identity matrix I(N)
     * by this matrix using SU::mult_by_2x2_left. After this we SU::reunitarize the matrix due to
     * numerical errors in the multiplication.
     *
     * @param nhits Number of times I(N) matrix is multiplied to generate random SU(N) matrix, only
     * relevant for N > 2s
     * @return const SU&
     */
    const SU &random(int nhits = 16) out_only {

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
            SU<2, T> m2;

            for (int h = 1; h <= nhits; h++) {
                for (int r = 0; r < N - 1; r++)
                    for (int q = r + 1; q < N; q++) {
                        m2.random();
                        this->mult_by_2x2_left(r, q, m2);
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

    ///
    /// Project matrix to antihermitean and traceless algebra
    /// of the group.
    /// suN generators, normalized as
    ///  Tr(\lambda_n \lambda_m) = -1/2 \delta_nm ,
    ///  or equivalently: Tr(\lambda^{\dagger}_n \lambda_m) = 1/2 \delta_nm
    ///
    /// off-diagonal are just that:
    ///    \lambda^od_nm = i/2 for elements nm and mn
    ///    \lambda^od_nm = 1/2 for nm; -1/2 for mn
    ///
    /// diagonals: su(N) has N-1 diag generators
    /// parametrize these recursively:
    ///  \lambda_1 = diag(i,-i,0,0,..)/sqrt(1)/2
    ///  \lambda_2 = diag(i,i,-2i,0,..)/sqrt(3)/2
    ///  \lambda_3 = diag(i,i,i,-3i,0,..)/sqrt(6)/2
    ///  \lambda_k = diag(i,.. i,-ki,0,..)/sqrt(k(k+1)/2)/2
    ///  ..
    ///  \lambda_N-1 = diag(i,.. i,-(N-1)i)/sqrt(N(N-1)/2)/2
    ///
    /// Define \lambda's so that diagonals come first
    ///
    /// Dividing U = U_ah + U_h, U_ah = 1/2 (U - U^+) = a_k \lambda_k + tr.im I/N
    /// =>  Tr(\lambda_k U_ah) = -1/2 a_k = 1/2 (\lambda_k)_lm (u_ml - u_lm^*)
    /// =>  a_k = - (\lambda_k)_lm (u_ml - u_lm^*)
    ///
    /// Thus, for diags,
    ///     a_k = (u_00 + ..u_(k-1)(k-1) - k*u_kk).im 2/(sqrt(2k(k+1)))
    ///
    /// and off-diags:
    /// symm: a_i = -i (u_kj - u_kj^* + u_jk - u_jk^*)/2 = -i i(u_kj.im + u_jk.im)
    ///           = (u_kj.im + u_jk.im)
    /// antisymm:  a_i = -i (i u_kj - i u_jk^* - i u_jk + i u_kj^*)/2
    ///                = (u_kj.re - u_jk.re)

    Algebra<SU<N, T>> project_to_algebra() const {
        // computes real vector a[] of Lie-algebra decomposition coefficients
        // of A[][]=(*this) , s.t.  a[i] = 2 ReTr( \lambda^{\dagger}_i A )

        Algebra<SU<N, T>> a;

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

    Algebra<SU<N, T>> project_to_algebra_scaled(T scf) const {
        // as above, but rescales the output vector by the factor scf
        Algebra<SU<N, T>> a;

        // diagonal generators
        T sum = this->e(0, 0).im;
        for (int i = 1; i < N; i++) {
            a.e(i - 1) = scf * (sum - i * this->e(i, i).im) / sqrt(0.5 * i * (i + 1));
            sum += this->e(i, i).im;
        }

        // Then off-diag bits
        int k = a.n_diag;
        for (int i = 0; i < N - 1; i++) {
            for (int j = i + 1; j < N; j++) {
                auto od = this->e(i, j) - this->e(j, i).conj();
                a.e(k) = scf * od.re;
                a.e(k + 1) = scf * od.im;
                k += 2;
            }
        }

        return a;
    }

    Algebra<SU<N, T>> project_to_algebra(out_only T &onenorm) const {
        // computes and returns vector a[] of real-valued Lie-algebra
        // projection coefficients and sets in addition onenorm
        // to be the 1-norm of a[]
        Algebra<SU<N, T>> a;

        onenorm = 0;
        // diagonal generators
        T sum = this->e(0, 0).im;
        for (int i = 1; i < N; i++) {
            a.e(i - 1) = (sum - i * this->e(i, i).im) / sqrt(0.5 * i * (i + 1));
            sum += this->e(i, i).im;
            onenorm += abs(a.e(i - 1));
        }

        // Then off-diag bits
        int k = a.n_diag;
        for (int i = 0; i < N - 1; i++) {
            for (int j = i + 1; j < N; j++) {
                auto od = this->e(i, j) - this->e(j, i).conj();
                a.e(k) = od.re;
                a.e(k + 1) = od.im;
                onenorm += abs(a.e(k)) + abs(a.e(k + 1));
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
class Algebra<SU<N, T>> : public Matrix_t<N * N - 1, 1, T, Algebra<SU<N, T>>> {
  public:
    // std incantation for field types
    using base_type = hila::arithmetic_type<T>;
    using argument_type = T;

    // storage for the diagonal and off-diag
    // components of the antihermitean traceless matrix
    static constexpr int n_offdiag = N * (N - 1);
    static constexpr int n_diag = N - 1;
    static constexpr int N_a = N * N - 1;

    /// std constructors and operators derived from vector
    using Matrix_t<N * N - 1, 1, T, Algebra<SU<N, T>>>::Matrix_t;
    using Matrix_t<N * N - 1, 1, T, Algebra<SU<N, T>>>::operator=;
    Algebra &operator=(const Algebra &a) && = delete;

    // suN generators, normalized as
    //  Tr(\lambda_i \lambda_j) = -1/2 \delta_ij ,
    // or equivalently: Tr(\lambda_i \lambda_j) = 1/2 \delta_ij
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

    static void generator_list(Alg_gen<N, T>(out_only &gen_list)[N_a]) {
        // set gen_list to contain N_a=N^2-1 generators of su(N) as sparse matrices
        // (sparse -> list of non-zero elements (ind1,ind2,val))
        for (int i = 1; i < N; ++i) {
            gen_list[i - 1].empty();
            for (int j = 0; j < i; ++j) {
                gen_list[i - 1].push(j, j, Complex<T>(0, 1.0 / sqrt(2.0 * (i * (i + 1)))));
            }
            gen_list[i - 1].push(i, i, Complex<T>(0, -i / sqrt(2.0 * (i * (i + 1)))));
        }

        int k = n_diag;
        for (int i = 0; i < N - 1; i++) {
            for (int j = i + 1; j < N; j++) {
                gen_list[k].empty();
                gen_list[k].push(i, j, Complex<T>(0.5, 0));
                gen_list[k].push(j, i, Complex<T>(-0.5, 0));
                ++k;
                gen_list[k].empty();
                gen_list[k].push(i, j, Complex<T>(0, 0.5));
                gen_list[k].push(j, i, Complex<T>(0, 0.5));
                ++k;
            }
        }
    }

    static void generator_product_list(const Alg_gen<N, T>(&gen_list)[N_a], Alg_gen<N, T>(out_only &gen_prod_list)[N_a][N_a]) {
        // set gen_prod_list[i][j]=\lambda_i \lambda_j as sparse matrices
        // (sparse -> list of non-zero elements (ind1,ind2,val))
        Matrix<N, N, Complex<T>> temp;
        int i1, i2, j1, j2;
        for (i1 = 0; i1 < N * N - 1; ++i1) {
            const auto &gen1 = gen_list[i1];
            for (i2 = 0; i2 < N * N - 1; ++i2) {
                temp = 0;
                const auto &gen2 = gen_list[i2];
                for (auto el1 = gen1.c; el1 != gen1.c + gen1.size; ++el1) {
                    for (auto el2 = gen2.c; el2 != gen2.c + gen2.size; ++el2) {
                        if(el1->ind2==el2->ind1) {
                            temp.e(el1->ind1, el2->ind2) += el1->val*el2->val;
                        }
                    }
                }
                for (j1 = 0; j1 < N; ++j1) {
                    for (j2 = 0; j2 < N; ++j2) {
                        if(temp.e(j1,j2)!=0) {
                            gen_prod_list[i1][i2].push(j1, j2, temp.e(j1, j2));
                        }
                    }
                }
            }
        }
    }

    static void generator_product_list(Alg_gen<N, T>(out_only &gen_prod_list)[N_a][N_a]) {
        // set gen_prod_list[i][j]=2 \lambda_i \lambda_j as sparse matrices
        // (sparse -> list of non-zero elements (ind1,ind2,val))
        Alg_gen<N, T> gen_list[N_a];
        generator_list(gen_list);
        generator_product_list(gen_list, gen_prod_list);
    }

    SU<N, T> expand() const {
        SU<N, T> m;

        Vector<N, T> d;

        d.e(0) = (T)0;
        for (int i = 1; i < N; i++) {
            T r = this->c[i - 1] * sqrt(0.5 / (i * (i + 1)));
            // the contributions from 1's above
            for (int j = 0; j < i; j++)
                d.e(j) += r;
            // and set the negative element - no sum here, first contrib
            d.e(i) = -i * r;
        }

        for (int i = 0; i < N; i++)
            m.e(i, i) = Complex<T>(0, d.e(i));

        int k = n_diag;
        T inv2 = 0.5;
        for (int i = 0; i < N - 1; i++)
            for (int j = i + 1; j < N; j++) {
                Complex<T> v(this->c[k] * inv2, this->c[k + 1] * inv2);
                m.e(i, j) = v;
                m.e(j, i) = -v.conj();
                k += 2;
            }

        return m;
    }

    /// expand algebra to scaled matrix rep - antihermitian
    SU<N, T> expand_scaled(T scf) const {
        SU<N, T> m;

        Vector<N, T> d;

        d.e(0) = (T)0;
        for (int i = 1; i < N; i++) {
            T r = this->c[i - 1] * sqrt(0.5 / (i * (i + 1)));
            // the contributions from 1's above
            for (int j = 0; j < i; j++)
                d.e(j) += r;
            // and set the negative element - no sum here, first contrib
            d.e(i) = -i * r;
        }

        for (int i = 0; i < N; i++)
            m.e(i, i) = Complex<T>(0, scf * d.e(i));

        int k = n_diag;
        T inv2 = 0.5 * scf;
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
    // Wrapper for base class' gaussian_random with appropriate
    // default gaussian width for chosen algebra normalization
    Algebra &gaussian_random(T width = sqrt(2.0)) out_only {
        Matrix_t<N * N - 1, 1, T, Algebra<SU<N, T>>>::gaussian_random(width);
        return *this;
    }

    // dot product of two Algebra vectors
    T dot(const Algebra &rhs) const {
        T r = 0.0;
        for (int i = 0; i < N_a; ++i) {
            r += this->e(i) * rhs.e(i);
        }
        return r * 0.5;
    }
};

template <int N, typename T>
SU<N, T> exp(const Algebra<SU<N, T>> &a) {
    SU<N, T> m = a.expand();
    return exp(m);

    // SU<N,T> m = a.expand() * (-I); // make hermitean
    // SquareMatrix<N,Complex<T>> D;
    // Vector<N,T> ev;
    // m.eigen_jacobi(ev,D);
    // Vector<N,Complex<T>> expv;

    // for (int i=0; i<N; i++) expv[i] = exp(I*ev[i]);
    // for (int i=0; i<N; i++) for (int j=0; j<N; j++) {
    //     m.e(i,j) = D.e(i,0) * expv[0] * D.e(j,0).conj();
    //     for (int k=1; k<N; k++)
    //         m.e(i,j) += D.e(i,k) * expv[k] * D.e(j,k).conj();
    // }
    // return m;
}


// overload of matrix exponential with iterative Cayley-Hamilton (ch) defined in matrix.h.
template <int N, typename T>
SU<N, T> chexp(const Algebra<SU<N, T>> &a) {
    SU<N, T> m = a.expand();
    return chexp(m);
}


// overload of matrix exponential with iterative Cayley-Hamilton using
// "chs" implementation (defined in matrix.h) which needs less temporary
// memory, but is a bit slower.
template <int N, typename T>
SU<N, T> chsexp(const Algebra<SU<N, T>> &a) {
    SU<N, T> m = a.expand();
    return chsexp(m);
}


// logarithm of SU(N) matrix with iterative Cayley-Hamilton
template <int N, typename T>
Algebra<SU<N, T>> log(const SU<N, T> &a) {
    int maxit = 5 * N;
    T fprec = fp<T>::epsilon * 10.0 * Algebra<SU<N, T>>::N_a;
    Matrix_t<N, N, Complex<T>, SU<N, T>> pl[N];

    SU<N, T> tmat = a, tmat2;
    Algebra<SU<N, T>> res = 0, tres;
    T trn, rn;
    int it, i;
    for (it = 0; it < maxit; ++it) {
        tres = tmat.project_to_algebra(trn);
        rn = 0;
        for (i = 0; i < Algebra<SU<N, T>>::N_a; ++i) {
            res.e(i) -= tres.e(i);
            rn += abs(res.e(i));
        }
        if (trn < fprec * (rn + 1.0)) {
            break;
        }
        tmat = res.expand();
        chexp(tmat, tmat2, pl);
        mult(a, tmat2, tmat);
    }

    return -res;
}

template <int N, typename T, typename MT, typename fT = hila::arithmetic_type<T>>
void project_to_algebra_bilinear(const Matrix_t<N, N, T, MT> &w1, const Matrix_t<N, N, T, MT> &w2,
                                 out_only Matrix<N * N - 1, N * N - 1, fT> &omat,
                                 const Alg_gen<N, fT> (&genlist)[N * N - 1]) {
    // computes real matrix outmat[i][j] = 2 * Tr(\lambda_i^{\dagger} * w1 * lambda_j * w2)
    // where the list of algebra generators is provided by genlist[]
    int i1, i2;
    fT temp;
    for (i1 = 0; i1 < N * N - 1; ++i1) {
        const auto &gen1 = genlist[i1];
        for (i2 = 0; i2 < N * N - 1; ++i2) {
            temp = 0;
            const auto &gen2 = genlist[i2];
            for (auto el1 = gen1.c; el1 != gen1.c + gen1.size; ++el1) {
                for (auto el2 = gen2.c; el2 != gen2.c + gen2.size; ++el2) {
                    temp += real(el1->val * w1.e(el1->ind2, el2->ind1) * el2->val *
                                 w2.e(el2->ind2, el1->ind1));
                }
            }
            omat.e(i1, i2) = -2.0 * temp;
        }
    }
}

template <int N, typename T, typename MT, typename fT = hila::arithmetic_type<T>>
void project_to_algebra_bilinear(const Matrix_t<N, N, T, MT> &w1, const Matrix_t<N, N, T, MT> &w2,
                                 out_only Matrix<N * N - 1, N * N - 1, fT> &omat) {
    // computes real matrix outmat[i][j] = 2 * Tr(\lambda_i^{\dagger} * w1 * lambda_j * w2)
    Alg_gen<N, fT> genlist[N * N - 1];
    Algebra<SU<N, fT>>::generator_list(genlist);
    project_to_algebra_bilinear(w1, w2, omat, genlist);
}

template <int N, typename T, typename MT, typename fT = hila::arithmetic_type<T>>
void project_to_algebra_bilinear(const Matrix_t<N, N, T, MT> &w1,
                                 out_only Matrix<N * N - 1, N * N - 1, fT> &omat,
                                 const Alg_gen<N, fT> (&genprodlist)[N * N - 1][N * N - 1]) {
    // computes real matrix outmat[i][j] = 2 * ReTr(\lambda_i^{\dagger} * w1 * lambda_j)
    // where the list of algebra generator products is provided by genprodlist[][]
    int i1, i2;
    fT temp;
    for (i1 = 0; i1 < N * N - 1; ++i1) {
        for (i2 = 0; i2 < N * N - 1; ++i2) {
            const auto &gen = genprodlist[i2][i1];
            temp = 0;
            for (auto el1 = gen.c; el1 != gen.c + gen.size; ++el1) {
                temp += real(el1->val * w1.e(el1->ind2, el1->ind1));
            }
            omat.e(i1, i2) = -2.0 * temp;
        }
    }
}

template <int N, typename T, typename MT, typename fT = hila::arithmetic_type<T>>
void project_to_algebra_bilinear(const Matrix_t<N, N, T, MT> &w1,
                                 out_only Matrix<N * N - 1, N * N - 1, fT> &omat) {
    // computes real matrix outmat[i][j] = 2 * ReTr(\lambda_i^{\dagger} * w1 * lambda_j)
    Alg_gen<N, fT> genprodlist[N * N - 1][N * N - 1];
    Algebra<SU<N, fT>>::generator_product_list(genprodlist);
    project_to_algebra_bilinear(w1, omat, genprodlist);
}

template <const int N, typename T, typename MT>
void project_to_algebra_bilinear(const MT (&w)[N][N],
                                 out_only Matrix<N * N - 1, N * N - 1, T> &omat,
                                 const Alg_gen<N, T> (&genlist)[N * N - 1]) {
    // computes real matrix outmat[i][j] = 2 * Re(Tr(\lambda_i^{\dagger} * w1[k][l]) *
    // \lambda_j[l][k]) where the list of algebra generators is provided by genlist[]
    int i1, i2;
    T temp;
    for (i1 = 0; i1 < N * N - 1; ++i1) {
        const auto &gen1 = genlist[i1];
        for (i2 = 0; i2 < N * N - 1; ++i2) {
            temp = 0;
            const auto &gen2 = genlist[i2];
            for (auto el1 = gen1.c; el1 != gen1.c + gen1.size; ++el1) {
                for (auto el2 = gen2.c; el2 != gen2.c + gen2.size; ++el2) {
                    temp +=
                        real(el1->val * w[el2->ind1][el2->ind2].e(el1->ind2, el1->ind1) * el2->val);
                }
            }
            omat.e(i1, i2) = -2.0 * temp;
        }
    }
}

template <int N, typename T, typename MT>
void project_to_algebra_bilinear(const MT (&w)[N][N],
                                 out_only Matrix<N * N - 1, N * N - 1, T> &omat) {
    // computes real matrix outmat[i][j] = 2 * Re(Tr(\lambda_i^{\dagger} * w1[k][l]) *
    // \lambda_j[l][k])
    Alg_gen<N, T> genlist[N * N - 1];
    Algebra<SU<N, T>>::generator_list(genlist);
    project_to_algebra_bilinear(w, omat, genlist);
}

namespace hila {

///
/// Function hila::random(SU<N,T> & m), equivalent to m.random()
// template <int N, typename T>
// void random(out_only SU<N, T> &m) {
//     m.random();
// }


} // namespace hila
#endif