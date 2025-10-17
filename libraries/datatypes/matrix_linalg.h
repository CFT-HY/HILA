#ifndef MATRIX_LINALG_H
#define MATRIX_LINALG_H

#include "datatypes/matrix.h"

#include <limits>

namespace hila {


/// @internal Find largest offdiag element of Hermitean matrix M

template <int n, typename Mtype>
inline double find_largest_offdiag(const SquareMatrix<n, Mtype> &M, int &p, int &q) {
    double abs_mpq = -1; // guarantees a result for trivial matrix

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
    return ::sqrt(abs_mpq);
}


/** @internal GivensMatrix is a support class for 2x2 rotations, in SU(2)
 *
 *  Matrix is  | c    s |
 *             | -s*  c |
 * s complex or real, c real, and c^2 + |s|^2 = 1.
 *
 * mult_by_Givens_left/right multiply matrices with Givens matrix left/right.
 * where the "Givens" is taken to be unity except on row/col = p,q
 */

template <typename Dtype>
struct GivensMatrix {
    Dtype s;
    double c;

    GivensMatrix dagger() {
        GivensMatrix res;
        res.c = c;
        res.s = -s;
        return res;
    }

    Vector<2, Dtype> mult_vector_left(Vector<2, Dtype> v) const {
        auto t0 = c * v.e(0) + s * v.e(1);
        v.e(1) = c * v.e(1) - ::conj(s) * v.e(0);
        v.e(0) = t0;
        return v;
    }

    RowVector<2, Dtype> mult_row_right(RowVector<2, Dtype> v) const {
        auto t0 = v.e(0) * c - ::conj(s) * v.e(1);
        v.e(1) = v.e(0) * s + v.e(1) * c;
        v.e(0) = t0;
        return v;
    }

    template <typename Mtype>
    void mult_by_Givens_left(Mtype &M, int p, int q) const {
        Vector<2, Dtype> a;
        for (int i = 0; i < M.columns(); i++) {
            a.e(0) = M.e(p, i);
            a.e(1) = M.e(q, i);

            a = this->mult_vector_left(a);

            M.e(p, i) = a.e(0);
            M.e(q, i) = a.e(1);
        }
    }

    template <typename Mtype>
    void mult_by_Givens_right(Mtype &M, int p, int q) const {
        RowVector<2, Dtype> a;
        for (int i = 0; i < M.rows(); i++) {
            a.e(0) = M.e(i, p);
            a.e(1) = M.e(i, q);

            a = this->mult_row_right(a);

            M.e(i, p) = a.e(0);
            M.e(i, q) = a.e(1);
        }
    }
};


/** @internal Do 2x2 eigenvalue analysis for hermitean matrix
 * return 2x2 unitary/orthogonal matrix which diagonalizes the input
 * input matrix is
 *   | mpp  mpq |
 *   | mqp  mqq |
 * where mpp and mqq are real, and mqp = mpg*
 * Name comes from the fact that this is meant to operate on rows p,q of a larger matrix
 *
 *     | c   s    |  p
 * P = |-s*  c*   |  q  = Ppq
 *
 * with P^+ = P^-1 ->  cc* + ss* = 1, det P = 1
 * Thus, P is SU(2) or O(2) matrix
 *
 *     | c*  -s | | mpp mpq | | c   s |
 * M = | s*   c | | mqp mqq | |-s*  c*|
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
 *
 * Parametrize the result as a Givens matrix c, s.
 */


template <typename Dtype>
GivensMatrix<Dtype> diagonalize_2x2(const double mpp, const double mqq, const Dtype mpq) {

    GivensMatrix<Dtype> res;
    double mpq2 = squarenorm(mpq);

    // leave / 2|mpq| away, avoid divby0
    double a = (mqq - mpp);

    // now t is above t / |mpq|
    double t = 2.0 / (abs(a) + sqrt(a * a + 4 * mpq2));
    if (a < 0.0)
        t = -t;
    res.c = 1.0 / sqrt(mpq2 * t * t + 1.0);
    res.s = mpq * (t * res.c);

    return res;
}


} // namespace hila

/**
 * @brief Calculate eigenvalues and vectors of an hermitean or real symmetric matrix
 *
 * Algorithm uses fully pivoted Jacobi rotations.
 *
 * Two interfaces:
 *
 *    H.eigen_hermitean(E, U, [optional: sort]);
 *    E: output is DiagnoalMatrix containing real eigenvalues
 *    U: nxn unitary matrix, columns are normalized eigenvectors
 *
 *    auto res = H.eigen_hermitean([optional: sort]);
 *    This saves the trouble of defining E and U (see below)
 *
 * Thus, H = U E U^*   and   H U = U E
 *
 * Both interfaces allow optional sorting according to eigenvalues:
 *    hila::sort::unsorted [default]
 *    hila::sort::ascending / descending
 *
 * Example:
 * @code {.cpp}
 *   SquareMatrix<n,Complex<double>> M;
 *   M = ... // make unitary
 *   DiagonalMatrix<n,double> eigenvalues;
 *   SquareMatrix<n,Complex<double>> eigenvectors;
 *
 *   int rotations = M.eigen_hermitean(eigenvalues,eigenvectors,hila::sort::ascending);
 * @endcode
 *
 * @note The matrix has to be hermitean/symmetric
 *
 * @param E   diagonal matrix of real eigenvalues
 * @param U   Unitary nxn matrix of eigenvectors
 * @param sorted  sorting of eigenvalues (default:unsorted)
 * @return int  number of jacobi rotations
 */

template <int n, int m, typename T, typename Mtype>
template <typename Et, typename Mt, typename MT>
int Matrix_t<n, m, T, Mtype>::eigen_hermitean(out_only DiagonalMatrix<n, Et> &E,
                                              out_only Matrix_t<n, n, Mt, MT> &U,
                                              enum hila::sort sorted) const {

    static_assert(!hila::contains_complex<T>::value || hila::contains_complex<Mt>::value,
                  "Eigenvector matrix must be complex with complex original matrix");

    static_assert(n == m, "Eigensystem can be solved only for square matrices");

    using Dtype =
        typename std::conditional<hila::contains_complex<T>::value, Complex<double>, double>::type;

    // std::numeric_limits does not exist on cuda/gpu, so use explicit value
    // constexpr double eps = std::numeric_limits<double>::epsilon();
    constexpr double eps = 2.22e-16;


    int rot;
    SquareMatrix<n, Dtype> M, V;
    DiagonalMatrix<n, double> eigenvalues;

    // Do it in double prec; copy fields
    V = 1;
    M = (*this);

    // don't need the imag. parts of diag (are zero)
    eigenvalues = M.diagonal().real();

    for (rot = 0; rot < 100 * n * n; rot++) {

        /* find the largest off-diag element */
        int p, q;
        double abs_mpq = hila::find_largest_offdiag(M, p, q);

        // if off-diag elements are tiny return

        if (abs_mpq <= eps * sqrt(::abs(eigenvalues.e(p)) * ::abs(eigenvalues.e(q)))) {
            break;
        }

        // Find diagonalizing matrix
        auto P = hila::diagonalize_2x2(eigenvalues.e(p), eigenvalues.e(q), M.e(p, q));

        P.dagger().mult_by_Givens_left(M, p, q);
        P.mult_by_Givens_right(M, p, q);

        eigenvalues.e(p) = ::real(M.e(p, p));
        eigenvalues.e(q) = ::real(M.e(q, q));

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

        P.mult_by_Givens_right(V, p, q);
    }

    // we usually enter here through break
    if (sorted == hila::sort::unsorted) {

        // return values and vectors as is
        E = eigenvalues;
        U = V;

    } else {
        // bubble sort eigenvalues to decreasing order
        Vector<n, int> perm;
        E = eigenvalues.sort(perm, sorted);
        U = V.permute_columns(perm);
    }
    return (rot);
}

/**
 * @brief eigenvalues and -vectors of hermitean/symmetric matrix, alternative interface
 *
 * result = eigen_hermitean(hila::sort [optional]) returns
 * struct eigen_result<Mtype>, with fields
 *   eigenvalues: DiagonalMatrix of eigenvalues
 *   eigenvectors: nxn unitary matrix of eigenvectors
 *
 * Example:
 * @code {.cpp}
 *   Field<SquareMatrix<n,double>> M, sn;
 *   M[ALL] = .. // make symmetric
 *
 *   onsites(ALL) {
 *       auto R = M[X].svd();
 *
 *       // Use EV decomposition to evaluate function: sin(M) = U sin(eigenvals) U^T
 *       sn[X] = R.eigenvectors * sin(R.eigenvalues) * R.eigenvectors.dagger();
 *   }
 * @endcode
 *
 * This interface saves the trouble of defining the eigenvalue and -vector variables.
 *
 */
template <int n, int m, typename T, typename Mtype>
hila::eigen_result<Mtype> Matrix_t<n, m, T, Mtype>::eigen_hermitean(enum hila::sort sorted) const {

    hila::eigen_result<Mtype> res;
    this->eigen_hermitean(res.eigenvalues, res.eigenvectors, sorted);
    return res;
}


/**
 * @brief Singular value decomposition: divide matrix A = U S V*, where U,V unitary and S diagonal
 * matrix of real singular values. Fully pivoted Jacobi rotations
 *
 * Use:
 *   M.svd_pivot(U, S, V, [sorted]);
 *
 * The algorithm works by one-sided Jacobi rotations.
 *
 * - U = V = 1
 * - Compute Gram matrix B = A^* A.  B is Hermitean, with B_ii >= 0.
 * - Loop:
 *   - Find largest B_ij, j > i   - full pivoting
 *        if |B_ij|^2 <= tol * B_ii B_jj we're done
 *   - find J to diagonalize J^* [ B_ii B_ij ] J = diag(D_ii, D_jj)
 *                               [ B_ji B_jj ]
 *   - A <- A J ,  Columns i and j change
 *   - V <- V J
 *   - Recompute B_ik and B_jk, k /= i,j  (B_ii = D_ii, B_jj = D_jj, B_ij = 0)  (n^2)
 * - Endloop
 * Now A = U S, A^* A = S* S = S^2
 * - S_ii = A_ii / |A_i|     (norm of column i)
 * - U = A / S
 *
 * @param U
 */

#pragma hila novector
template <int n, int m, typename T, typename Mtype>
template <typename Et, typename Mt, typename MT>
int Matrix_t<n, m, T, Mtype>::svd_pivot(out_only Matrix_t<n, n, Mt, MT> &_U,
                                        out_only DiagonalMatrix<n, Et> &_S,
                                        out_only Matrix_t<n, n, Mt, MT> &_V,
                                        enum hila::sort sorted) const {

    static_assert(!hila::contains_complex<T>::value || hila::contains_complex<Mt>::value,
                  "SVD: diagonalizing matrix must be complex with complex original matrix");

    static_assert(n == m, "SVD can be solved only for square matrices");

    using Dtype =
        typename std::conditional<hila::contains_complex<T>::value, Complex<double>, double>::type;

    // std::numeric_limits does not exist on cuda/gpu, so use explicit value
    // constexpr double eps = std::numeric_limits<double>::epsilon();
    constexpr double eps = 2.22e-16;


    int rot;
    SquareMatrix<n, Dtype> M, V, B;
    DiagonalMatrix<n, double> S;

    // Do it in double prec; copy fields
    V = 1;
    M = (*this);

    // calc Gram matrix - only upper triangle really needed
    B = (M.dagger() * M);

    for (rot = 0; rot < 100 * sqr(n); rot++) {

        /* find the largest element */
        int p, q;
        double abs_pq = hila::find_largest_offdiag(B, p, q);

        // check if we're done - only very small off-diags
        if (abs_pq <= eps * sqrt(::abs(::real(B.e(p, p)) * ::real(B.e(q, q)))))
            break;

        auto P = hila::diagonalize_2x2(::real(B.e(p, p)), ::real(B.e(q, q)), B.e(p, q));

        // now do p,q rotation
        P.mult_by_Givens_right(M, p, q); // only columns p,q change
        P.mult_by_Givens_right(V, p, q);

        // update also B - could rotate B directly but accuracy suffers. Explicit computation
        // expensive, n^2
        auto col = M.column(p);
        for (int i = 0; i < n; i++) {
            if (i != q) {
                B.e(p, i) = col.dot(M.column(i));
                B.e(i, p) = ::conj(B.e(p, i));
            }
        }
        col = M.column(q);
        for (int i = 0; i < n; i++) {
            if (i != p) {
                B.e(q, i) = col.dot(M.column(i));
                B.e(i, q) = ::conj(B.e(q, i));
            }
        }
        B.e(q, p) = 0;
        B.e(p, q) = 0;
    }

    // Now M = U S. Normalize columns

    for (int i = 0; i < n; i++) {
        auto col = M.column(i);
        S.e(i) = col.norm();
        M.set_column(i, col / S.e(i));
    }

    if (sorted == hila::sort::unsorted) {

        // return values and vectors as is
        _S = S;
        _V = V;
        _U = M;

    } else {
        // bubble sort eigenvalues to decreasing order
        Vector<n, int> perm;
        _S = S.sort(perm, sorted);
        _V = V.permute_columns(perm);
        _U = M.permute_columns(perm);
    }
    return (rot);
}

/**
 * @brief Singular value decomposition: divide matrix A = U S V*, where U,V unitary and S diagonal
 * matrix of real singular values. Unpivoted Jacobi rotations.
 *
 * @details Unpivoted rotation is generally faster than pivoted (esp. on gpu), but a bit less
 * accurate
 *
 * Use:
 *   M.svd(U, S, V, [sorted]);
 *
 * Result satisfies M = U S V.dagger()
 *
 * The algorithm works by one-sided Jacobi rotations.
 *
 * - U = V = 1
 * - Loop
 *   - for i=0..(n-2); j=(i+1)..(n-1)
 *     - calculate B_ii, B_ij, B_jj, (B_ij = (A^* A)_ij)
 *     - If B_ij > 0
 *        - diagonalize J^* [ B_ii B_ij ] J = diag(D_ii, D_jj)
 *                          [ B_ji B_jj ]
 *        - A <- A J ,  Columns i and j change
 *        - V <- V J
 *   - endfor
 *   - if nothing changed, break Loop
 * - Endloop
 * Now A = U S, A^* A = S* S = S^2
 * - S_ii = A_ii / |A_i|     (norm of column i)
 * - U = A / S
 */


#pragma hila novector
template <int n, int m, typename T, typename Mtype>
template <typename Et, typename Mt, typename MT>
int Matrix_t<n, m, T, Mtype>::svd(out_only Matrix_t<n, n, Mt, MT> &_U,
                                  out_only DiagonalMatrix<n, Et> &_S,
                                  out_only Matrix_t<n, n, Mt, MT> &_V,
                                  enum hila::sort sorted) const {

    static_assert(!hila::contains_complex<T>::value || hila::contains_complex<Mt>::value,
                  "SVD: diagonalizing matrix must be complex with complex original matrix");

    static_assert(n == m, "SVD can be solved only for square matrices");

    using Dtype =
        typename std::conditional<hila::contains_complex<T>::value, Complex<double>, double>::type;

    // std::numeric_limits does not exist on cuda/gpu, so use explicit value
    // constexpr double eps = std::numeric_limits<double>::epsilon();
    constexpr double eps = 2.22e-16;

    int rot;
    SquareMatrix<n, Dtype> M, V;
    DiagonalMatrix<n, double> S;
    Dtype Bpq;
    double Bpp, Bqq;

    // Do it in double prec; copy fields
    V = 1;
    M = (*this);

    bool cont = true;
    for (rot = 0; cont && rot < 100 * sqr(n);) {

        cont = false;
        for (int p = 0; p < n - 1; p++) {
            bool need_pp = true;
            Vector<n, Dtype> colp;
            for (int q = p + 1; q < n; q++) {

                double Bpp;
                if (need_pp) {
                    colp = M.column(p);
                    Bpp = colp.squarenorm();
                    need_pp = false;
                }

                auto colq = M.column(q);
                double Bqq = colq.squarenorm();

                Dtype Bpq = colp.dot(colq);

                rot++;

                // check if need to rotate
                if (::abs(Bpq) > eps * sqrt(Bpp * Bqq)) {
                    cont = true;
                    need_pp = true;

                    auto P = hila::diagonalize_2x2(Bpp, Bqq, Bpq);

                    // now do p,q rotation
                    P.mult_by_Givens_right(M, p, q); // only columns p,q change
                    P.mult_by_Givens_right(V, p, q);
                }
            }
        }
    }

    // Now M = U S. Normalize columns

    for (int i = 0; i < n; i++) {
        auto col = M.column(i);
        S.e(i) = col.norm();
        M.set_column(i, col / S.e(i));
    }

    if (sorted == hila::sort::unsorted) {

        // return values and vectors as is
        _S = S;
        _V = V;
        _U = M;

    } else {
        // bubble sort eigenvalues to decreasing order
        Vector<n, int> perm;
        _S = S.sort(perm, sorted);
        _V = V.permute_columns(perm);
        _U = M.permute_columns(perm);
    }
    return (rot);
}


/**
 * @brief svd and svd_pivot - alternative interface
 *
 * svd(hila::sort [optional]) returns
 * struct svd_result<Mtype>, with fields
 *   U: nxn unitary matrix
 *   singularvalues: DiagonalMatrix of singular values (real, > 0)
 *   V: nxn unitary matrix
 *
 * auto res = M.svd();
 *       now res.U : nxn unitary
 *           res.singularvalues : DiagonalMatrix of singular values
 *           res.V : nxn unitary
 *
 * result satifies  M = res.U res.singularvalues res.V.dagger()
 */
template <int n, int m, typename T, typename Mtype>
hila::svd_result<Mtype> Matrix_t<n, m, T, Mtype>::svd(enum hila::sort sorted) const {

    hila::svd_result<Mtype> res;
    this->svd(res.U, res.singularvalues, res.V, sorted);
    return res;
}

template <int n, int m, typename T, typename Mtype>
hila::svd_result<Mtype> Matrix_t<n, m, T, Mtype>::svd_pivot(enum hila::sort sorted) const {

    hila::svd_result<Mtype> res;
    this->svd_pivot(res.U, res.singularvalues, res.V, sorted);
    return res;
}


/**
 * @brief LU decomposition
 */
template <int n, int m, typename T, typename Mtype>
hila::LU_result<Mtype> Matrix_t<n, m, T, Mtype>::LU_decompose() const {
    // we know matrix is square here, LU_result type demands it

    hila::LU_result<Mtype> res;
    // Work with res.LU matrix
    res.LU = *this;
    // numeric_limits does not work on cuda
    // constexpr double tol = 5 * 2.22e-16;


    for (int i = 0; i < n; i++)
        res.P[i] = i;

    // loop over each column
    for (int i = 0; i < n; i++) {
        double max = 0.0;
        int imax = 0;

        // find max row elem on column i
        for (int k = i; k < n; k++) {
            auto norm = ::squarenorm(res.LU.e(k, i));
            if (max < norm) {
                max = norm;
                imax = k;
            }
        }
        // if (max < tol)
        //     return 0; // failure  this makes no sense here

        if (imax != i) {
            hila::swap(res.P[i], res.P[imax]);
            for (int j = 0; j < n; j++)
                hila::swap(res.LU.e(i, j), res.LU.e(imax, j));
        }

        for (int j = i + 1; j < n; j++) {
            res.LU.e(j, i) /= res.LU.e(i, i);
            for (int k = i + 1; k < n; k++)
                res.LU.e(j, k) -= res.LU.e(j, i) * res.LU.e(i, k);
        }
    }
    return res;
}


template <typename Mat>
template <int n, typename Vt, typename Rt>
Vector<n, Rt> hila::LU_result<Mat>::solve(const Vector<n, Vt> &rhs) const {

    static_assert(n == Mat::rows(),
                  "Vector size in LU_result::solve() must be the same as the matrix size");

    Vector<n, Rt> res;
    for (int i = 0; i < n; i++) {
        res[i] = rhs[P[i]];

        for (int k = 0; k < i; k++)
            res[i] -= LU.e(i, k) * res[k];
    }

    for (int i = n - 1; i >= 0; i--) {
        for (int k = i + 1; k < n; k++) {
            res[i] -= LU.e(i, k) * res[k];
        }
        res[i] /= LU.e(i, i);
    }

    return res;
}

template <typename Mat>
Mat hila::LU_result<Mat>::invert() const {
    Mat r;
    constexpr int n = Mat::rows();

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            r.e(i, j) = (P[i] == j) ? 1.0 : 0.0;

            for (int k = 0; k < i; k++)
                r.e(i, j) -= LU.e(i, k) * r.e(k, j);
        }

        for (int i = n - 1; i >= 0; i--) {
            for (int k = i + 1; k < n; k++)
                r.e(i, j) -= LU.e(i, k) * r.e(k, j);

            r.e(i, j) /= LU.e(i, i);
        }
    }
    return r;
}

namespace hila {

// templates needed for naive calculation of determinants
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

} // namespace hila

/**
 * @brief Determinant of matrix, using Laplace method
 * @details Defined only for square matrices
 *
 * For perfomance overloads exist for \f$ 1 \times 1 \f$, \f$ 2 \times 2 \f$  and \f$ 3 \times 3 \f$
 * matrices.
 *
 * @param mat matrix to compute determinant for
 * @return T result determinant
 */


#pragma hila novector
template <int n, int m, typename T, typename Mtype>
T Matrix_t<n, m, T, Mtype>::det_laplace() const {

    static_assert(n == m, "determinants defined only for square matrices");

    if constexpr (n == 1) {
        return this->e(0.0);
    }
    if constexpr (n == 2) {
        return this->e(0, 0) * this->e(1, 1) - this->e(1, 0) * this->e(0, 1);
    }
    if constexpr (n == 3) {
        return this->e(0, 0) * (this->e(1, 1) * this->e(2, 2) - this->e(2, 1) * this->e(1, 2)) -
               this->e(0, 1) * (this->e(1, 0) * this->e(2, 2) - this->e(1, 2) * this->e(2, 0)) +
               this->e(0, 2) * (this->e(1, 0) * this->e(2, 1) - this->e(1, 1) * this->e(2, 0));
    }

    if constexpr (n > 3) {
        T result(0);
        hila::arithmetic_type<T> parity = 1, opposite = -1;
        for (int i = 0; i < n; i++) {
            Matrix<n - 1, m - 1, T> minor = hila::Minor(*this, 0, i);
            result += parity * minor.det_laplace() * (*this).e(0, i);
            parity *= opposite;
        }
        return result;
    }
}


/**
 * @brief Matrix determinant with LU decomposition
 * @details Algorithm: numerical Recipes, 2nd ed. p. 47 ff
 * Works for Real and Complex matrices
 * Defined only for square matrices
 *
 * @return Complex<radix> Determinant
 */

#pragma hila novector
template <int n, int m, typename T, typename Mtype>
T Matrix_t<n, m, T, Mtype>::det_lu() const {

    static_assert(n == m, "Determinant is defined only for square matrix");

    Mtype a(*this); // copy this matrix to a, modify a in place
    Vector<n, hila::arithmetic_type<T>> vv;

    hila::arithmetic_type<T> d = 1, big, tmp;
    int imax = -1;
    T ret = 0;

    for (int i = 0; i < n; i++) {
        big = 0;
        for (int j = 0; j < n; j++) {
            if ((tmp = ::squarenorm(a.e(i, j))) > big)
                big = tmp;
        }

        // If full row is 0 det must be 0 too
        if (big == 0)
            return ret;

        // scaling saved
        vv[i] = 1.0 / sqrt(big);
    }

    // loop over rows
    for (int j = 0; j < n; j++) {

        // build the lower diagonal (L)
        for (int i = 0; i < j; i++) {
            for (int k = 0; k < i; k++) {
                a.e(i, j) -= a.e(i, k) * a.e(k, j);
            }
        }

        big = 0;
        // search for the pivot, and start upper triangle
        for (int i = j; i < n; i++) {
            auto csum = a.e(i, j);
            for (int k = 0; k < j; k++) {
                csum -= a.e(i, k) * a.e(k, j);
            }
            a.e(i, j) = csum;
            auto dum = vv[i] * ::abs(csum);
            if (dum >= big) {
                imax = i;
                big = dum;
            }
        }
        // and swap rows if needed
        if (j != imax) {
            for (int k = 0; k < n; k++) {
                hila::swap(a.e(imax, k), a.e(j, k));
            }
            d = -d;
            vv[imax] = vv[j];
        }

        // det must be zero if diag is zero now
        if (::abs(a.e(j, j)) == 0.0)
            return ret; // ret still 0

        // and divide by the pivot
        if (j != n - 1) {
            auto dum = 1 / a.e(j, j);
            for (int i = j + 1; i < n; i++) {
                a.e(i, j) *= dum;
            }
        }
    }

    // form the det
    ret = d;
    for (int j = 0; j < n; j++) {
        ret *= a.e(j, j);
    }

    return (ret);
}


/**
 * @brief determinant function - if matrix size is < 5, use Laplace, otherwise LU
 */

#pragma hila novector
template <int n, int m, typename T, typename Mtype>
T Matrix_t<n, m, T, Mtype>::det() const {
    static_assert(n == m, "Determinant only for square matrix");
    if constexpr (n < 5)
        return this->det_laplace();
    else
        return this->det_lu();
}


/**
 * @brief function (as opposed to method) interfaces to det-functions
 */

template <int n, int m, typename T, typename Mtype>
T det_laplace(const Matrix_t<n, m, T, Mtype> &mat) {
    return mat.det_laplace();
}

template <int n, int m, typename T, typename Mtype>
T det_lu(const Matrix_t<n, m, T, Mtype> &mat) {
    return mat.det_lu();
}

template <int n, int m, typename T, typename Mtype>
T det(const Matrix_t<n, m, T, Mtype> &mat) {
    return mat.det();
}

namespace hila {
/**
 * @brief Inversed diagnal + const. matrix using Sherman-Morrison formula
 */

/**
 * @details Sherman-Morrison formula (generalized to complex) is
 * \f[
 *     (A + u v^{\dagger})^{-1} = A^{-1} - \frac{A^{-1} u v^{\dagger} A^{-1}}{(1 + v^{\dagger}
 A^{-1} u)}, \f]
 * where \f$A\f$ is invertible matrix and \f$u\f$,\f$v\f$ are vectors with outer product
 * \f$u v^{\dagger}\f$.
 * Let's specialize this here for the case where \f$A\f$ is diagonal and
 * \f[
 *     u = v = \sqrt{c} [1, 1, 1, ...]^{T}
 * \f]
 * i.e. the inversed matrix \f$M^{-1} = (A + C)^{-1}\f$, where \f$C = c uv^{\dagger}\f$ is constant
 matrix. The inversed matrix \f$ M^{-1}\f$ exists iff \f$(1 + v^{\dagger} A^{-1} u) \neq 0\f$.
 */
template <int N, typename T, typename C,
          std::enable_if_t<hila::is_complex_or_arithmetic<C>::value, int> = 0>
auto invert_diagonal_plus_constant_matrix(const DiagonalMatrix<N, T> &D, const C c) {

    DiagonalMatrix<N, T> B = 1 / D;
    auto tmul = c / (1 + c * trace(B));

    return B - tmul * B.asVector() * B.asVector().transpose();
}

} // namespace hila

#endif // MATRIX_LINALG_H
