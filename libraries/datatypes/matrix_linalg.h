#ifndef MATRIX_LINALG_H
#define MATRIX_LINALG_H

#include "datatypes/matrix.h"

#include <limits>

namespace hila {
namespace linalg {


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
 */

template <typename Dtype>
SquareMatrix<2, Dtype> diagonalize_2x2(const double mpp, const double mqq, const Dtype mpq) {

    double abs_mpq = ::abs(mpq);
    double a = (mqq - mpp) / (2 * abs_mpq);
    double t = 1.0 / (std::abs(a) + std::sqrt(a * a + 1.0));
    if (a < 0.0)
        t = -t;
    double c = 1.0 / std::sqrt(t * t + 1.0);

    Dtype s = mpq * (t * c / abs_mpq);
    SquareMatrix<2, Dtype> P;
    P.e(0, 0) = P.e(1, 1) = 1.0 / std::sqrt(t * t + 1.0);
    P.e(0, 1) = s;
    P.e(1, 0) = -::conj(s);
    return P;
}


} // namespace linalg
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
        double abs_mpq = hila::linalg::find_largest_offdiag(M, p, q);

        // if off-diag elements are tiny return

        if (abs_mpq <= eps * sqrt(::abs(eigenvalues.e(p)) * ::abs(eigenvalues.e(q)))) {
            break;
        }

        // Find diagonalizing matrix
        Matrix<2, 2, Dtype> P =
            hila::linalg::diagonalize_2x2(eigenvalues.e(p), eigenvalues.e(q), M.e(p, q));

        M.mult_by_2x2_left(p, q, P.dagger());
        M.mult_by_2x2_right(p, q, P);

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

        V.mult_by_2x2_right(p, q, P);
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
eigen_result<Mtype> Matrix_t<n, m, T, Mtype>::eigen_hermitean(enum hila::sort sorted) const {

    eigen_result<Mtype> res;
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
        double abs_pq = hila::linalg::find_largest_offdiag(B, p, q);

        // check if we're done - only very small off-diags
        if (abs_pq <= eps * sqrt(::abs(::real(B.e(p, p)) * ::real(B.e(q, q)))))
            break;

        auto P = hila::linalg::diagonalize_2x2(::real(B.e(p, p)), ::real(B.e(q, q)), B.e(p, q));

        // now do p,q rotation
        M.mult_by_2x2_right(p, q, P); // only columns p,q change
        V.mult_by_2x2_right(p, q, P);

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
 * Unpivoted rotation is generally faster than pivoted (esp. on gpu), but a bit less accurate
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

                    auto P = hila::linalg::diagonalize_2x2(Bpp, Bqq, Bpq);

                    // now do p,q rotation
                    M.mult_by_2x2_right(p, q, P); // only columns p,q change
                    V.mult_by_2x2_right(p, q, P);
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
svd_result<Mtype> Matrix_t<n, m, T, Mtype>::svd(enum hila::sort sorted) const {

    svd_result<Mtype> res;
    this->svd(res.U, res.singularvalues, res.V, sorted);
    return res;
}

template <int n, int m, typename T, typename Mtype>
svd_result<Mtype> Matrix_t<n, m, T, Mtype>::svd_pivot(enum hila::sort sorted) const {

    svd_result<Mtype> res;
    this->svd_pivot(res.U, res.singularvalues, res.V, sorted);
    return res;
}


namespace hila {
namespace linalg {

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

} // namespace linalg
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
        hila::scalar_type<T> parity = 1, opposite = -1;
        for (int i = 0; i < n; i++) {
            Matrix<n - 1, m - 1, T> minor = hila::linalg::Minor(*this, 0, i);
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

template <int n, int m, typename T, typename Mtype>
T Matrix_t<n, m, T, Mtype>::det_lu() const {

    static_assert(n == m, "Determinant is defined only for square matrix");

    Mtype a(*this); // copy this matrix to a, modify a in place
    Vector<n, hila::scalar_type<T>> vv;

    hila::scalar_type<T> d = 1, big, tmp;
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



#endif // MATRIX_LINALG_H
