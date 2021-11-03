#ifndef SUN_MATRIX_OPS_H_
#define SUN_MATRIX_OPS_H_

#include "matrix.h"

/// Define type SU<N,type>
template <int n, typename T = double>
using SU = SquareMatrix<n, Complex<T>>;

/// Make the matrix unitary by orthogonalizing the rows
/// There must be a faster way to do this, but this is simple
///  i = 0 ... n-1
///     normalize row i
///     make rows i+1 .. (n-1) orthogonal to row i
template <int n, int m, typename T>
inline void Matrix<n, m, T>::make_unitary() {

    static_assert(n == m, "make_unitary only for square matrices");

    using radix = hila::number_type<T>;

    for (int r = 0; r < n; r++) {

        // normalize row r
        radix n2 = 0;
        // use here function instead of method, works for double/float too
        for (int c = 0; c < n; c++) n2 += ::squarenorm(this->e(r, c));
        n2 = 1.0 / sqrt(n2);
        for (int c = 0; c < n; c++) this->e(r, c) *= n2;

        // Now make rows r+1 .. n-1 orthogonal to row r,
        // by doing j = j - (r^* j) r

        T d;
        for (int j = r + 1; j < n; j++) {
            // dot productof r^* j
            d = 0;
            for (int i = 0; i < n; i++) {
                d += ::conj(this->e(r, i)) * this->e(j, i);
            }
            // and j -= d * r
            for (int i = 0; i < n; i++) {
                this->e(j, i) -= d * this->e(r, i);
            }
        }
    }
}

/// Set the determinant of the SU(N) matrix to 1
template <int n, int m, typename T>
void Matrix<n,m,T>::fix_det() {

    static_assert(n == m && hila::is_complex<T>::value,
                  "fix_det only for complex square matrices");

    T d, factor;
    hila::number_type<T> t;

    d = det(*(this));
    t = d.arg() / n;
    factor = T(cos(-t), sin(-t));
    this->asArray() *= factor;
}

/// Make the matrix special unitary.  Needs to be complex nxn
template <int n, int m, typename T>
void Matrix<n,m,T>::make_SU() {
    static_assert(n == m && hila::is_complex<T>::value,
                  "make_SU only for complex square matrices");

    make_unitary();
    fix_det();
}

///  Calculate exp of a square matrix
///  Go to  order ORDER in the exponential of the matrices
///  matrices, since unitarity seems important.			
///  Evaluation is done as:					
/// 	exp(H) = 1 + H + H^2/2! + H^3/3! ..-
///	           = 1 + H*( 1 + (H/2)*( 1 + (H/3)*( ... )))
///  Do it backwards in order to reduce accumulation of errors

template <int n, int m, typename T>
inline Matrix<n,m,T> exp(const Matrix<n,m,T> & mat, const int order = 30) {
    static_assert(n == m,
                  "exp() only for square matrices");

    Matrix<n,m,T> r;
    constexpr hila::number_type<T> one = 1.0;

    r = mat*(one/order) + one;

    // doing the loop step-by-step should reduce temporaries
    for (int k=order-1; k>1; k--) {
        r *= mat;
        r *= (one/k);
        r += one;
    }
    r *= mat;
    r += one;

    return r;
}


/// Take traceless antihermitean part
/// More precisely, 0.5*( U - U.dagger() ) - 1/n tr(.)
/// Not the most efficient version, but this routine 
/// should not be critical
template <int n, int m, typename T>
inline Matrix<n,m,T> project_antihermitean(const Matrix<n,m,T> & U) {
    static_assert(n == m,
                  "project_antihermitean() only for square matrices");

    Matrix<n,m,T> r;
    r = 0.5*(U - U.dagger());
    r -= r.trace() / n;

    return r;
}




#endif