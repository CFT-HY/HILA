#ifndef WILSON_VECTOR_H
#define WILSON_VECTOR_H

#include "datatypes/cmplx.h"
#include "datatypes/matrix.h"
#include "datatypes/sun.h"
#include "plumbing/coordinates.h"
#include "plumbing/random.h"
#include "datatypes/compoundvector.h"

#if (NDIM == 5 || NDIM == 4)
#define Gammadim 4
#define NGamma 5
enum class GammaMatrix : int { g0 = 0, g1, g2, g3, g5 };

constexpr GammaMatrix gamma0 = GammaMatrix::g0;
constexpr GammaMatrix gamma1 = GammaMatrix::g1;
constexpr GammaMatrix gamma2 = GammaMatrix::g2;
constexpr GammaMatrix gamma3 = GammaMatrix::g3;
constexpr GammaMatrix gamma5 = GammaMatrix::g5;

// List order of directions for conveniance
GammaMatrix gamma_matrix[5] = {gamma1, gamma2, gamma3, gamma0, gamma5};

#elif (NDIM == 3 || NDIM == 2)
#define Gammadim 2
#define NGamma 3
enum class GammaMatrix : unsigned { g0 = 0, g1, g2 };

constexpr GammaMatrix gamma0 = GammaMatrix::g0;
constexpr GammaMatrix gamma1 = GammaMatrix::g1;
constexpr GammaMatrix gamma2 = GammaMatrix::g2;

// List order of directions for conveniance
GammaMatrix gamma_matrix[3] = {gamma0, gamma1, gamma2};

#else
static_assert(false, "Wilson fermions only implemented for 1 < NDIM < 6");
#endif


/// Wilson_vector contains the data for a single pseudofermion with the
/// Wilson action. For each gamma dimension, contains a SU(N) vector.
/// Wilson fermions can be multiplied by the gamma matrices gamma0, gamma1
/// ... gamma5 (or up to gamma 2 for 2 and 3 dimensions).
///
/// A gamma matrix can also be projected to a half_Wilson_vector using the
/// constructor half_Wilson_vector(wilson_vector, dir, sign). This multiplies
/// the vector with 0.5*(1-sign*gamma_dir), leaving half the degrees of freedom.
/// A full Wilson vector can be recovered using half_Wilson_vector::grow.
///
/// use type WilsonVector_t, which can be used to define WilsonVector
/// and HalfWilsonVector -types


template <int Nvectors, int N, typename T>
class WilsonVector_t {

  public:
    using base_type = T;
    using argument_type = Vector<N, Complex<T>>;

    Vector<N, Complex<T>> c[Nvectors];

    WilsonVector_t() = default;
    WilsonVector_t(const WilsonVector_t &m) = default;
    ~WilsonVector_t() = default;

    // construct from 0
    WilsonVector_t(std::nullptr_t n) out_only {
        for (int i = 0; i < Nvectors; i++)
            c[i] = 0;
    }

    // and different type WilsonVector
    template <typename A>
    WilsonVector_t(const WilsonVector_t<Nvectors, N, A> &m) out_only {
        for (int i = 0; i < Nvectors; i++) {
            c[i] = m.c[i];
        }
    }

    // Define constant methods rows(), columns() - may be useful in template code
    static constexpr int rows() const {
        return N * Nvectors;
    }
    static constexpr int columns() const {
        return 1;
    }

    /// unary -
    inline WilsonVector_t operator-() const {
        WilsonVector_t res;
        for (int i = 0; i < Nvectors; i++) {
            res.c[i] = -c[i];
        }
        return res;
    }

    /// unary +
    inline const WilsonVector_t &operator+() const {
        return *this;
    }

    /// assign from 0
    inline WilsonVector_t &operator=(const std::nullptr_t &z) out_only {
        for (int i = 0; i < Nvectors; i++) {
            c[i] = 0;
        }
        return *this;
    }

    /// Assign from WilsonVector
    template <typename S>
    inline WilsonVector_t &operator=(const WilsonVector_t<Nvectors, N, S> &rhs) out_only {
        for (int i = 0; i < Nvectors; i++) {
            c[i] = rhs.c[i];
        }
        return *this;
    }

    /// Add assign from Wvec
    template <typename S>
    inline WilsonVector_t &operator+=(const WilsonVector_t<Nvectors, N, S> &rhs) {
        for (int i = 0; i < Nvectors; i++) {
            c[i] += rhs.c[i];
        }
        return *this;
    }

    /// Sub assign from Wvec
    template <typename S>
    inline WilsonVector_t &operator-=(const WilsonVector_t<Nvectors, N, S> &rhs) {
        for (int i = 0; i < Nvectors; i++) {
            c[i] -= rhs.c[i];
        }
        return *this;
    }

    /// Mul assign by scalar
    template <typename S, std::enable_if_t<hila::is_complex_or_arithmetic<S>::value> = 0>
    inline WilsonVector_t &operator*=(const S rhs) {
        for (int i = 0; i < Nvectors; i++) {
            c[i] *= rhs;
        }
        return *this;
    }

    /// Div assign by scalar
    template <typename S, std::enable_if_t<hila::is_complex_or_arithmetic<S>::value> = 0>
    inline WilsonVector_t &operator/=(const S rhs) {
        for (int i = 0; i < Nvectors; i++) {
            c[i] /= rhs;
        }
        return *this;
    }

    /// fill method
    template <typename S, std::enable_if_t<hila::is_complex_or_arithmetic<S>::value> = 0>
    inline WilsonVector_t &fill(const S rhs) out_only {
        for (int i = 0; i < Nvectors; i++) {
            c[i].fill(S);
        }
        return *this;
    }

    /// cast to Array
    inline Array<N * Nvectors, 1, T> &asArray() const_function {
        return *reinterpret_cast<Array<N * Nvectors, 1, T> *>(this);
    }
    const Array<N * Nvectors, 1, T> &asArray() const {
        return *reinterpret_cast<const Array<N * Nvectors, 1, T> *>(this);
    }

    // do we need dagger?
    // /// dagger() casts WilsonVector to WilsonVector_dagger type
    // inline const WilsonVector_dagger<N,T> &dagger() const {
    //     return *reinterpret_cast<const WilsonVector_dagger<N,T> *>(this);
    // }

    /// complex conjugate
    inline WilsonVector_t conj() const {
        WilsonVector_t res;
        for (int i = 0; i < Nvectors; i++) {
            res.c[i] = c[i].conj();
        }
        return res;
    }

    // /// transpose = conj() + dagger()
    // inline const WilsonVector_dagger<N,T> & transpose() const {
    //     return (this->conj()).dagger();
    // }

    /// gaussian random
    void gaussian_random(T width = 1.0) {
        for (int i = 0; i < Nvectors; i++) {
            c[i].gaussian_random(width);
        }
    }

    /// norms
    inline T squarenorm() const {
        T r = 0;
        for (int i = 0; i < Nvectors; i++) {
            r += c[i].squarenorm();
        }
        return r;
    }

    inline T norm() const {
        return sqrt(this->squarenorm());
    }

    // dot is this.dagger() * (rhs)
    template <typename S>
    inline auto dot(const Wilson_vector_t<Nvectors, N, S> &rhs) const {
        auto r = c[0].dot(rhs.c[0]);
        for (int i = 1; i < Nvectors; i++) {
            r += c[i].dot(rhs.c[i]);
        }
        return r;
    }

    /// Returns a square matrix, cast into the Matrix<N,N,Complex<T>> -type,
    /// which is the sum of the outer products of the colour vectors
    /// in this Wilson vector and the argument
    template <typename S>
    inline auto outer_product(const Wilson_vector_t<Nvectors, N, S> &rhs) const {
        auto r = c[0].outer_product(rhs.c[0]);
        for (int i = 1; i < Nvectors; i++) {
            r += c[i].outer_product(rhs.c[i]);
        }
        return r;
    }

    std::string str() const {
        std::string text = "";
        for (int i = 0; i < Nvectors; i++) {
            text += c[i].str() + "\n";
        }
        text += "\n";
        return text;
    }
};


/// Define WilsonVector and HalfWilsonVector as aliases
template <int N, typename T>
using WilsonVector = WilsonVector_t<Gammadim, N, T>;

template <int N, typename T>
using HalfWilsonVector = WilsonVector_t<Gammadim / 2, N, T>;


/// lhs * Wvec = Wvec
/// Multiplying with a matrix should multiply each element, not the gamma-
/// dimension.  Last template parameter finds only types which can be multiplied
/// (we'll keep Wvec type nevertheless)
template <int Nv, int N, typename T, typename M, typename R = hila::type_mul<M, T>>
inline WilsonVector_t<Nv, N, T> operator*(const T &lhs, WilsonVector_t<Nv, N, T> rhs) {
    for (int i = 0; i < Nv; i++) {
        rhs.c[i] = lhs * rhs.c[i];
    }
    return rhs;
}

/// Mult with a scalar from right
template <int Nv, int N, typename T, typename M,
          std::enable_if_t<hila::is_complex_or_arithmetic<M>::value, int> = 0>
inline WilsonVector_t<Nv, N, T> operator*(const WilsonVector_t<Nv, N, T> lhs, const M rhs) {
    lhs *= rhs;
    return lhs;
}

/// Div by scalar
template <int Nv, int N, typename T, typename M,
          std::enable_if_t<hila::is_complex_or_arithmetic<M>::value, int> = 0>
inline WilsonVector_t<Nv, N, T> operator/(const WilsonVector_t<Nv, N, T> lhs, const M rhs) {
    lhs /= rhs;
    return lhs;
}

/// add 2 wvects
template <int Nv, int N, typename T, typename M, typename R = hila::type_plus<T, M>>
inline WilsonVector_t<Nv, N, R> operator+(const WilsonVector_t<Nv, N, T> &lhs,
                                          const WilsonVector_t<Nv, N, M> &rhs) {
    WilsonVector_t<Nv, N, R> r;
    for (int i = 0; i < Nv; i++) {
        r.c[i] = lhs.c[i] + rhs.c[i];
    }
    return r;
}

/// sub wvects
template <int Nv, int N, typename T, typename M, typename R = hila::type_minus<T, M>>
inline WilsonVector_t<Nv, N, R> operator-(const WilsonVector_t<Nv, N, T> &lhs,
                                          const WilsonVector_t<Nv, N, M> &rhs) {
    WilsonVector_t<Nv, N, R> r;
    for (int i = 0; i < Nv; i++) {
        r.c[i] = lhs.c[i] - rhs.c[i];
    }
    return r;
}


/// Multiplication with gamma matrices

#if (Gammadim == 4)

template <int N, typename T>
WilsonVector<N, T> operator*(const GammaMatrix gamma, const WilsonVector<N, T> &rhs) {
    WilsonVector<N, T> r;
    switch (gamma) {
    case gamma0:
        r.c[0] = rhs.c[2];
        r.c[1] = rhs.c[3];
        r.c[2] = rhs.c[0];
        r.c[3] = rhs.c[1];
        break;
    case gamma1:
        r.c[0] = I * rhs.c[3];
        r.c[1] = I * rhs.c[2];
        r.c[2] = -I * rhs.c[1];
        r.c[3] = -I * rhs.c[0];
        break;
    case gamma2:
        r.c[0] = -rhs.c[3];
        r.c[1] = rhs.c[2];
        r.c[2] = rhs.c[1];
        r.c[3] = -rhs.c[0];
        break;
    case gamma3:
        r.c[0] = I * rhs.c[2];
        r.c[1] = -I * rhs.c[3];
        r.c[2] = -I * rhs.c[0];
        r.c[3] = I * rhs.c[1];
        break;
    case gamma5:
        r.c[0] = rhs.c[0];
        r.c[1] = rhs.c[1];
        r.c[2] = -rhs.c[2];
        r.c[3] = -rhs.c[3];
        break;
    }
    return r;
}

#elif (Gammadim == 2)

template <int N, typename T>
WilsonVector<N, T> operator*(const GammaMatrix gamma, const WilsonVector<N, T> &rhs) {
    WilsonVector<N, T> r;
    switch (gamma) {
    case gamma0:
        r.c[0] = rhs.c[1];
        r.c[1] = rhs.c[0];
        break;
    case gamma1:
        r.c[0] = Complex<T>(0, -1) * rhs.c[1];
        r.c[1] = Complex<T>(0, 1) * rhs.c[0];
        break;
    case gamma2:
        r.c[0] = rhs.c[0];
        r.c[1] = -rhs.c[1];
        break;
    }
    return r;
}

#endif

///////////////////////////////////////////////////////////////////////////////////
/// half_Wilson_vector is a Wilson vector projected by
/// 1 +- gamma_j and contains half the degrees of freedom
///
/// (1 +- gamma_j) is a projection operator. We will apply the projection
/// to a Wilson_vector and store the result in a half_Wilson_vector. This
/// will store all necessary information in half the amount of data
/// and reduce the effort of multiplying with gauge matrices by half.
///
/// The constructor half_Wilson_vector(Wilson_vector<n, radix> w, Direction d)
/// will take the projection automatically. A positive Direction d corresponds
/// to 1+gamma_j and a negative Direction d corresponds to 1-gamma_j.
///
/// A half_Wilson_vector can be expanded to a full Wilson_vector using the
/// method expand(Direction d). Notice that this requires knowing the Direction
/// that was used to project it. The Direction is not stored.
///
/// Note that the eigenvectors below are normalized to sqrt(2), or |v*v| = 2.
/// This is why we don't explicitly multiply by 2 when expanding to full
/// Wilson_vector.
///
/// gamma(XUP) 			eigenvectors	eigenvalue
///  0  0  0  i		( 1, 0, 0,-i)	  +1
///  0  0  i  0		( 0, 1,-i, 0)	  +1
///  0 -i  0  0		( 1, 0, 0,+i)	  -1
/// -i  0  0  0		( 0, 1,+i, 0)	  -1
///
/// gamma(YUP)			eigenvectors	eigenvalue
///  0  0  0 -1		( 1, 0, 0,-1)	  +1
///  0  0  1  0		( 0, 1, 1, 0)	  +1
///  0  1  0  0		( 1, 0, 0, 1)	  -1
/// -1  0  0  0		( 0, 1,-1, 0)	  -1
///
/// gamma(ZUP)			eigenvectors	eigenvalue
///  0  0  i  0		( 1, 0,-i, 0)	  +1
///  0  0  0 -i		( 0, 1, 0,+i)	  +1
/// -i  0  0  0		( 1, 0,+i, 0)	  -1
///  0  i  0  0		( 0, 1, 0,-i)	  -1
///
/// gamma(TUP)			eigenvectors	eigenvalue
///  0  0  1  0		( 1, 0, 1, 0)	  +1
///  0  0  0  1		( 0, 1, 0, 1)	  +1
///  1  0  0  0		( 1, 0,-1, 0)	  -1
///  0  1  0  0		( 0, 1, 0,-1)	  -1
///
/// gamma(FIVE) 			eigenvectors	eigenvalue
///  1  0  0  0    sq2( 1, 0, 0, 0)   +1
///  0  1  0  0    sq2( 0, 1, 0, 0)   +1
///  0  0 -1  0    sq2( 0, 0, 1, 0)   -1
///  0  0  0 -1    sq2( 0, 0, 0, 1)   -1
//////////////////////////////////////////////////////////////////////////////////

template <int N, typename T>
class HalfWilsonVector {
  public:
    using base_type = T;
    using argument_type = Vector<N, Complex<T>>;

    Vector<N, Complex<T>> c[Gammadim / 2];

    HalfWilsonVector() = default;
    HalfWilsonVector(const HalfWilsonVector &v) = default;
    ~HalfWilsonVector() = default;

    // This will take the projection 1 +- gamma_j
#if (Gammadim == 4)
    HalfWilsonVector(WilsonVector<N, T> w, Direction dir, int sign) {
        switch (dir) {
        case e_x:
            if (sign == 1) {
                c[0] = w.c[0] + I * w.c[3];
                c[1] = w.c[1] + I * w.c[2];
            } else {
                c[0] = w.c[0] - I * w.c[3];
                c[1] = w.c[1] - I * w.c[2];
            }
            break;
        case e_y:
            if (sign == 1) {
                c[0] = w.c[0] - w.c[3];
                c[1] = w.c[1] + w.c[2];
            } else {
                c[0] = w.c[0] + w.c[3];
                c[1] = w.c[1] - w.c[2];
            }
            break;
        case e_z:
            if (sign == 1) {
                c[0] = w.c[0] + I * w.c[2];
                c[1] = w.c[1] - I * w.c[3];
            } else {
                c[0] = w.c[0] - I * w.c[2];
                c[1] = w.c[1] + I * w.c[3];
            }
            break;
        case e_t:
            if (sign == 1) {
                c[0] = w.c[0] + w.c[2];
                c[1] = w.c[1] + w.c[3];
            } else {
                c[0] = w.c[0] - w.c[2];
                c[1] = w.c[1] - w.c[3];
            }
            break;
#if NDIM == 5
        case 4:
            if (sign == 1) {
                c[0] = sqrt(2.0) * w.c[0];
                c[1] = sqrt(2.0) * w.c[1];
            } else {
                c[0] = sqrt(2.0) * w.c[2];
                c[1] = sqrt(2.0) * w.c[3];
            }
            break;
#endif
        default:
            assert(false && "ERROR: Half Wilson vector projection called incorrectly \n");
        }
    }

    WilsonVector<N, T> expand(Direction dir, int sign) const {
        WilsonVector<N, T> r;
        switch (dir) {
        case e_x:
            if (sign == 1) {
                r.c[0] = c[0];
                r.c[1] = c[1];
                r.c[2] = -I * c[1];
                r.c[3] = -I * c[0];
            } else {
                r.c[0] = c[0];
                r.c[1] = c[1];
                r.c[2] = I * c[1];
                r.c[3] = I * c[0];
            }
            break;
        case e_y:
            if (sign == 1) {
                r.c[0] = c[0];
                r.c[1] = c[1];
                r.c[2] = c[1];
                r.c[3] = -c[0];
            } else {
                r.c[0] = c[0];
                r.c[1] = c[1];
                r.c[2] = -c[1];
                r.c[3] = c[0];
            }
            break;
        case e_z:
            if (sign == 1) {
                r.c[0] = c[0];
                r.c[1] = c[1];
                r.c[2] = -I * c[0];
                r.c[3] = I * c[1];
            } else {
                r.c[0] = c[0];
                r.c[1] = c[1];
                r.c[2] = I * c[0];
                r.c[3] = -I * c[1];
            }
            break;
        case e_t:
            if (sign == 1) {
                r.c[0] = c[0];
                r.c[1] = c[1];
                r.c[2] = c[0];
                r.c[3] = c[1];
            } else {
                r.c[0] = c[0];
                r.c[1] = c[1];
                r.c[2] = -c[0];
                r.c[3] = -c[1];
            }
            break;
#if NDIM == 5
        case 4:
            if (sign == 1) {
                r.c[0] = sqrt(2.0) * c[0];
                r.c[1] = sqrt(2.0) * c[1];
                r.c[2] = 0;
                r.c[3] = 0;
            } else {
                r.c[0] = 0;
                r.c[1] = 0;
                r.c[2] = sqrt(2.0) * c[0];
                r.c[3] = sqrt(2.0) * c[1];
            }
            break;
#endif
        default:
            assert(false && "ERROR: Half Wilson vector projection called incorrectly \n");
        }
        return r;
    }

#elif (Gammadim == 2)
    /*
     gamma(e_x) 	 eigenvectors	 eigenvalue
       0  1		      ( 1, 1)	       +1
       1  0		      ( 1,-1)	       -1

     gamma(e_y)		 eigenvectors	 eigenvalue
       0  i	        ( 1, i)	       +1
      -i  0	  	    ( 1,-i)	       -1

     gamma(e_z)		 eigenvectors  eigenvalue
       1  0	        ( 1, 0)	       +1
       0 -1	  	    ( 0, 1)	       -1
    */
    HalfWilsonVector(WilsonVector<N, T> w, Direction dir, int sign) {
        Complex<T> I(0, 1);
        switch (dir) {
        case e_x:
            if (sign == 1) {
                c[0] = w.c[0] + w.c[1];
            } else {
                c[0] = w.c[0] - w.c[1];
            }
            break;
        case e_y:
            if (sign == 1) {
                c[0] = w.c[0] - I * w.c[1];
            } else {
                c[0] = w.c[0] + I * w.c[1];
            }
            break;
#if NDIM == 3
        case e_z:
            if (sign == 1) {
                c[0] = sqrt(2.0) * w.c[0];
            } else {
                c[0] = sqrt(2.0) * w.c[1];
            }
            break;
#endif
        default:
            assert(false && "ERROR: Half Wilson vector projection called incorrectly \n");
        }

        WilsonVector<N, T> expand(Direction dir, int sign) const {
            WilsonVector<N, T> r;
            Complex<T> I(0, 1);
            switch (dir) {
            case e_x:
                if (sign == 1) {
                    r.c[0] = c[0];
                    r.c[1] = c[0];
                } else {
                    r.c[0] = c[0];
                    r.c[1] = -c[0];
                }
                break;
            case e_y:
                if (sign == 1) {
                    r.c[0] = c[0];
                    r.c[1] = I * c[0];
                } else {
                    r.c[0] = c[0];
                    r.c[1] = -I * c[0];
                }
                break;
#if NDIM == 3
            case e_z:
                if (sign == 1) {
                    r.c[0] = sqrt(2.0) * c[0];
                    r.c[1] = 0;
                } else {
                    r.c[0] = 0;
                    r.c[1] = sqrt(2.0) * c[0];
                }
                break;
#endif
            default:
                assert(false && "ERROR: Half Wilson vector projection called incorrectly \n");
            }
            return r;
        }

#endif

    /// Returns the norm squared of (1+-gamma_j) * wilson_vector.
    /// Thus the factor 2.
    inline T squarenorm() {
        T r = 0;
        for (int i = 0; i < Gammadim / 2; i++) {
            r += c[i].squarenorm();
        }
        return r;
    }

    HalfWilsonVector &operator+=(const HalfWilsonVector &rhs) {
        for (int i = 0; i < Gammadim / 2; i++) {
            c[i] += rhs.c[i];
        }
        return *this;
    }

    HalfWilsonVector &operator-=(const HalfWilsonVector &rhs) {
        for (int i = 0; i < Gammadim / 2; i++) {
            c[i] -= rhs.c[i];
        }
        return *this;
    }

    HalfWilsonVector operator-() const {
        HalfWilsonVector r;
        for (int i = 0; i < Gammadim / 2; i++) {
            r.c[i] = -c[i];
        }
        return r;
    }

    std::string str() const {
        std::string text = "";
        for (int i = 0; i < Gammadim / 2; i++) {
            text += c[i].str() + "\n";
        }
        text += "\n";
        return text;
    }
};

template <int N, typename S, typename T>
HalfWilsonVector<N, T> operator*(const S lhs, HalfWilsonVector<N, T> rhs) {
    for (int i = 0; i < Gammadim / 2; i++) {
        rhs.c[i] *= lhs;
    }
    return rhs;
}

template <int N, typename T, typename S>
HalfWilsonVector<N, T> operator*(HalfWilsonVector<N, T> lhs, const S rhs) {
    for (int i = 0; i < Gammadim / 2; i++) {
        lhs.c[i] *= rhs;
    }
    return lhs;
}

template <int N, typename T, typename S>
HalfWilsonVector<N, T> operator/(HalfWilsonVector<N, T> lhs, const S rhs) {
    for (int i = 0; i < Gammadim / 2; i++) {
        lhs.c[i] /= rhs;
    }
    return lhs;
}


template <int N, typename T, typename P, typename R = hila::type_plus<T, P>>
HalfWilsonVector<N, R> operator+(const HalfWilsonVector<N, T> &lhs,
                                 const HalfWilsonVector<N, P> &rhs) {
    HalfWilsonVector<N, R> r;
    for (int i = 0; i < Gammadim / 2; i++) {
        r.c[i] = lhs.c[i] + rhs.c[i];
    }
    return r;
}

template <int N, typename T, typename P, typename R = hila::type_plus<T, P>>
HalfWilsonVector<N, R> operator-(const HalfWilsonVector<N, T> &lhs,
                                 const HalfWilsonVector<N, P> &rhs) {
    HalfWilsonVector<N, R> r;
    for (int i = 0; i < Gammadim / 2; i++) {
        r.c[i] = lhs.c[i] - rhs.c[i];
    }
    return r;
}


#endif