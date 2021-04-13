#ifndef WILSON_VECTOR_H
#define WILSON_VECTOR_H

#include "datatypes/cmplx.h"
#include "datatypes/matrix.h"
#include "datatypes/sun.h"
#include "plumbing/coordinates.h"
#include "plumbing/random.h"

#if (NDIM == 5 || NDIM == 4)
#define Gammadim 4
#define NGamma 5
enum class gamma_matrix_type : unsigned { g0 = 0, g1, g2, g3, g5 };

constexpr gamma_matrix_type gamma0 = gamma_matrix_type::g0;
constexpr gamma_matrix_type gamma1 = gamma_matrix_type::g1;
constexpr gamma_matrix_type gamma2 = gamma_matrix_type::g2;
constexpr gamma_matrix_type gamma3 = gamma_matrix_type::g3;
constexpr gamma_matrix_type gamma5 = gamma_matrix_type::g5;

// List order of directions for conveniance
gamma_matrix_type gamma_matrix[5] = {gamma1, gamma2, gamma3, gamma0, gamma5};

#elif (NDIM == 3 || NDIM == 2)
#define Gammadim 2
#define NGamma 3
enum class gamma_matrix_type : unsigned { g0 = 0, g1, g2 };

constexpr gamma_matrix_type gamma0 = gamma_matrix_type::g0;
constexpr gamma_matrix_type gamma1 = gamma_matrix_type::g1;
constexpr gamma_matrix_type gamma2 = gamma_matrix_type::g2;

// List order of directions for conveniance
gamma_matrix_type gamma_matrix[5] = {gamma0, gamma1, gamma2};

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
template <int N, typename radix> class Wilson_vector {
  public:

    using base_type = hila::number_type<SU_vector<N, radix>>;
    using argument_type = SU_vector<N,radix>;

    SU_vector<N, radix> c[Gammadim];


    Wilson_vector() = default;

    Wilson_vector(SU_vector<N, radix> m) {
        for (int i = 0; i < Gammadim; i++) {
            c[i] = m;
        }
    }

    template <typename scalart, std::enable_if_t<hila::is_arithmetic<scalart>::value, int> = 0>
    Wilson_vector(const Wilson_vector<N, scalart> m) {
        for (int i = 0; i < Gammadim; i++) {
            c[i] = m.c[i];
        }
    }

    template <typename scalart, std::enable_if_t<hila::is_arithmetic<scalart>::value, int> = 0>
    Wilson_vector(const scalart rhs) {
        for (int i = 0; i < Gammadim; i++) {
            c[i] = rhs;
        }
    }

    void gaussian() {
        for (int i = 0; i < Gammadim; i++) {
            c[i].gaussian();
        }
    }

    template <typename scalart, std::enable_if_t<hila::is_arithmetic<scalart>::value, int> = 0>
    Wilson_vector &operator=(const scalart rhs) {
        for (int i = 0; i < Gammadim; i++) {
            c[i] = rhs;
        }
        return *this;
    }

    Wilson_vector &operator+=(const Wilson_vector &rhs) {
        for (int i = 0; i < Gammadim; i++) {
            c[i] += rhs.c[i];
        }
        return *this;
    }

    Wilson_vector &operator-=(const Wilson_vector &rhs) {
        for (int i = 0; i < Gammadim; i++) {
            c[i] -= rhs.c[i];
        }
        return *this;
    }

    Wilson_vector operator-() const {
        Wilson_vector r;
        for (int i = 0; i < Gammadim; i++) {
            r.c[i] = -c[i];
        }
        return r;
    }

    inline radix squarenorm() {
        radix r = 0;
        for (int i = 0; i < Gammadim; i++) {
            r += c[i].squarenorm();
        }
        return r;
    }

    inline Complex<radix> dot(const Wilson_vector &rhs) const {
        Complex<radix> r = 0;
        for (int i = 0; i < Gammadim; i++) {
            r += c[i].dot(rhs.c[i]);
        }
        return r;
    }

    inline radix rdot(const Wilson_vector &rhs) const {
        radix r = (0.0);
        for (int i = 0; i < Gammadim; i++) {
            r += c[i].rdot(rhs.c[i]);
        }
        return r;
    }

    /// Returns a square matrix, cast into the SU(N) matrix type,
    /// which is the sum of the outer products of the SUN vectors
    /// in this Wilson vector and the argument
    inline auto outer_product(const Wilson_vector rhs) const {
        auto r = c[0].outer_product(rhs.c[0]);
        for (int i = 1; i < Gammadim; i++) {
            r += c[i].outer_product(rhs.c[i]);
        }
        return r;
    }

    std::string str() const {
        std::string text = "";
        for (int i = 0; i < Gammadim; i++) {
            text += c[i].str() + "\n";
        }
        text += "\n";
        return text;
    }
};

/// Replace multiplication with any most types by an element-wise multiplication.
/// Multiplying with an SU(N) matrix should multiply each element, not the gamma-
/// dimension
template <int N, typename radix, typename T>
Wilson_vector<N, radix> operator*(const T lhs, const Wilson_vector<N, radix> rhs) {
    Wilson_vector<N, radix> r;
    for (int i = 0; i < Gammadim; i++) {
        r.c[i] = lhs * rhs.c[i];
    }
    return r;
}

/// Replace multiplication with any most types by an element-wise multiplication.
/// Multiplying with an SU(N) matrix should multiply each element, not the gamma-
/// dimension
template <int N, typename radix, typename T>
Wilson_vector<N, radix> operator*(const Wilson_vector<N, radix> lhs, const T rhs) {
    Wilson_vector<N, radix> r;
    for (int i = 0; i < Gammadim; i++) {
        r.c[i] = lhs.c[i] * rhs;
    }
    return r;
}

template <int N, typename radix>
Wilson_vector<N, radix> operator+(const Wilson_vector<N, radix> lhs,
                                  const Wilson_vector<N, radix> rhs) {
    Wilson_vector<N, radix> r;
    for (int i = 0; i < Gammadim; i++) {
        r.c[i] = lhs.c[i] + rhs.c[i];
    }
    return r;
}

template <int N, typename radix>
Wilson_vector<N, radix> operator-(const Wilson_vector<N, radix> lhs,
                                  const Wilson_vector<N, radix> rhs) {
    Wilson_vector<N, radix> r;
    for (int i = 0; i < Gammadim; i++) {
        r.c[i] = lhs.c[i] - rhs.c[i];
    }
    return r;
}

/// Multiplication with gamma matrices

#if (Gammadim == 4)

template <int N, typename radix>
Wilson_vector<N, radix> operator*(const gamma_matrix_type gamma,
                                  const Wilson_vector<N, radix> rhs) {
    Wilson_vector<N, radix> r;
    switch (gamma) {
    case gamma0:
        r.c[0] = rhs.c[2];
        r.c[1] = rhs.c[3];
        r.c[2] = rhs.c[0];
        r.c[3] = rhs.c[1];
        break;
    case gamma1:
        r.c[0] = Complex<radix>(0, 1) * rhs.c[3];
        r.c[1] = Complex<radix>(0, 1) * rhs.c[2];
        r.c[2] = Complex<radix>(0, -1) * rhs.c[1];
        r.c[3] = Complex<radix>(0, -1) * rhs.c[0];
        break;
    case gamma2:
        r.c[0] = -rhs.c[3];
        r.c[1] = rhs.c[2];
        r.c[2] = rhs.c[1];
        r.c[3] = -rhs.c[0];
        break;
    case gamma3:
        r.c[0] = Complex<radix>(0, 1) * rhs.c[2];
        r.c[1] = Complex<radix>(0, -1) * rhs.c[3];
        r.c[2] = Complex<radix>(0, -1) * rhs.c[0];
        r.c[3] = Complex<radix>(0, 1) * rhs.c[1];
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

template <int N, typename radix>
Wilson_vector<N, radix> operator*(const gamma_matrix_type gamma,
                                  const Wilson_vector<N, radix> rhs) {
    Wilson_vector<N, radix> r;
    switch (gamma) {
    case gamma0:
        r.c[0] = rhs.c[1];
        r.c[1] = rhs.c[0];
        break;
    case gamma1:
        r.c[0] = Complex<radix>(0, -1) * rhs.c[1];
        r.c[1] = Complex<radix>(0, 1) * rhs.c[0];
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

template <int N, typename radix> class half_Wilson_vector {
  public:

    using base_type = hila::number_type<SU_vector<N, radix>>;
    using argument_type = SU_vector<N, radix>;

    SU_vector<N, radix> c[Gammadim / 2];

    half_Wilson_vector() = default;

    // This will take the projection 1 +- gamma_j
#if (Gammadim == 4)
    half_Wilson_vector(Wilson_vector<N, radix> w, Direction dir, int sign) {
        Complex<radix> I(0, 1);
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

    Wilson_vector<N, radix> expand(Direction dir, int sign) const {
        Wilson_vector<N, radix> r;
        Complex<radix> I(0, 1);
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
    half_Wilson_vector(Wilson_vector<N, radix> w, Direction dir, int sign) {
        Complex<radix> I(0, 1);
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

        Wilson_vector<N, radix> expand(Direction dir, int sign) const {
            Wilson_vector<N, radix> r;
            Complex<radix> I(0, 1);
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
                assert(false &&
                       "ERROR: Half Wilson vector projection called incorrectly \n");
            }
            return r;
        }

#endif

    /// Returns the norm squared of (1+-gamma_j) * wilson_vector.
    /// Thus the factor 2.
    inline radix squarenorm() {
        radix r = 0;
        for (int i = 0; i < Gammadim / 2; i++) {
            r += c[i].squarenorm();
        }
        return r;
    }

    half_Wilson_vector &operator+=(const half_Wilson_vector &rhs) {
        for (int i = 0; i < Gammadim / 2; i++) {
            c[i] += rhs.c[i];
        }
        return *this;
    }

    half_Wilson_vector &operator-=(const half_Wilson_vector &rhs) {
        for (int i = 0; i < Gammadim / 2; i++) {
            c[i] -= rhs.c[i];
        }
        return *this;
    }

    half_Wilson_vector operator-() const {
        half_Wilson_vector r;
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

template <int N, typename radix, typename T>
half_Wilson_vector<N, radix> operator*(const T lhs,
                                       const half_Wilson_vector<N, radix> rhs) {
    half_Wilson_vector<N, radix> r;
    for (int i = 0; i < Gammadim / 2; i++) {
        r.c[i] = lhs * rhs.c[i];
    }
    return r;
}

template <int N, typename radix, typename T>
half_Wilson_vector<N, radix> operator*(const half_Wilson_vector<N, radix> lhs,
                                       const T rhs) {
    half_Wilson_vector<N, radix> r;
    for (int i = 0; i < Gammadim / 2; i++) {
        r.c[i] = lhs.c[i] * rhs;
    }
    return r;
}

template <int N, typename radix>
half_Wilson_vector<N, radix> operator+(const half_Wilson_vector<N, radix> lhs,
                                       const half_Wilson_vector<N, radix> rhs) {
    half_Wilson_vector<N, radix> r;
    for (int i = 0; i < Gammadim / 2; i++) {
        r.c[i] = lhs.c[i] + rhs.c[i];
    }
    return r;
}

template <int N, typename radix>
half_Wilson_vector<N, radix> operator-(const half_Wilson_vector<N, radix> lhs,
                                       const half_Wilson_vector<N, radix> rhs) {
    half_Wilson_vector<N, radix> r;
    for (int i = 0; i < Gammadim / 2; i++) {
        r.c[i] = lhs.c[i] - rhs.c[i];
    }
    return r;
}

#endif