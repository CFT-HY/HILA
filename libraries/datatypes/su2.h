#ifndef SU2_H_
#define SU2_H_

#include "matrix.h"

template <typename T>
class Algebra;

/// This implementation represents the matrices as
/// $U = d + a i\sigma_1 + b i\sigma_2 + c i\sigma_3$
/// This is in SU(2) if $a^2 + b^2 + c^2 + d^2 = 1$

template <typename T>
class SU2 {
    static_assert(hila::is_floating_point<T>::value, "SU2 requires a floating point type");

  public: // public on purpose
    T a, b, c, d;

  public:
    using base_type = T;
    using argument_type = T;

    SU2() = default;
    ~SU2() = default;
    SU2(const SU2 &) = default;

    /// construct from 'scalar'
    template <typename B, std::enable_if_t<hila::is_assignable<T &, B>::value, int> = 0>
    SU2(const B z) {
        a = 0;
        b = 0;
        c = 0;
        d = z;
    }
    /// initializer list constructor
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    SU2(std::initializer_list<S> rhs) {
        assert(rhs.size() == 4 && "SU2 initializer list size must be 4");
        auto it = rhs.begin();
        a = *(it++);
        b = *(it++);
        c = *(it++);
        d = *(it);
    }
    /// Normalize det = 1 to make sure it's an element of SU2
    inline const SU2<T> &normalize() {
        T len = sqrt(this->det());
        // assert(len != 0);
        // if (len <= 0)
        //     return *this;
        a /= len;
        b /= len;
        c /= len;
        d /= len;
        return *this;
    }

    /// Normalize det = 1 to make sure it's an element of SU2
    inline const SU2<T> &reunitarize() {
        return this->normalize();
    }
    /// complex conjugate transpose
    inline SU2<T> dagger() const {
        SU2<T> ret;
        ret.a = -a;
        ret.b = -b;
        ret.c = -c;
        ret.d = d;
        return ret;
    }
    /// for SU2 same as .dagger()
    // if matrix is not normalized?
    // SU2<T> inverse() const {
    //     return this->dagger();
    // };

    inline T trace() const {
        return 2.0 * d;
    }

    inline T det() const {
        return a * a + b * b + c * c + d * d;
    }
    static constexpr int size() {
        return 2;
    }
    /// unary -
    inline SU2<T> operator-() const {
        SU2<T> res;
        res.a = -a;
        res.b = -b;
        res.c = -c;
        res.d = -d;
        return res;
    }
    /// unary +
    inline SU2<T> operator+() const {
        return *this;
    }
#pragma hila loop_function
    /// assign from another SU2
    template <typename A, std::enable_if_t<hila::is_assignable<T &, A>::value, int> = 0>
    inline SU2<T> &operator=(const SU2<A> &rhs) out_only {
        a = rhs.a;
        b = rhs.b;
        c = rhs.c;
        d = rhs.d;
        return *this;
    }
#pragma hila loop_function
    /// assign from initializer list
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    inline SU2<T> &operator=(std::initializer_list<S> rhs) out_only {
        assert(rhs.size() == 4 && "SU2 initializer list size must be 4");
        auto it = rhs.begin();
        a = *(it++);
        b = *(it++);
        c = *(it++);
        d = *(it);
    }
#pragma hila loop_function
    /// assign from 'scalar'
    template <typename A, std::enable_if_t<hila::is_assignable<T &, A>::value, int> = 0>
    inline SU2<T> &operator=(A rhs) out_only {
        a = 0;
        b = 0;
        c = 0;
        d = rhs;
        return *this;
    }
#pragma hila loop_function
    /// add assign another SU2
    template <typename A,
              std::enable_if_t<hila::is_assignable<T &, hila::type_plus<T &, A>>::value, int> = 0>
    inline SU2<T> &operator+=(const SU2<A> &rhs) {
        a += rhs.a;
        b += rhs.b;
        c += rhs.c;
        d += rhs.d;
        return *this;
    }
#pragma hila loop_function
    /// subtract assign another SU2
    template <typename A,
              std::enable_if_t<hila::is_assignable<T &, hila::type_minus<T &, A>>::value, int> = 0>
    inline SU2<T> &operator-=(const SU2<A> &rhs) {
        a -= rhs.a;
        b -= rhs.b;
        c -= rhs.c;
        d -= rhs.d;
        return *this;
    }
#pragma hila loop_function
    /// add assign ('scalar' * identity matrix)
    template <typename A,
              std::enable_if_t<hila::is_assignable<T &, hila::type_plus<T &, A>>::value, int> = 0>
    inline SU2<T> &operator+=(const A rhs) {
        d += rhs;
        return *this;
    }
#pragma hila loop_function
    /// subtract assign ('scalar' * identity matrix)
    template <typename A,
              std::enable_if_t<hila::is_assignable<T &, hila::type_minus<T &, A>>::value, int> = 0>
    inline SU2<T> &operator-=(const A rhs) {
        d -= rhs;
        return *this;
    }
#pragma hila loop_function
    /// multiply assign scalar
    template <typename A,
              std::enable_if_t<hila::is_assignable<T &, hila::type_mul<T &, A>>::value, int> = 0>
    inline SU2<T> &operator*=(const A rhs) {
        a *= rhs;
        b *= rhs;
        c *= rhs;
        d *= rhs;
        return *this;
    }
#pragma hila loop_function
    /// divide assign scalar
    template <typename A,
              std::enable_if_t<hila::is_assignable<T &, hila::type_div<T &, A>>::value, int> = 0>
    inline SU2<T> &operator/=(const A rhs) {
        a /= rhs;
        b /= rhs;
        c /= rhs;
        d /= rhs;
        return *this;
    }

    /// make random SU2
    SU2<T> &random() out_only {
        double one, two;
        one = hila::gaussrand2(two);
        a = one;
        b = two;
        one = hila::gaussrand2(two);
        c = one;
        d = two;
        return this->normalize();
    }
    /// project SU2 to generators $\lambda_a = 1/2 \sigma_a$
    inline Algebra<SU2<T>> project_to_algebra() const {
        Algebra<SU2<T>> ret;
        ret.a = a;
        ret.b = b;
        ret.c = c;
        return 2.0 * ret; // factor of 2 from normalization, $\lambda_a = 1/2 \sigma_a$
    }
    /// SU2 matrix exp
    inline SU2<T> exp() const {
        // $exp(U) = e^d*(cos(r) + sin(r)/r *(a i\sigma_1 + b i\sigma_2 + c i\sigma_3))$
        // r = sqrt(a^2+b^2+c^2)
        SU2<T> ret;
        T r = sqrt(a * a + b * b + c * c);
        if (r <= 0) { // TODO: c++20 [[unlikely]] / [[likely]] ?
            ret = 1;
            return ret;
        }
        T ed = exp(d);
        T sr = ed * sin(r) / r;
        ret.d = ed * cos(r);
        ret.a *= sr;
        ret.b *= sr;
        ret.c *= sr;
        return ret;
    }
    /// SU2 matrix log, returns SU2 algebra
    inline Algebra<SU2<T>> log() const {
        //$ exp(U) = A => U = log(A)$
        // (a,b,c) -> (a,b,c)*arcCos(r)/r
        Algebra<SU2<T>> ret;
        T r = sqrt(a * a + b * b + c * c);
        if (r <= 0) { // TODO: c++20 [[unlikely]] / [[likely]] ?
            ret = 0;
            return ret;
        }
        r = acos(d) / sqrt(r);
        ret.a = r * this->a;
        ret.b = r * this->b;
        ret.c = r * this->c;
        return ret;
    }

    SquareMatrix<2, Complex<T>> convert_to_2x2_matrix() const {
        SquareMatrix<2, Complex<T>> res;
        res.e(0, 0) = Complex<T>(this->d, this->c);
        res.e(0, 1) = Complex<T>(this->b, this->a);
        res.e(1, 0) = Complex<T>(-this->b, this->a);
        res.e(1, 1) = Complex<T>(this->d, -this->c);
        return res;
    }
};

/// add two SU2's
template <typename A, typename B, typename R = hila::type_plus<A, B>>
inline SU2<R> operator+(const SU2<A> &rhs, const SU2<B> &lhs) {
    SU2<R> ret;
    ret.a = rhs.a + lhs.a;
    ret.b = rhs.b + lhs.b;
    ret.c = rhs.c + lhs.c;
    ret.d = rhs.d + lhs.d;
    return ret;
}
/// SU2 + ('scalar' * identity matrix)
template <typename A, typename B,
          std::enable_if_t<hila::is_assignable<A &, hila::type_plus<A, B>>::value, int> = 0>
inline SU2<A> operator+(const SU2<A> &rhs, const B lhs) {
    SU2<A> ret = rhs;
    ret.d += lhs;
    return ret;
}
/// ('scalar' * identity matrix) + SU2
template <typename A, typename B,
          std::enable_if_t<hila::is_assignable<B &, hila::type_plus<A, B>>::value, int> = 0>
inline SU2<B> operator+(const A rhs, const SU2<B> &lhs) {
    SU2<B> ret = lhs;
    ret.d += rhs;
    return ret;
}
/// subtract two SU2's
template <typename A, typename B, typename R = hila::type_minus<A, B>>
inline SU2<R> operator-(const SU2<A> &rhs, const SU2<B> &lhs) {
    SU2<R> ret;
    ret.a = rhs.a - lhs.a;
    ret.b = rhs.b - lhs.b;
    ret.c = rhs.c - lhs.c;
    ret.d = rhs.d - lhs.d;
    return ret;
}
/// SU2 - ('scalar' * identity matrix)
template <typename A, typename B,
          std::enable_if_t<hila::is_assignable<A &, hila::type_minus<A, B>>::value, int> = 0>
inline SU2<A> operator-(const SU2<A> &rhs, const B lhs) {
    SU2<A> ret = rhs;
    ret.d -= lhs;
    return ret;
}
/// ('scalar' * identity matrix) - SU2
template <typename A, typename B,
          std::enable_if_t<hila::is_assignable<B &, hila::type_minus<A, B>>::value, int> = 0>
inline SU2<B> operator-(const A rhs, const SU2<B> &lhs) {
    SU2<B> ret = lhs;
    ret.d -= rhs;
    return ret;
}
/// multiply two SU2's
template <typename A, typename B, typename R = hila::type_mul<A, B>>
inline SU2<R> operator*(const SU2<A> &x, const SU2<B> &y) {
    SU2<R> ret;
    ret.a = (x.d * y.a + x.a * y.d - x.b * y.c + x.c * y.b);
    ret.b = (x.d * y.b + x.b * y.d - x.c * y.a + x.a * y.c);
    ret.c = (x.d * y.c + x.c * y.d - x.a * y.b + x.b * y.a);
    ret.d = (x.d * y.d - x.a * y.a - x.b * y.b - x.c * y.c);
    return ret;
}
/// SU2 * scalar
template <typename A, typename B,
          std::enable_if_t<hila::is_assignable<A &, hila::type_mul<A, B>>::value, int> = 0>
inline SU2<A> operator*(const SU2<A> &x, const B y) {
    SU2<A> ret;
    ret.a = x.a * y;
    ret.b = x.b * y;
    ret.c = x.c * y;
    ret.d = x.d * y;
    return ret;
}
/// scalar * SU2
template <typename A, typename B,
          std::enable_if_t<hila::is_assignable<A &, hila::type_mul<A, B>>::value, int> = 0>
inline SU2<A> operator*(const B y, const SU2<A> &x) {
    SU2<A> ret;
    ret.a = x.a * y;
    ret.b = x.b * y;
    ret.c = x.c * y;
    ret.d = x.d * y;
    return ret;
}
/// SU2 / scalar
template <typename A, typename B,
          std::enable_if_t<hila::is_assignable<A &, hila::type_div<A, B>>::value, int> = 0>
inline SU2<A> operator/(const SU2<A> &x, const B y) {
    SU2<A> ret;
    ret.a = x.a / y;
    ret.b = x.b / y;
    ret.c = x.c / y;
    ret.d = x.d / y;
    return ret;
}

template <typename T>
inline T trace(const SU2<T> &U) {
    return U.trace();
}
template <typename T>
inline T det(const SU2<T> &U) {
    return U.det();
}
template <typename T>
inline T squarenorm(const SU2<T> &U) {
    return U.squarenorm();
}
template <typename T>
inline SU2<T> exp(const SU2<T> &U) {
    return U.exp();
}
/// SU2 matrix log, returns SU2 Algebra
template <typename T>
inline Algebra<SU2<T>> log(const SU2<T> &U) {
    return U.log();
}
/// std stream op <<
template <typename T>
inline std::ostream &operator<<(std::ostream &strm, const SU2<T> &U) {
    strm << U.a << " " << U.b << " " << U.c << " " << U.d;
    return strm;
}

/// extract SU2 from NxN complex matrix from elements (i,i), (i,j), (j,i), (j,j)
/// i < j should be valid here!  Return matrix is unnormalized
template <typename T, int N, typename Mtype>
inline SU2<T> project_from_matrix(const Matrix_t<N, N, Complex<T>, Mtype> &m, int i, int j) {
    SU2<T> u;
    u.d = (m.e(i, i).re + m.e(j, j).re) * 0.5;
    u.c = (m.e(i, i).im - m.e(j, j).im) * 0.5;
    u.a = (m.e(i, j).im + m.e(j, i).im) * 0.5;
    u.b = (m.e(i, j).re - m.e(j, i).re) * 0.5;
    return u;
}


// SU2 * vector (vec can be complex or real)
template <typename A, typename B>
inline auto operator*(const SU2<A> &lhs, const Vector<2,B> &rhs) {
    return lhs.convert_to_2x2_matrix() * rhs;
}

// horizontalvector * SU2 (vec can be complex or real)
template <typename A, typename B>
inline auto operator*(const HorizontalVector<2,B> &lhs, const SU2<A> &rhs) {
    return rhs * lhs.convert_to_2x2_matrix();
}


/// This implementation represents algebra as
/// $ a i\sigma_1 + b i\sigma_2 + c i\sigma_3$
template <typename T>
class Algebra<SU2<T>> {
  public: // public on purpose
    T a, b, c;

  public:
    using base_type = T;
    using argument_type = T;

    Algebra() = default;
    ~Algebra() = default;
    Algebra(const Algebra &) = default;
    /// construct from 0
    Algebra(const std::nullptr_t &z) {
        a = 0;
        b = 0;
        c = 0;
    }
    /// unary -
    inline Algebra<SU2<T>> operator-() const {
        Algebra<SU2<T>> res;
        res.a = -a;
        res.b = -b;
        res.c = -c;
        return res;
    }
    /// unary +
    inline Algebra<SU2<T>> operator+() const {
        return *this;
    }
#pragma hila loop_function
    /// assign from another Algebra<SU2>
    template <typename A, std::enable_if_t<hila::is_assignable<T &, A>::value, int> = 0>
    inline Algebra<SU2<T>> &operator=(const Algebra<SU2<A>> &rhs) out_only {
        a = rhs.a;
        b = rhs.b;
        c = rhs.c;
        return *this;
    }
#pragma hila loop_function
    /// assign from initializer list
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    inline Algebra<SU2<T>> &operator=(std::initializer_list<S> rhs) out_only {
        assert(rhs.size() == 3 && "Algebra<SU2> initializer list size must be 3");
        auto it = rhs.begin();
        a = *(it++);
        b = *(it++);
        c = *(it);
    }
#pragma hila loop_function
    /// assign from zero
    inline Algebra<SU2<T>> &operator=(const std::nullptr_t &z) out_only {
        a = 0;
        b = 0;
        c = 0;
        return *this;
    }
#pragma hila loop_function
    /// add assign another Algebra<SU2>
    template <typename A,
              std::enable_if_t<hila::is_assignable<T &, hila::type_plus<T &, A>>::value, int> = 0>
    inline Algebra<SU2<T>> &operator+=(const Algebra<SU2<A>> &rhs) {
        a += rhs.a;
        b += rhs.b;
        c += rhs.c;
        return *this;
    }
#pragma hila loop_function
    /// subtract assign another Algebra<SU2>
    template <typename A,
              std::enable_if_t<hila::is_assignable<T &, hila::type_minus<T &, A>>::value, int> = 0>
    inline Algebra<SU2<T>> &operator-=(const Algebra<SU2<A>> &rhs) {
        a -= rhs.a;
        b -= rhs.b;
        c -= rhs.c;
        return *this;
    }
#pragma hila loop_function
    /// multiply assign scalar
    template <typename A,
              std::enable_if_t<hila::is_assignable<T &, hila::type_mul<T &, A>>::value, int> = 0>
    inline Algebra<SU2<T>> &operator*=(const A rhs) {
        a *= rhs;
        b *= rhs;
        c *= rhs;
        return *this;
    }
#pragma hila loop_function
    /// divide assign scalar
    template <typename A,
              std::enable_if_t<hila::is_assignable<T &, hila::type_div<T &, A>>::value, int> = 0>
    inline Algebra<SU2<T>> &operator/=(const A rhs) {
        a /= rhs;
        b /= rhs;
        c /= rhs;
        return *this;
    }
    /// a^2 + b^2 + c^2
    inline T squarenorm() const {
        return a * a + b * b + c * c;
    }
    /// Expand algebra to SU2
    inline SU2<T> expand() const {
        SU2<T> ret;
        ret.a = a;
        ret.b = b;
        ret.c = c;
        ret.d = 0;
        return ret * 0.5; // factor of 1/2 from normalization, $\lambda_a = 1/2 \sigma_a$
    }
    /// SU2 Algebra $exp( E ) = exp( i 1/2 a_n\sigma_n )$ , returns SU2
    inline SU2<T> exp() const {
        // $U = exp(E) = (cos(r) + sin(r)/r *(a i\sigma_1 + b i\sigma_2 + c i\sigma_3))$
        // r = sqrt(a^2+b^2+c^2)
        Algebra<SU2<T>> tmp =
            (*this) * 0.5; // factor of 1/2 from normalization, $\lambda_a = 1/2 \sigma_a$
        SU2<T> ret;
        T r = sqrt(tmp.squarenorm());
        if (r <= 0) { // TODO: c++20 [[unlikely]] / [[likely]] ?
            ret = 1;
            return ret;
        }
        T sr = sin(r) / r;
        ret.d = cos(r);
        ret.a = sr * tmp.a;
        ret.b = sr * tmp.b;
        ret.c = sr * tmp.c;
        return ret;
    }
    ///
    inline Algebra<SU2<T>> &gaussian_random(double width = 1.0) out_only {
        T one, two;
        one = hila::gaussrand2(two);
        a = width * one;
        b = width * two;
        c = width * hila::gaussrand();
        return *this;
    }
};

/// add two Algebra<SU2>'s
template <typename A, typename B, typename R = hila::type_plus<A, B>>
inline Algebra<SU2<R>> operator+(const Algebra<SU2<A>> &rhs, const Algebra<SU2<B>> &lhs) {
    Algebra<SU2<R>> ret;
    ret.a = rhs.a + lhs.a;
    ret.b = rhs.b + lhs.b;
    ret.c = rhs.c + lhs.c;
    return ret;
}
/// subtract two Algebra<SU2>'s
template <typename A, typename B, typename R = hila::type_minus<A, B>>
inline Algebra<SU2<R>> operator-(const Algebra<SU2<A>> &rhs, const Algebra<SU2<B>> &lhs) {
    Algebra<SU2<R>> ret;
    ret.a = rhs.a - lhs.a;
    ret.b = rhs.b - lhs.b;
    ret.c = rhs.c - lhs.c;
    return ret;
}
/// multiply two Algebra<SU2>'s
template <typename A, typename B, typename R = hila::type_mul<A, B>>
inline Algebra<SU2<R>> operator*(const Algebra<SU2<A>> &x, const Algebra<SU2<B>> &y) {
    Algebra<SU2<R>> ret;
    ret.a = (-x.b * y.c + x.c * y.b);
    ret.b = (-x.c * y.a + x.a * y.c);
    ret.c = (-x.a * y.b + x.b * y.a);
    return ret;
}
/// Algebra<SU2> * scalar
template <typename A, typename B,
          std::enable_if_t<hila::is_assignable<A &, hila::type_mul<A, B>>::value, int> = 0>
inline Algebra<SU2<A>> operator*(const Algebra<SU2<A>> &x, const B y) {
    Algebra<SU2<A>> ret;
    ret.a = x.a * y;
    ret.b = x.b * y;
    ret.c = x.c * y;
    return ret;
}
/// scalar * Algebra<SU2>
template <typename A, typename B,
          std::enable_if_t<hila::is_assignable<A &, hila::type_mul<A, B>>::value, int> = 0>
inline Algebra<SU2<A>> operator*(const B y, const Algebra<SU2<A>> &x) {
    Algebra<SU2<A>> ret;
    ret.a = x.a * y;
    ret.b = x.b * y;
    ret.c = x.c * y;
    return ret;
}
/// Algebra<SU2> / scalar
template <typename A, typename B,
          std::enable_if_t<hila::is_assignable<A &, hila::type_div<A, B>>::value, int> = 0>
inline Algebra<SU2<A>> operator/(const Algebra<SU2<A>> &x, const B y) {
    Algebra<SU2<A>> ret;
    ret.a = x.a / y;
    ret.b = x.b / y;
    ret.c = x.c / y;
    return ret;
}

template <typename T>
inline T squarenorm(const Algebra<SU2<T>> &E) {
    return E.squarenorm();
}

template <typename T>
inline SU2<T> exp(const Algebra<SU2<T>> &E) {
    return E.exp();
}

/// std stream op <<
template <typename T>
inline std::ostream &operator<<(std::ostream &strm, const Algebra<SU2<T>> &E) {
    strm << E.a << " " << E.b << " " << E.c;
    return strm;
}

namespace hila {
template <typename T>
std::string prettyprint(const SU2<T> &A, int prec = 8) {
    std::stringstream strm;
    strm.precision(prec);

    strm << "[ " << A.d << u8" 𝟙 + " << A.a << u8" iσ₁ + " << A.b << u8" iσ₂ + " << A.c << u8" iσ₃ ]";
    return strm.str();
}
} // namespace hila


// SU2 matrix multiplication routines ------------------------
//#define nn_a(x, y) (x.d * y.a + x.a * y.d - x.b * y.c + x.c * y.b)
//#define nn_b(x, y) (x.d * y.b + x.b * y.d - x.c * y.a + x.a * y.c)
//#define nn_c(x, y) (x.d * y.c + x.c * y.d - x.a * y.b + x.b * y.a)
//#define nn_d(x, y) (x.d * y.d - x.a * y.a - x.b * y.b - x.c * y.c)
//
//#define na_a(x, y) (-x.d * y.a + x.a * y.d + x.b * y.c - x.c * y.b)
//#define na_b(x, y) (-x.d * y.b + x.b * y.d + x.c * y.a - x.a * y.c)
//#define na_c(x, y) (-x.d * y.c + x.c * y.d + x.a * y.b - x.b * y.a)
//#define na_d(x, y) (x.d * y.d + x.a * y.a + x.b * y.b + x.c * y.c)
//
//#define an_a(x, y) (x.d * y.a - x.a * y.d + x.b * y.c - x.c * y.b)
//#define an_b(x, y) (x.d * y.b - x.b * y.d + x.c * y.a - x.a * y.c)
//#define an_c(x, y) (x.d * y.c - x.c * y.d + x.a * y.b - x.b * y.a)
//#define an_d(x, y) (x.d * y.d + x.a * y.a + x.b * y.b + x.c * y.c)
//
//#define aa_a(x, y) (-x.d * y.a - x.a * y.d - x.b * y.c + x.c * y.b)
//#define aa_b(x, y) (-x.d * y.b - x.b * y.d - x.c * y.a + x.a * y.c)
//#define aa_c(x, y) (-x.d * y.c - x.c * y.d - x.a * y.b + x.b * y.a)
//#define aa_d(x, y) (x.d * y.d - x.a * y.a - x.b * y.b - x.c * y.c)

#endif // SU2_H_
