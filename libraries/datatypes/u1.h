#ifndef U1_H_
#define U1_H_

#include <type_traits>
#include <sstream>
#include "plumbing/defs.h"
#include "datatypes/cmplx.h"

////////////////////////////////////////////////////////////////
/// Type U1 as a phase angle - interface is analogous to (square) matrix,
/// in order to facilitate compatible templates
////////////////////////////////////////////////////////////////
template <typename T>
class U1 {
  public:
    static_assert(hila::is_arithmetic<T>::value, "U1 requires arithmetic type");
    // content is just the phase angle
    T phase;

    /// std incantation for field types
    using base_type = hila::arithmetic_type<T>;
    using argument_type = T;

    /// define default constructors to ensure std::is_trivial
    U1() = default;
    ~U1() = default;
    U1(const U1 &v) = default;

    // and make non-explicit constructor from 0
    inline U1(const std::nullptr_t &z) {
        phase = 0;
    }

    /// unary -
    inline U1 operator-() const {
        U1 res;
        res.phase = phase + M_PI;
        return res;
    }

    /// unary +
    inline U1 operator+() const {
        return *this;
    }

    /// multiply assign
    template <
        typename S,
        std::enable_if_t<std::is_convertible<hila::type_plus<T, S>, T>::value, int> = 0>
    U1<T> &operator*=(const U1<S> &rhs) {
        phase += rhs.phase;
        return *this;
    }

    /// U1 -> complex number
    inline Complex<T> complex() const {
        return Complex<T>(cos(phase), sin(phase));
    }

    template <typename S, std::enable_if_t<hila::is_arithmetic<S>::value, int> = 0>
    U1 &set_phase(const S val) out_only {
        phase = val;
        return *this;
    }

    template <typename S, std::enable_if_t<hila::is_arithmetic<S>::value, int> = 0>
    U1 &set_phase(const Complex<S> val) out_only {
        phase = val.arg();
        return *this;
    }

    inline U1 conj() const {
        U1 res;
        res.phase = -phase;
        return res;
    }

    inline U1 dagger() const {
        return conj();
    }

    inline T real() const {
        return cos(phase);
    }

    inline T imag() const {
        return sin(phase);
    }


    /// Generate random elements
    U1 &random() out_only {
        phase = M_PI * (2.0 * hila::random() - 1.0);
        return *this;
    }

    U1 &gaussian_random(double width=1.0) out_only {
        phase = hila::gaussrand() * width;
        return *this;
    }
};

/// conjugate
template <typename T>
inline U1<T> conj(const U1<T> arg) {
    return arg.conj();
}
/// real part
template <typename T>
inline T real(const U1<T> arg) {
    return arg.real();
}
/// imaginary part
template <typename T>
inline T imag(const U1<T> arg) {
    return arg.imag();
}

/// and U1*U1
template <typename T>
inline U1<T> operator*(U1<T> a, const U1<T> b) {
    a += b;
    return a;
}


/// Multiply complex number
template <typename T, typename S>
inline Complex<hila::type_mul<T, S>> operator*(const U1<T> a, const Complex<S> &b) {
    Complex<hila::type_mul<T, S>> r;
    r = a.complex() * b;
    return r;
}

/// Multiply complex number
template <typename T, typename S>
inline Complex<hila::type_mul<T, S>> operator*(const Complex<S> &b, const U1<T> a) {
    Complex<hila::type_mul<T, S>> r;
    r = a.complex() * b;
    return r;
}


/// Stream operator
template <typename T>
std::ostream &operator<<(std::ostream &strm, const U1<T> A) {
    return operator<<(strm, A.phase);
}


namespace hila {
// Cast operators to different number type
// cast_to<double>(a);

template <typename Ntype, typename T,
          std::enable_if_t<hila::is_arithmetic_or_extended<T>::value, int> = 0>
U1<Ntype> cast_to(const U1<T> m) {
    U1<Ntype> res;
    res.phase = m.phase;
    return res;
}

} // namespace hila

#endif