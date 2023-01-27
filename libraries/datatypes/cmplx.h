#ifndef CMPLX_H_
#define CMPLX_H_

// let's not include the std::complex
// #include <complex>
// #include <cmath>

#include "plumbing/defs.h"

/// TEMPORARY location for vector intrinsic analogues -- result obvious

template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline T mul_add(T a, T b, T c) {
    return a * b + c;
}

template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline T mul_sub(T a, T b, T c) {
    return a * b - c;
}

template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline T nmul_add(T a, T b, T c) {
    return c - a * b;
}


///////////////////////////////////////////////////////////////////////////////
/// main cmpx type definition
/// Define complex type as a template. This allows Hilapp to replace the internal
/// type with a vector. The datatype T must be an arithmetic type.
///////////////////////////////////////////////////////////////////////////////

template <typename T>
class Complex {

    static_assert((hila::is_arithmetic<T>::value && !std::is_integral<T>::value),
                  "Complex can be used only with floating point type");
    // This incantation is needed to make Field<Complex<>> vectorized

  public:
    using base_type = hila::number_type<T>;
    using argument_type = T;

    // and the content of the complex number
    T re, im;

    Complex<T>() = default;
    ~Complex<T>() = default;
    Complex<T>(const Complex<T> &a) = default;

    // constructor from single complex --IS THIS NEEDED?
    // template <typename A,
    //           std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0 >
    // constexpr Complex<T>(const Complex<A> a) : re(static_cast<T>(a.re)),
    // im(static_cast<T>(a.im)) {}

    // constructor from single scalar value
    // Remember to mark this explicit, we do not want this to be invoked
    // in automatic conversions (there should be methods)

#pragma hila loop_function
    template <typename S, std::enable_if_t<hila::is_arithmetic<S>::value, int> = 0>
    explicit constexpr Complex<T>(const S val) : re(val), im(0) {}

    // allow construction from 0 (nullptr)
    constexpr Complex<T>(const std::nullptr_t n) : re(0), im(0) {}

    // constructor c(a,b)
    template <typename A, typename B, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0,
              std::enable_if_t<hila::is_arithmetic<B>::value, int> = 0>
#pragma hila novector loop_function
    explicit constexpr Complex<T>(const A &a, const B &b) : re(a), im(b) {
    }

    // make also std accessors real() and imag() - don't return reference, because
    // v.real() would not then work inside loops!
    inline T real() const {
        return re;
    }

    inline T imag() const {
        return im;
    }

    inline T &real() const_function {
        return re;
    }

    inline T &imag() const_function {
        return im;
    }
    // automatic casting from Complex<T> -> Complex<A>
    // TODO: ensure this works if A is vector type!
    template <typename A>
#pragma hila loop_function // TODO
    operator Complex<A>() const {
        return Complex<A>(re, im);
    }

    inline Complex<T> &operator=(const Complex<T> &s) = default;

    // Assignment from Complex<A>
    template <typename A>
    inline Complex<T> &operator=(const Complex<A> &s) {
        re = s.re;
        im = s.im;
        return *this;
    }

#pragma hila loop_function
    template <typename S, std::enable_if_t<hila::is_arithmetic<S>::value, int> = 0>
    inline Complex<T> &operator=(S s) {
        re = s;
        im = 0;
        return *this;
    }

    template <typename S>
    bool operator==(const Complex<S> &rhs) const {
        return (re == rhs.re && im == rhs.im);
    }

    template <typename S>
    bool operator!=(const Complex<S> &rhs) const {
        return (re != rhs.re || im != rhs.im);
    }

    inline T squarenorm() const {
        return re * re + im * im;
    }

    inline T abs() const {
        return sqrt(squarenorm());
    }
    inline T arg() const {
        return atan2(im, re);
    }

    inline Complex<T> conj() const {
        return Complex<T>(re, -im);
    }

    // alias dagger to conjugate
    inline Complex<T> dagger() const {
        return Complex<T>(re, -im);
    }

    inline Complex<T> polar(const T r, const T theta) out_only {
        re = r * cos(theta);
        im = r * sin(theta);
        return *this;
    }

    inline Complex<T> &random() out_only {
        re = hila::random();
        im = hila::random();
        return *this;
    }

    inline Complex<T> &gaussian_random(double width = 1.0) out_only {
        double d;
        re = hila::gaussrand2(d) * width;
        im = d * width;
        return *this;
    }

    // unary + and -
    inline Complex<T> operator+() const {
        return *this;
    }
    inline Complex<T> operator-() const {
        return Complex<T>(-re, -im);
    }

    // mark += and *= as loop_functions, because these can be used
    // in reductions implicitly
#pragma hila loop_function
    template <typename A>
    inline Complex<T> &operator+=(const Complex<A> &lhs) {
        re += lhs.re;
        im += lhs.im;
        return *this;
    }

    template <typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
    inline Complex<T> &operator+=(const A &a) {
        re += a;
        return *this;
    }

    template <typename A>
    inline Complex<T> &operator-=(const Complex<A> &lhs) {
        re -= lhs.re;
        im -= lhs.im;
        return *this;
    }

    template <typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
    inline Complex<T> &operator-=(const A &a) {
        re -= a;
        return *this;
    }

    // inline Complex<T> & operator*= (const Complex<T> & lhs) {
    //   T r = re * lhs.re - im * lhs.im;
    //   im  = im * lhs.re + re * lhs.im;
    //   re = r;
    //   return *this;
    // }

#pragma hila loop_function
    template <typename A>
    inline Complex<T> &operator*=(const Complex<A> &lhs) {
        T r = mul_sub(re, lhs.re, im * lhs.im); // a*b-c
        im = mul_add(im, lhs.re, re * lhs.im);  // a*b+c
        re = r;
        return *this;
    }

    template <typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
    inline Complex<T> &operator*=(const A a) {
        re *= a;
        im *= a;
        return *this;
    }

    // a/b = a b*/|b|^2 = (a.re*b.re + a.im*b.im + i(a.im*b.re - a.re*b.im))/|b|^2
    // inline Complex<T> & operator/= (const Complex<T> & lhs) {
    //   T n = lhs.squarenorm();
    //   T r = (re * lhs.re + im * lhs.im)/n;
    //   im  = (im * lhs.re - re * lhs.im)/n;
    //   re = r;
    //   return *this;
    // }
    template <typename A>
    inline Complex<T> &operator/=(const Complex<A> &lhs) {
        T n = lhs.squarenorm();
        T r = mul_add(re, lhs.re, im * lhs.im) / n; // a*b+c
        im = mul_sub(im, lhs.re, re * lhs.im) / n;  // a*b-c
        re = r;
        return *this;
    }

    template <typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
    inline Complex<T> &operator/=(const A &a) {
        re /= a;
        im /= a;
        return *this;
    }

    // define also increment and decrement operators (though not so intuitive for complex)

    inline Complex<T> &operator++() {
        this->re++;
        return *this;
    }

    inline Complex<T> &operator--() {
        this->re--;
        return *this;
    }

    inline Complex<T> operator++(int) {
        Complex<T> a = *this;
        this->re++;
        return a;
    }

    inline Complex<T> operator--(int) {
        Complex<T> a = *this;
        this->re--;
        return a;
    }

    std::string str(int prec = 8, char separator = ' ') const {
        std::stringstream ss;
        ss.precision(prec);
        ss << re << separator << im;
        return ss.str();
    }


    // Convenience method a.conj_mul(b) == a^* b
    template <typename A>
    inline Complex<T> conj_mul(const Complex<A> &b) const {
        return Complex<T>(re * b.re + im * b.im, re * b.im - im * b.re);
    }

    // Convenience method a.mul_conj(b) == a * b^*
    template <typename A>
    inline Complex<T> mul_conj(const Complex<A> &b) const {
        return Complex<T>(re * b.re + im * b.im, im * b.re - re * b.im);
    }

    // cast to another number type (IS THIS NEEDED FOR COMPLEX?)
    template <typename Ntype>
    Complex<Ntype> cast_to() const {
        Complex<Ntype> res;
        res.re = re;
        res.im = im;
        return res;
    }
};


namespace hila {
////////////////////////////////////////////////////////////////////////
/// Section for adding complex utility templates (in namespace hila)
/// hila::is_complex_or_arithmetic<T>::value
/// hila::contains_complex<T>::value
/// hila::underlying_type<T>   returns Complex or arithmetic type
/// hila::ntype_op<A,B>        returns the conventionally upgraded complex or scalar number type
/// hila::complex_x_scalar_type<A,B>  type of operation Complex<A> * scalar<B>


////////////////////////////////////////////////////////////////////////
/// Define hila::is_complex<T>::value -template, using specialization
template <typename T>
struct is_complex : std::integral_constant<bool, false> {};

template <typename T>
struct is_complex<Complex<T>> : std::integral_constant<bool, true> {};

/// hila::is_complex_or_arithmetic<T>::value
template <typename T>
struct is_complex_or_arithmetic
    : std::integral_constant<bool, hila::is_arithmetic<T>::value || hila::is_complex<T>::value> {};

/////////////////////////////////////////////////////////////////////////
/// Utility to check that the type contains complex numbers
/// Use as contains_complex<T>::value

template <typename T, typename Enable = void>
struct contains_complex : std::integral_constant<bool, false> {};

template <typename T>
struct contains_complex<T, typename std::enable_if_t<hila::is_field_class_type<T>::value>>
    : std::integral_constant<bool, hila::contains_type<T, Complex<hila::number_type<T>>>::value> {};

/////////////////////////////////////////////////////////////////////////
/// Utility hila::underlying_type<T>  returns complex or arithmetic type
/// depending on the type

template <typename T, typename Enable = void, typename N = hila::number_type<T>>
struct complex_or_arithmetic_type_struct {
    using type = N;
};

template <typename T>
struct complex_or_arithmetic_type_struct<
    T, typename std::enable_if_t<hila::contains_complex<T>::value>> {
    using type = Complex<hila::number_type<T>>;
};

template <typename T>
using underlying_type = typename complex_or_arithmetic_type_struct<T>::type;

/////////////////////////////////////////////////////////////////////////
/// Utility hila::ntype_op<T1,T2>  returns real or complex type,
/// "conventionally" defined double/float etc. upgrading.
/// Note:  hila::type_plus<T1,T2> gives e.g. Complex<float>,double -> Complex<float>

template <typename T1, typename T2, typename Enable = void>
struct ntype_op_s {
    using type = hila::type_plus<hila::number_type<T1>, hila::number_type<T2>>;
};

// if one type is complex?
template <typename T1, typename T2>
struct ntype_op_s<T1, T2,
                  typename std::enable_if_t<(hila::contains_complex<T1>::value ||
                                             hila::contains_complex<T2>::value)>> {
    using type = Complex<hila::type_plus<hila::number_type<T1>, hila::number_type<T2>>>;
};

template <typename T1, typename T2>
using ntype_op = typename ntype_op_s<T1, T2>::type;

////////////////////////////////////////////////////////////////////////////
/// utility complex_x_scalar, used in Complex<A> + scalar<B> ops
/// - if result is convertible to Complex<A>, return Complex<A>
/// - otherwise return Complex<type_plus<A,B>>
/// This is done in order to keep vector types as vectors

template <typename A, typename B, typename Enable = void>
struct complex_x_scalar_s {
    using type = Complex<hila::type_plus<A, B>>;
};

template <typename A, typename B>
struct complex_x_scalar_s<
    A, B, typename std::enable_if_t<std::is_assignable<A &, hila::type_plus<A, B>>::value>> {
    using type = Complex<A>;
};

template <typename A, typename B>
using complex_x_scalar_type = typename complex_x_scalar_s<A, B>::type;

////////////////////////////////////////////////////////////////////////
/// as_complex_array(T var)
/// casts the var to a pointer to complex<number_type<T>>
/// assuming var contains complex type.  This enables access
/// of complex elements as
///  as_complex_array(var)[i]
/// comment out as hilapp gets confused at the moment
////////////////////////////////////////////////////////////////////////

// #pragma hila loop_function
// template <typename T, std::enable_if_t<hila::contains_complex<T>::value, int> = 0>
// inline const Complex<hila::number_type<T>> * as_complex_array(const T &var) {
//     return (const Complex<hila::number_type<T>> *)(void *)&var;
// }

// #pragma hila loop_function
// template <typename T, std::enable_if_t<hila::contains_complex<T>::value, int> = 0>
// inline Complex<hila::number_type<T>> * as_complex_array(T &var) {
//     return (Complex<hila::number_type<T>> *)(void *)&var;
// }

////////////////////////////////////////////////////////////////////////
/// get_complex_element(var,i)
///    returns the i:th complex number embedded in variable v
/// set_complex_element(var,i,value)
///    sets the i:th element in var
////////////////////////////////////////////////////////////////////////

template <typename T, std::enable_if_t<hila::contains_complex<T>::value, int> = 0>
inline Complex<hila::number_type<T>> get_complex_in_var(const T &var, int i) {
    return *(reinterpret_cast<const Complex<hila::number_type<T>> *>(&var) + i);
}

template <typename T, std::enable_if_t<hila::contains_complex<T>::value, int> = 0>
inline void set_complex_in_var(T &var, int i, const Complex<hila::number_type<T>> &val) {
    *(reinterpret_cast<Complex<hila::number_type<T>> *>(&var) + i) = val;
}

} // namespace hila

// generic Complex constructor - type from arguments
template <typename T, std::enable_if_t<hila::is_floating_point<T>::value, int> = 0>
Complex<T> complex(const T re, const T im) {
    return Complex<T>(re, im);
}


// functions real(), imag()

template <typename T>
inline T real(const Complex<T> &a) {
    return a.re;
}


template <typename T>
inline T imag(const Complex<T> &a) {
    return a.im;
}

template <typename T>
Complex<T> polar(T r, T arg) {
    Complex<T> res(r * cos(arg), r * sin(arg));
    return res;
}


// template <typename T>
// inline Complex<T> operator+(const Complex<T> & a, const Complex<T> & b) {
//   return Complex<T>(a.re + b.re, a.im + b.im);
// }

template <typename T1, typename T2, typename Tr = hila::type_plus<T1, T2>>
inline Complex<Tr> operator+(const Complex<T1> &a, const Complex<T2> &b) {
    return Complex<Tr>(a.re + b.re, a.im + b.im);
}

// TODO: for avx vector too -- #define new template macro
template <typename T, typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline auto operator+(const Complex<T> &c, const A &a) {
    return hila::complex_x_scalar_type<T, A>(c.re + a, c.im);
}

template <typename T, typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline auto operator+(const A &a, const Complex<T> &c) {
    return hila::complex_x_scalar_type<T, A>(c.re + a, c.im);
}

// -
// template <typename T>
// inline Complex<T> operator-(const Complex<T> & a, const Complex<T> & b) {
//   return Complex<T>(a.re - b.re, a.im - b.im);
// }
template <typename T1, typename T2, typename Tr = hila::type_plus<T1, T2>>
inline Complex<Tr> operator-(const Complex<T1> &a, const Complex<T2> &b) {
    return Complex<Tr>(a.re - b.re, a.im - b.im);
}

// complex - scalar
template <typename T, typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline auto operator-(const Complex<T> &c, const A &a) {
    return hila::complex_x_scalar_type<T, A>(c.re - a, c.im);
}

// scalar - complex
template <typename T, typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline auto operator-(const A &a, const Complex<T> &c) {
    return hila::complex_x_scalar_type<T, A>(a - c.re, -c.im);
}

///
/// complex * complex
template <typename T1, typename T2, typename Tr = hila::type_mul<T1, T2>>
inline Complex<Tr> operator*(const Complex<T1> &a, const Complex<T2> &b) {
    return Complex<Tr>(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}

/// complex * scalar
template <typename T, typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline auto operator*(const Complex<T> &c, const A &a) {
    return hila::complex_x_scalar_type<T, A>(c.re * a, c.im * a);
}

template <typename T, typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline auto operator*(const A &a, const Complex<T> &c) {
    return hila::complex_x_scalar_type<T, A>(a * c.re, a * c.im);
}

// /   a/b = ab*/|b|^2
// template <typename T>
// inline Complex<T> operator/(const Complex<T> & a, const Complex<T> & b) {
//   T n = b.squarenorm();
//   return Complex<T>( (a.re*b.re + a.im*b.im)/n, (a.im*b.re - a.re*b.im)/n );
// }
template <typename T1, typename T2, typename Tr = hila::type_mul<T1, T2>>
inline Complex<Tr> operator/(const Complex<T1> &a, const Complex<T2> &b) {
    T2 n = 1.0 / b.squarenorm();
    return Complex<Tr>((a.re * b.re + a.im * b.im) * n, (a.im * b.re - a.re * b.im) * n);
}

template <typename T, typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline auto operator/(const Complex<T> &c, const A &a) {
    return hila::complex_x_scalar_type<T, A>(c.re / a, c.im / a);
}

// a/c = ac*/|c|^2
template <typename T, typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline auto operator/(const A &a, const Complex<T> &c) {
    T n = c.squarenorm();
    return hila::complex_x_scalar_type<T, A>((a * c.re) / n, -(a * c.im) / n);
}

// write also multiply-add directly with complex numbers
template <typename T>
inline Complex<T> mul_add(const Complex<T> &a, const Complex<T> &b, const Complex<T> &c) {
    // a*b + c
    Complex<T> r;
    T t1 = mul_add(a.re, b.re, c.re);
    T t2 = mul_add(a.re, b.im, c.im);
    r.re = mul_add(a.im, b.im, t1); // -a.im*b.im + a.re*b.re + c.re
    r.im = mul_add(a.im, b.re, t2); // a.im*b.re + a.re*b.im + c.im
    return r;
}

template <typename A, typename B>
inline bool operator==(const Complex<A> &a, const Complex<B> &b) {
    return (a.re == b.re && a.im == b.im);
}

template <typename A, typename B, std::enable_if_t<hila::is_arithmetic<B>::value, int> = 0>
inline bool operator==(const Complex<A> &a, const B b) {
    return (a.re == b && a.im == 0);
}

template <typename A, typename B, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
inline bool operator==(const A a, const Complex<B> &b) {
    return b == a;
}


// Cast operators to different number or Complex type
// cast_to<double>(a);
// cast_to<float>(b);
// Cast from number->number, number->Complex, Complex->Complex OK,
//     Complex->number not.

template <typename Ntype, typename T, std::enable_if_t<hila::is_arithmetic<Ntype>::value, int> = 0>
Complex<Ntype> cast_to(const Complex<T> &m) {
    Complex<Ntype> res;
    res = m;
    return res;
}

//////////////////////////////////////////////////////////////////////////////////
// Some operations in function form.  Useful in templates when the arg type is not known

/// abs
template <typename T>
inline T abs(const Complex<T> &a) {
    return a.abs();
}

/// arg
template <typename T>
inline T arg(const Complex<T> &a) {
    return a.arg();
}

/// Conjugate
template <typename T>
inline Complex<T> conj(const Complex<T> &val) {
    return val.conj();
}

/// Alias dagger to conjugate
template <typename T>
inline Complex<T> dagger(const Complex<T> &val) {
    return val.conj();
}

/// norm_squared
template <typename T>
inline auto squarenorm(const Complex<T> &val) {
    return val.squarenorm();
}

/// random() : set argument to random vals [0,1]
template <typename T>
inline void random(out_only Complex<T> &c) {
    ::random(c.re);
    ::random(c.im);
}

template <typename T>
inline void gaussian_random(out_only Complex<T> &c, double width = 1.0) {
    c.gaussian_random(width);
}

//////////////////////////////////////////////////////////////////////////////////
/// Print a complex value as (re,im)
//////////////////////////////////////////////////////////////////////////////////

template <typename T>
std::ostream &operator<<(std::ostream &strm, const Complex<T> &A) {
    return strm << A.real() << ' ' << A.imag();
}

/// //////////////////////////////////////////////////////////////////////////////
/// Function hila::to_string

namespace hila {
template <typename T>
std::string to_string(const Complex<T> &A, int prec = 8, char separator = ' ') {
    return A.str(prec, separator);
}

template <typename T>
std::string prettyprint(const Complex<T> &A, int prec = 8) {
    std::stringstream ss;
    ss.precision(prec);
    ss << "( " << A.real() << ", " << A.imag() << " )";
    return ss.str();
}

} // namespace hila


//////////////////////////////////////////////////////////////////////////////////
/// Operators to implement imaginary unit 1_i, enablig expressions  3 + 2_i  etc.
/// This is defined as an user-defined literal, which requires an underscore.
////////////////////////////////////////////////////////////////////////////////

constexpr Complex<double> operator""_i(long double a) {
    return Complex<double>{0.0, a};
}

constexpr Complex<double> operator""_i(unsigned long long a) {
    return Complex<double>(0.0, static_cast<double>(a));
}


/////////////////////////////////////////////////////////////////////////////
/// Define Imaginary type, to represent purely imaginary numbers
/// Useful for short-circuiting im * complex -type operations
/// Derived from complex class, so all complex operations should be valid
/////////////////////////////////////////////////////////////////////////////

template <typename T>
class Imaginary_t : public Complex<T> {
  public:
    constexpr Imaginary_t() = default;
    ~Imaginary_t() = default;
    constexpr Imaginary_t(const Imaginary_t &i) = default;

    // construct from scalar
#pragma hila loop_function
    template <typename A, std::enable_if_t<hila::is_arithmetic<A>::value, int> = 0>
    explicit constexpr Imaginary_t(const A v) : Complex<T>(0, v) {}

    constexpr Imaginary_t operator-() const {
        return Imaginary_t<T>(-this->im);
    }
};


template <typename A, typename T, std::enable_if_t<hila::contains_complex<T>::value, int> = 0>
inline auto operator*(const Imaginary_t<A> &iv, const T &c) {
    Complex<hila::number_type<T>> ca, cb;
    T res;
    constexpr int n_cmplx = sizeof(T) / sizeof(Complex<hila::number_type<T>>);
    for (int i = 0; i < n_cmplx; i++) {
        ca = hila::get_complex_in_var(c, i);
        cb.re = -ca.im * iv.imag();
        cb.im = ca.re * iv.imag();
        hila::set_complex_in_var(res, i, cb);
    }
    return res;
}

template <typename A, typename T, std::enable_if_t<hila::contains_complex<T>::value, int> = 0>
inline auto operator*(const T &c, const Imaginary_t<A> &iv) {
    return iv * c;
}

/// Imag * scalar = imag -- note: using std::is_arithmetic, not done for vector types
template <typename A, typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
inline Imaginary_t<A> operator*(Imaginary_t<A> iv, const T &c) {
    iv.imag() *= c;
    return iv;
}

template <typename A, typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
inline Imaginary_t<A> operator*(const T &c, Imaginary_t<A> iv) {
    iv.imag() *= c;
    return iv;
}


////////////////////////////////////////////////////////////////////////////
/// Define imaginary unit class, for imaginary I
/// This is a derived special class of Complex<float>.  Define separately
/// I*cmplx and cmplx*I
/// Real micro-optimization, could have used Imaginary_t(1)
///

class Imaginaryunit_t : public Imaginary_t<double> {
  public:
    constexpr Imaginaryunit_t() : Imaginary_t(1.0) {}
    ~Imaginaryunit_t() = default;
    Imaginaryunit_t(const Imaginaryunit_t &i) = default;

    constexpr Imaginary_t<double> operator-() const {
        return Imaginary_t(-1);
    }

    template <typename T>
    operator Imaginary_t<T>() const {
        return Imaginary_t<T>(1);
    }
};

/// Define constexpr imaginary unit I
#if defined(CUDA) || defined(HIP)
__device__
#endif
    constexpr Imaginaryunit_t I =
        Imaginaryunit_t(); // this fails on GPUs without additional support
// #define I Imaginaryunit_t()

template <typename T, std::enable_if_t<hila::contains_complex<T>::value, int> = 0>
inline auto operator*(const Imaginaryunit_t &iv, T c) {
    Complex<hila::number_type<T>> ca, cb;
    for (int i = 0; i < sizeof(T) / sizeof(Complex<hila::number_type<T>>); i++) {
        ca = hila::get_complex_in_var(c, i);
        cb.re = -ca.im;
        cb.im = ca.re;
        hila::set_complex_in_var(c, i, cb);
    }
    return c;
}

template <typename T, std::enable_if_t<hila::contains_complex<T>::value, int> = 0>
inline auto operator*(T c, const Imaginaryunit_t &iv) {
    return iv * c;
}

template <typename T, std::enable_if_t<hila::is_floating_point<T>::value, int> = 0>
inline auto operator*(const Imaginaryunit_t &iv, T c) {
    return Imaginary_t<T>(c);
}

template <typename T, std::enable_if_t<hila::is_floating_point<T>::value, int> = 0>
inline auto operator*(T c, const Imaginaryunit_t &iv) {
    return Imaginary_t<T>(c);
}

// convert ints to double
template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0>
inline auto operator*(T c, const Imaginaryunit_t &iv) {
    return Imaginary_t<double>(c);
}

template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0>
inline auto operator*(const Imaginaryunit_t &iv, T c) {
    return c * iv;
}


///////////////////////////////////////////////////////////////////////////////
/// Set of complex functions
///////////////////////////////////////////////////////////////////////////////


/// exp(z)
template <typename T>
inline Complex<T> exp(const Complex<T> z) {
    return exp(z.re) * Complex<T>(cos(z.im), sin(z.im));
}

/// exp(i x)
template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline Complex<T> expi(T a) {
    return Complex<T>(cos(a), sin(a));
}

/// log(z)
template <typename T>
inline Complex<T> log(Complex<T> z) {
    return Complex<T>(static_cast<T>(0.5) * log(z.squarenorm()), z.arg());
}

/// sqrt(z) branch cut at -x axis
template <typename T>
inline Complex<T> sqrt(Complex<T> z) {
    T r = z.squarenorm();
    T a = z.arg();
    return pow(r, 0.25) * expi(0.5 * a);
}

/// cbrt(z)
template <typename T>
inline Complex<T> cbrt(Complex<T> z) {
    T r = z.squarenorm();
    T a = z.arg();
    return pow(r, (1.0 / 6.0)) * expi((1.0 / 3.0) * a);
}

/// pow(z.p) = z^p = exp(p*log(z))
template <typename A, typename B>
inline auto pow(Complex<A> z, Complex<B> p) {
    return exp(p * log(z));
}

/// pow(z.p) with scalar power
template <typename T, typename S, std::enable_if_t<hila::is_arithmetic<S>::value, int> = 0>
inline Complex<T> pow(Complex<T> z, S p) {
    return exp(p * log(z));
}

/// pow(z.p) with scalar base
template <typename T, typename S, std::enable_if_t<hila::is_arithmetic<S>::value, int> = 0>
inline Complex<T> pow(S z, Complex<T> p) {
    return exp(p * log(z));
}

/// sin(z)
/// = sin(re + i im) = sin(re)cos(i im) + cos(re)sin(i im)
/// = sin(re) cosh(im) + i cos(re) sinh(im)
template <typename T>
inline Complex<T> sin(Complex<T> z) {
    return Complex<T>(sin(z.re) * cosh(z.im), cos(z.re) * sinh(z.im));
}

/// cos(z)
/// = cos(re)cos(i im) - sin(re)sin(i im) = cos(re)cosh(im) - i sin(re)sinh(im)
template <typename T>
inline Complex<T> cos(Complex<T> z) {
    return Complex<T>(cos(z.re) * cosh(z.im), -sin(z.re) * sinh(z.im));
}

/// tan(z) - rely on optimizer to simplify
template <typename T>
inline Complex<T> tan(Complex<T> z) {
    return sin(z) / cos(z);
}

/// sinh(z) = sinh(re)cosh(i im) + cosh(re)sinh(i im)
/// = sinh(re)cos(im) + i cosh(re)sin(im)
template <typename T>
inline Complex<T> sinh(Complex<T> z) {
    return Complex<T>(sinh(z.re) * cos(z.im), cosh(z.re) * sin(z.im));
}

/// cosh(z) = cosh(re)cosh(i im) - sinh(re)sinh(i im)
/// = cosh(re)cos(im) - i sinh(re)sin(im)
template <typename T>
inline Complex<T> cosh(Complex<T> z) {
    return Complex<T>(cosh(z.re) * cos(z.im), sinh(z.re) * sin(z.im));
}

/// tanh(z)
template <typename T>
inline Complex<T> tanh(Complex<T> z) {
    return sinh(z) / cosh(z);
}

/// arctan(z)
template <typename T>
inline Complex<T> atan(Complex<T> z) {
    return -0.5 * I * (log((1_i - z) / (1_i + z)));
}

/// arcsin(z)
template <typename T>
inline Complex<T> asin(Complex<T> z) {
    return -I * (log(I * z + sqrt(1 - z * z)));
}

/// arccos(z)
template <typename T>
inline Complex<T> acos(Complex<T> z) {
    return -I * (log(z + I * (sqrt(1 - z * z))));
}

/// artanh(z)
template <typename T>
inline Complex<T> atanh(Complex<T> z) {
    return 0.5 * log((1 + z) / (1 - z));
}

/// arsinh(z)
template <typename T>
inline Complex<T> asinh(Complex<T> z) {
    return log(z + sqrt(1 + z * z));
}

/// arcosh(z)
template <typename T>
inline Complex<T> acosh(Complex<T> z) {
    return log(z + sqrt(z * z - 1));
}


#endif
