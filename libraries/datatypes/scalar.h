#ifndef SCALAR_H
#define SCALAR_H

///////////////////////////////////////////////////////////
/// A scalar type. Templated, so that hilapp can convert
/// to a vector type
///////////////////////////////////////////////////////////
template <typename T = double> struct scalar {
    using base_type =  hila::number_type<T>;
    using argument_type = T;

    T value;

    // assignment is automatically OK, by c-standard
    //   scalar operator=(scalar rhs) {
    //     value = rhs.value;
    //     return *this;
    //   }
    scalar<T>() = default;

    scalar<T>(const scalar<T> &a) = default;

    // constructor from single scalar value
    template <typename scalar_t,
              std::enable_if_t<std::is_arithmetic<scalar_t>::value, int> = 0>
    constexpr scalar<T>(const scalar_t val) : value(static_cast<T>(val)) {}

    ~scalar<T>() = default;

    // automatic casting from scalar<T> -> scalar<A>
    // TODO: ensure this works if A is vector type!
    template <typename A> operator scalar<A>() const {
        return scalar<A>({static_cast<A>(value)});
    }

    // Conversion to T
    operator T() { return value; }

    template <typename scalar_t,
              std::enable_if_t<std::is_arithmetic<scalar_t>::value, int> = 0>
    scalar<T> &operator=(scalar_t s) {
        value = static_cast<T>(s);
        return *this;
    }

    T real() const { return value; }
    T imag() const { return 0; }

    T norm() const { return value * value; }
    // TODO: make this work for vector type!  Not double

    // currently this gives a compilation error
    double abs() const { return sqrt(static_cast<double>(norm())); }

    scalar<T> conj() const { return scalar<T>({value}); }

    // unary + and -
    scalar<T> operator+() const { return *this; }
    scalar<T> operator-() const { return scalar<T>(-value); }

    scalar<T> &operator+=(const scalar<T> &lhs) {
        value += lhs.value;
        return *this;
    }

    // TODO: for avx vector too -- #define new template macro
    template <typename A, std::enable_if_t<std::is_arithmetic<A>::value, int> = 0>
    scalar<T> &operator+=(const A &a) {
        value += static_cast<T>(a);
        return *this;
    }

    scalar<T> &operator-=(const scalar<T> &lhs) {
        value -= lhs.value;
        return *this;
    }

    // TODO: for vector too
    template <typename A, std::enable_if_t<std::is_arithmetic<A>::value, int> = 0>
    scalar<T> &operator-=(const A &a) {
        value -= static_cast<T>(a);
        return *this;
    }

    scalar<T> &operator*=(const scalar<T> &lhs) {
        value = value * lhs.value;
        return *this;
    }

    // TODO: for vector too
    template <typename A, std::enable_if_t<std::is_arithmetic<A>::value, int> = 0>
    scalar<T> &operator*=(const A &a) {
        value *= static_cast<T>(a);
        return *this;
    }

    // a/b
    scalar<T> &operator/=(const scalar<T> &lhs) {
        value = value / lhs.value;
        return *this;
    }

    // TODO: for vector too
    template <typename A, std::enable_if_t<std::is_arithmetic<A>::value, int> = 0>
    scalar<T> &operator/=(const A &a) {
        value /= static_cast<T>(a);
        return *this;
    }
};

template <typename T> scalar<T> operator+(const scalar<T> &a, const scalar<T> &b) {
    return scalar<T>(a.value + b.value);
}

// TODO: for avx vector too -- #define new template macro
template <typename T, typename A, std::enable_if_t<std::is_arithmetic<A>::value, int> = 0>
scalar<T> operator+(const scalar<T> &c, const A &a) {
    return scalar<T>(c.value + a);
}

template <typename T, typename A, std::enable_if_t<std::is_arithmetic<A>::value, int> = 0>
scalar<T> operator+(const A &a, const scalar<T> &c) {
    return scalar<T>(c.value + a);
}

// -
template <typename T> scalar<T> operator-(const scalar<T> &a, const scalar<T> &b) {
    return scalar<T>(a.value - b.value);
}

// TODO: for avx vector too -- #define new template macro
template <typename T, typename A, std::enable_if_t<std::is_arithmetic<A>::value, int> = 0>
scalar<T> operator-(const scalar<T> &c, const A &a) {
    return scalar<T>(c.value - a);
}

template <typename T, typename A, std::enable_if_t<std::is_arithmetic<A>::value, int> = 0>
scalar<T> operator-(const A &a, const scalar<T> &c) {
    return scalar<T>(a - c.value);
}

// *
template <typename T> scalar<T> operator*(const scalar<T> &a, const scalar<T> &b) {
    return scalar<T>(a.value * b.value);
}

// TODO: for avx vector too -- #define new template macro
template <typename T, typename A, std::enable_if_t<std::is_arithmetic<A>::value, int> = 0>
scalar<T> operator*(const scalar<T> &c, const A &a) {
    return scalar<T>(c.value * a);
}

template <typename T, typename A, std::enable_if_t<std::is_arithmetic<A>::value, int> = 0>
scalar<T> operator*(const A &a, const scalar<T> &c) {
    return scalar<T>(a * c.value);
}

// /   a/b = ab*/|b|^2
template <typename T> scalar<T> operator/(const scalar<T> &a, const scalar<T> &b) {
    return scalar<T>(a.value / b.value);
}

// TODO: for avx vector too -- #define new template macro
template <typename T, typename A, std::enable_if_t<std::is_arithmetic<A>::value, int> = 0>
scalar<T> operator/(const scalar<T> &c, const A &a) {
    return scalar<T>(c.value / a);
}

// a/c = ac*/|c|^2
template <typename T, typename A, std::enable_if_t<std::is_arithmetic<A>::value, int> = 0>
scalar<T> operator/(const A &a, const scalar<T> &c) {
    return scalar<T>(a / c.value);
}

#endif
