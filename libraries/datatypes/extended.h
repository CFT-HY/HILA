/**
 * @file extended.h
 * @brief This files containts definitions for the extended precision class that allows for high
 * precision reductions.
 */

#ifndef _HILA_EXTENDED_H
#define _HILA_EXTENDED_H

#include "plumbing/defs.h"

class ExtendedPrecision {

  public:
    double value;
    double compensation;
    ExtendedPrecision() = default;
    ~ExtendedPrecision() = default;
    ExtendedPrecision(const ExtendedPrecision &rhs)
        : value(rhs.value), compensation(rhs.compensation) {}
    ExtendedPrecision(double v) : value(v), compensation(0) {}
    ExtendedPrecision(double v, double c) : value(v), compensation(c) {}

// Assignments are used in generated code, mark as loop functions
#pragma hila loop_function
    ExtendedPrecision &operator=(const ExtendedPrecision &rhs) {
        if (this != &rhs) {
            value = rhs.value;
            compensation = rhs.compensation;
        }
        return *this;
    }
#pragma hila loop_function
    template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
    ExtendedPrecision &operator=(const T &v) {
        value = v;
        compensation = 0;
        return *this;
    }


    ExtendedPrecision operator+() const {
        return *this;
    }

    template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
    ExtendedPrecision operator+(const T &rhs) const {
        ExtendedPrecision temp(rhs);
        ExtendedPrecision result = *this;
        result += temp;
        return result;
    }

    ExtendedPrecision operator-() const {
        return ExtendedPrecision(-value, -compensation);
    }

    static ExtendedPrecision fast_two_sum(double a, double b) {
        double sum = a + b;
        double compensation;
        if (abs(a) >= abs(b)) {
            compensation = b - (sum - a);
        } else {
            compensation = a - (sum - b);
        }
        return ExtendedPrecision(sum, compensation);
    }

    /**
     * @brief += addition assignment operator
     * @details Impmlements addition assignment with second order Neumaier algorithm.
     *
     * @param rhs
     * @return ExtendedPrecision&
     */
    inline const ExtendedPrecision &operator+=(const ExtendedPrecision &rhs) {
        ExtendedPrecision temp = fast_two_sum(this->value, rhs.value);
        ExtendedPrecision temp2 = fast_two_sum(this->compensation, rhs.compensation);
        this->compensation = temp.compensation + temp2.value;
        this->value = temp.value;
        return *this;
    }

    inline const ExtendedPrecision &operator-=(const ExtendedPrecision &rhs) {
        return *this += ExtendedPrecision(-rhs.value, -rhs.compensation);
    }

    template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
    inline const ExtendedPrecision &operator+=(const T &rhs) {
        return *this += ExtendedPrecision(rhs);
    }

    template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
    inline const ExtendedPrecision &operator-=(const T &rhs) {
        return *this += -rhs;
    }

    template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
    inline const ExtendedPrecision &operator*=(const T &rhs) {
        this->value *= rhs;
        this->compensation *= rhs;
        return *this;
    }

    inline const ExtendedPrecision &operator*=(const ExtendedPrecision &rhs) {
        ExtendedPrecision a(*this);
        *this *= rhs.value;
        a *= rhs.compensation;
        *this += a;
        return *this;
    }

    template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
    inline const ExtendedPrecision &operator/=(const T &rhs) {
        this->value /= rhs;
        this->compensation /= rhs;
        return *this;
    }

    // double double sum
    // ExtendedPrecision &operator+=(const ExtendedPrecision &rhs) {

    //     ExtendedPrecision<T> temp = fast_two_sum((*this).value, rhs.value);
    //     (*this) = fast_two_sum(temp.value, temp.compensation + (*this).compensation);

    //     return *this;

    // }

    /**
     * @brief Conversion back to double from ExtendedPrecision
     * @details Converts the ExtendedPrecision object back to double by summing the compensation and
     * returning. Marked explicit in order to avoid mistakenly losing precision
     * @return double
     */
    explicit operator double() const {
        return value + compensation;
    }

    /**
     * @brief Returns the compensated value as double
     */
    double get_value() const {
        return value + compensation;
    }
};

/// template utilities: hila::is_extended<T>::value and hila::is_arithmetic_or_extended<T>::value

namespace hila {
template <typename T>
struct is_extended : std::integral_constant<bool, std::is_same<T, ExtendedPrecision>::value> {};

template <typename T>
struct is_arithmetic_or_extended
    : std::integral_constant<bool, std::is_arithmetic<T>::value || hila::is_extended<T>::value> {};

template <typename A, typename B>
struct AorB_extended
    : std::integral_constant<
          bool, (hila::is_extended<A>::value && hila::is_arithmetic_or_extended<B>::value) ||
                    (hila::is_extended<B>::value && std::is_arithmetic<A>::value)> {};

} // namespace hila


/// ExtendedPrecision precision operators + - * / -- also ext * arithm.


// +
template <typename A, typename B, std::enable_if_t<hila::AorB_extended<A, B>::value, int> = 0>
inline ExtendedPrecision operator+(const A &a, const B &b) {
    ExtendedPrecision t(a);
    return t += b;
}

// -
template <typename A, typename B, std::enable_if_t<hila::AorB_extended<A, B>::value, int> = 0>
inline ExtendedPrecision operator-(const A &a, const B &b) {
    ExtendedPrecision t(a);
    return t -= b;
}

// *
template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
inline ExtendedPrecision operator*(const ExtendedPrecision &a, T v) {
    return ExtendedPrecision(a.value * v, a.compensation * v);
}
template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
inline ExtendedPrecision operator*(T v, const ExtendedPrecision &a) {
    return ExtendedPrecision(a.value * v, a.compensation * v);
}
inline ExtendedPrecision operator*(const ExtendedPrecision &a, const ExtendedPrecision &b) {
    return a * b.value + a * b.compensation;
}

// /
template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
inline ExtendedPrecision operator/(const ExtendedPrecision &a, T v) {
    return ExtendedPrecision(a.value / v, a.compensation / v);
}
inline ExtendedPrecision operator/(const ExtendedPrecision &a, const ExtendedPrecision &b) {
    auto divb = b.compensation / b.value;
    auto s = 1.0 - sqr(divb);
    return a / (b.value * s) * ExtendedPrecision(1.0, -divb);
}
template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
inline ExtendedPrecision operator/(T v, const ExtendedPrecision &b) {
    auto divb = b.compensation / b.value;
    auto s = 1.0 - sqr(divb);
    return v / (b.value * s) * ExtendedPrecision(1.0, -divb);
}

/// comparison operators

// Values can be equal even if the .value and .compensation are not identical!
// Identical values will satisfy the result below

template <typename A, typename B, std::enable_if_t<hila::AorB_extended<A, B>::value, int> = 0>
inline bool operator==(const A &a, const B &b) {
    ExtendedPrecision tmp = a - b;
    return (tmp.get_value() == 0.0);
}
template <typename A, typename B, std::enable_if_t<hila::AorB_extended<A, B>::value, int> = 0>
inline bool operator!=(const A &a, const B &b) {
    ExtendedPrecision tmp = a - b;
    return (tmp.get_value() != 0.0);
}
template <typename A, typename B, std::enable_if_t<hila::AorB_extended<A, B>::value, int> = 0>
inline bool operator>(const A &a, const B &b) {
    ExtendedPrecision tmp = a - b;
    return (tmp.get_value() > 0.0);
}
template <typename A, typename B, std::enable_if_t<hila::AorB_extended<A, B>::value, int> = 0>
inline bool operator<(const A &a, const B &b) {
    ExtendedPrecision tmp = a - b;
    return (tmp.get_value() < 0.0);
}
template <typename A, typename B, std::enable_if_t<hila::AorB_extended<A, B>::value, int> = 0>
inline bool operator>=(const A &a, const B &b) {
    ExtendedPrecision tmp = a - b;
    return (tmp.get_value() >= 0.0);
}
template <typename A, typename B, std::enable_if_t<hila::AorB_extended<A, B>::value, int> = 0>
inline bool operator<=(const A &a, const B &b) {
    ExtendedPrecision tmp = a - b;
    return (tmp.get_value() <= 0.0);
}


inline std::ostream &operator<<(std::ostream &strm, const ExtendedPrecision &var) {
    return strm << var.value + var.compensation;
}

namespace hila {
inline std::string prettyprint(const ExtendedPrecision & val, int prec = 8) {
    return prettyprint(val.get_value(), prec);
}

} // namespace hila


#endif // EXTENDED_H