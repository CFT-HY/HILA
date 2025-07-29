#ifndef HILA_EXTENDED_H
#define HILA_EXTENDED_H

/**
 * @file extended.h
 * @brief This files containts definitions for the extended precision class that allows for high
 * precision reductions using Kahan-Babushka-Klein summation.
 *
 * @details
 * Uses internally three double precision variables to store the value and compensation. Precision
 * can be > 15 orders of magnitude better than plain doble.
 *
 * Use:
 * @code {.cpp}
 *
 * ExtendedPrecision s1 = 0, s2 = 0;
 * onsites(ALL) {
 *     s1 += <double/float/int -valued expression>;
 *     s2 += <some other expression>
 * }
 *
 * // below subtraction and comparison is done in extended precision
 * if (s1 == s2)
 *     hila::out0 << "s1 equals s2\n";
 * else
 *     hila::out0 << "s1 is " << (s1 > s2) ? "greater" : "less"
 *                << " than s2, with difference " << s1 - s2 << '\n';
 * @endcode
 *
 * Class implements basic arithmetics ( + - * / ), assignments (= += *=) and
 * comparison ops for ExtendedPrecision types. ExtendedPrecision is not
 * automatically downgraded to double in order to avoid accidental loss of accuracy.
 *
 * To get a double approximation of the value use
 * @code {.cpp}
 *    s1.get_value()
 *    // or equvalently explicit casting
 *    static_cast<double>(s1)
 * @endcode
 *
 */


#include "plumbing/defs.h"

class ExtendedPrecision {

  public:
    double value;
    double compensation;
    double compensation2;

    ExtendedPrecision() = default;
    ~ExtendedPrecision() = default;
    ExtendedPrecision(const ExtendedPrecision &rhs)
        : value(rhs.value), compensation(rhs.compensation), compensation2(rhs.compensation2) {}
    ExtendedPrecision(double v) : value(v), compensation(0), compensation2(0) {}
    ExtendedPrecision(double v, double c, double c2)
        : value(v), compensation(c), compensation2(c2) {}

// Assignments are used in generated code, mark as loop functions
#pragma hila loop_function
    ExtendedPrecision &operator=(const ExtendedPrecision &rhs) {
        if (this != &rhs) {
            value = rhs.value;
            compensation = rhs.compensation;
            compensation2 = rhs.compensation2;
        }
        return *this;
    }
#pragma hila loop_function
    template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
    ExtendedPrecision &operator=(const T &v) {
        value = v;
        compensation = compensation2 = 0;
        return *this;
    }


    ExtendedPrecision operator+() const {
        return *this;
    }

    // template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
    // ExtendedPrecision operator+(const T &rhs) const {
    //     ExtendedPrecision temp(rhs);
    //     ExtendedPrecision result = *this;
    //     result += temp;
    //     return result;
    // }

    inline ExtendedPrecision operator-() const {
        return ExtendedPrecision(-value, -compensation, -compensation2);
    }

    /**
     * @brief += addition assignment operator
     * @details Impmlements the Kahan-Babushka-Klein summation
     *
     * @param rhs - value to be summed
     * @return ExtendedPrecision&
     */

#pragma hila loop_function
    template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
    inline const ExtendedPrecision &operator+=(const T &rhs_) {
        double rhs = rhs_;
        double t = value + rhs;
        double c, cc;
        if (abs(value) >= abs(rhs)) {
            c = (value - t) + rhs;
        } else {
            c = (rhs - t) + value;
        }
        value = t;
        t = compensation + c;
        if (abs(compensation) >= abs(c)) {
            cc = (compensation - t) + c;
        } else {
            cc = (c - t) + compensation;
        }
        compensation = t;
        compensation2 += cc;

        return *this;
    }


    /**
     * @brief += addition assignment operator for two EPs
     */

#pragma hila loop_function
    inline const ExtendedPrecision &operator+=(const ExtendedPrecision &rhs) {
        *this += rhs.compensation2;
        *this += rhs.compensation;
        *this += rhs.value;
        return *this;
    }

    inline const ExtendedPrecision &operator-=(const ExtendedPrecision &rhs) {
        return *this += -rhs;
    }

    template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
    inline const ExtendedPrecision &operator-=(const T &rhs) {
        return *this += -rhs;
    }

    template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
    inline const ExtendedPrecision &operator*=(const T &rhs) {
        this->value *= rhs;
        this->compensation *= rhs;
        this->compensation2 *= rhs;
        return *this;
    }

    inline const ExtendedPrecision &operator*=(const ExtendedPrecision &rhs) {
        ExtendedPrecision a(*this);
        ExtendedPrecision b(*this);
        *this *= rhs.value;
        a *= rhs.compensation;
        *this += a;
        b *= rhs.compensation2;
        *this += b;
        return *this;
    }

    template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
    inline const ExtendedPrecision &operator/=(const T &rhs) {
        *this *= 1.0 / rhs;
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
        return value + (compensation + compensation2);
    }

    /**
     * @brief Returns the compensated value as double
     */
    double get_value() const {
        return value + (compensation + compensation2);
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

// use here copy constructor on purpose
// +
template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
inline ExtendedPrecision operator+(ExtendedPrecision a, const T b) {
    return a += b;
}
template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
inline ExtendedPrecision operator+(const T b, ExtendedPrecision a) {
    return a += b;
}
inline ExtendedPrecision operator+(ExtendedPrecision a, const ExtendedPrecision &b) {
    return a += b;
}

// -
template <typename A, typename B, std::enable_if_t<hila::AorB_extended<A, B>::value, int> = 0>
inline ExtendedPrecision operator-(const A &a, const B &b) {
    return a + (-b);
}

// *
template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
inline ExtendedPrecision operator*(ExtendedPrecision a, T v) {
    return a *= v;
}
template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
inline ExtendedPrecision operator*(T v, ExtendedPrecision a) {
    return a *= v;
}
inline ExtendedPrecision operator*(ExtendedPrecision a, const ExtendedPrecision &b) {
    return a *= b;
}

// /
template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
inline ExtendedPrecision operator/(ExtendedPrecision a, T v) {
    return a /= v;
}
// division by EP may lose some accuracy
template <typename T, std::enable_if_t<hila::is_arithmetic_or_extended<T>::value, int> = 0>
inline ExtendedPrecision operator/(const T &a, const ExtendedPrecision &b) {
    auto divb = (b.compensation + b.compensation2) / b.value;
    auto s = 1.0 - sqr(divb);
    return a / (b.value * s) *
           ExtendedPrecision(1.0, -b.compensation / b.value, -b.compensation2 / b.value);
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
    return strm << var.get_value();
}

namespace hila {
inline std::string prettyprint(const ExtendedPrecision &val, int prec = 8) {
    return prettyprint(val.get_value(), prec);
}

} // namespace hila


#endif // EXTENDED_H