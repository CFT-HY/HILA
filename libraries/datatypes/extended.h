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
 *    s1.to_double()
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
    
    #pragma hila loop_function
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
    // need to take a copy to avoid problems with a += a
    inline const ExtendedPrecision &operator+=(ExtendedPrecision rhs) {
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
    inline const ExtendedPrecision &operator*=(T rhs) {
        static_assert(std::is_integral<T>::value,
                      "ExtendedPrecision can be multiplied only by integral value without losing "
                      "precision! Convert to double for more general operations");

        if (rhs == 0) {
            *this = 0;
        } else {
            bool negative = (!std::is_unsigned<T>::value && (rhs < 0));
            if (negative)
                rhs = -rhs;

            T mult = 1;
            while (rhs % 2 == 0) {
                rhs /= 2;
                mult *= 2;
            }

            // it's safe to multiply by power of 2
            this->value *= mult;
            this->compensation *= mult;
            this->compensation2 *= mult;

            auto a = *this;
            // the rest is done by addition
            // we add a rhs-1 times to *this, multiplying it by rhs
            // in any case, factorize rhs as 2^n + 2^m + .. + 1  (rhs is odd)
            while (rhs > 1) {
                mult = 1;
                while (2 * mult < rhs) {
                    mult *= 2;
                }
                ExtendedPrecision b;
                // again power of 2
                b.value = a.value * mult;
                b.compensation = a.compensation * mult;
                b.compensation2 = a.compensation2 * mult;

                *this += b;
                rhs -= mult;
            }
            if (negative) {
                this->value = -this->value;
                this->compensation = -this->compensation;
                this->compensation2 = -this->compensation2;
            }
        }
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
     * @details Converts the ExtendedPrecision object back to double by summing the compensation
     * and returning. Marked explicit in order to avoid mistakenly losing precision
     * @return double
     */
    explicit operator double() const {
        return value + (compensation + compensation2);
    }

    /**
     * @brief Returns the compensated value as double
     */
    double to_double() const {
        return value + (compensation + compensation2);
    }
};

/// template utilities: hila::is_extended<T>::value and
/// hila::is_arithmetic_or_extended<T>::value

namespace hila {
template <typename T>
struct is_extended : std::integral_constant<bool, std::is_same<T, ExtendedPrecision>::value> {};

template <typename T>
struct is_arithmetic_or_extended
    : std::integral_constant<bool, hila::is_arithmetic<T>::value || hila::is_extended<T>::value> {};

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
    static_assert(std::is_integral<T>::value,
                  "ExtendedPrecision can be multiplied only by integral value without losing "
                  "precision! Convert to double for more general operations");
    return a *= v;
}
template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
inline ExtendedPrecision operator*(T v, ExtendedPrecision a) {
    static_assert(std::is_integral<T>::value,
                  "ExtendedPrecision can be multiplied only by integral value without losing "
                  "precision! Convert to double for more general operations");
    return a *= v;
}


/// comparison operators

// Values can be equal even if the .value and .compensation are not identical!
// Identical values will satisfy the result below

template <typename A, typename B, std::enable_if_t<hila::AorB_extended<A, B>::value, int> = 0>
inline bool operator==(const A &a, const B &b) {
    ExtendedPrecision tmp = a - b;
    return (tmp.to_double() == 0.0);
}
template <typename A, typename B, std::enable_if_t<hila::AorB_extended<A, B>::value, int> = 0>
inline bool operator!=(const A &a, const B &b) {
    ExtendedPrecision tmp = a - b;
    return (tmp.to_double() != 0.0);
}
template <typename A, typename B, std::enable_if_t<hila::AorB_extended<A, B>::value, int> = 0>
inline bool operator>(const A &a, const B &b) {
    ExtendedPrecision tmp = a - b;
    return (tmp.to_double() > 0.0);
}
template <typename A, typename B, std::enable_if_t<hila::AorB_extended<A, B>::value, int> = 0>
inline bool operator<(const A &a, const B &b) {
    ExtendedPrecision tmp = a - b;
    return (tmp.to_double() < 0.0);
}
template <typename A, typename B, std::enable_if_t<hila::AorB_extended<A, B>::value, int> = 0>
inline bool operator>=(const A &a, const B &b) {
    ExtendedPrecision tmp = a - b;
    return (tmp.to_double() >= 0.0);
}
template <typename A, typename B, std::enable_if_t<hila::AorB_extended<A, B>::value, int> = 0>
inline bool operator<=(const A &a, const B &b) {
    ExtendedPrecision tmp = a - b;
    return (tmp.to_double() <= 0.0);
}


inline std::ostream &operator<<(std::ostream &strm, const ExtendedPrecision &var) {
    return strm << var.to_double();
}

namespace hila {

template <typename Ntype>
Ntype cast_to(const ExtendedPrecision &ep) {
    if constexpr (std::is_same<Ntype, ExtendedPrecision>::value)
        return ep;
    else
        return static_cast<Ntype>(ep.to_double());
}

inline std::string prettyprint(const ExtendedPrecision &val, int prec = 8) {
    return prettyprint(val.to_double(), prec);
}

} // namespace hila


#endif // EXTENDED_H