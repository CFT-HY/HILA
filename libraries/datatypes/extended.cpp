
#include "plumbing/defs.h"
#include "datatypes/extended.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>


/**
 * @brief += addition assignment operator
 * @details Impmlements the Kahan-Babushka-Klein summation
 *
 * @param rhs - value to be summed
 * @return ExtendedPrecision&
 */

#pragma hila loop_function
const ExtendedPrecision &ExtendedPrecision::operator+=(double rhs) {
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


/////////////////////////////////////////////////////////////////////////////////////////
/// Decimal output in double-double precision.  These rely on exact IEEE rounding,
/// thus must not be compiled with reassociating optimizations (see above).
/////////////////////////////////////////////////////////////////////////////////////////

namespace {

struct dd_real {
    double hi, lo;
};

// s + e == a + b exactly
inline void two_sum(double a, double b, double &s, double &e) {
    s = a + b;
    double bb = s - a;
    e = (a - (s - bb)) + (b - bb);
}

// as above, assuming |a| >= |b|
inline void quick_two_sum(double a, double b, double &s, double &e) {
    s = a + b;
    e = b - (s - a);
}

// p + e == a * b exactly
inline void two_prod(double a, double b, double &p, double &e) {
    p = a * b;
    e = std::fma(a, b, -p);
}

dd_real dd_add(dd_real a, dd_real b) {
    double s, e, t, f;
    two_sum(a.hi, b.hi, s, e);
    two_sum(a.lo, b.lo, t, f);
    e += t;
    quick_two_sum(s, e, s, e);
    e += f;
    quick_two_sum(s, e, s, e);
    return {s, e};
}

dd_real dd_mul(dd_real a, dd_real b) {
    double p, e;
    two_prod(a.hi, b.hi, p, e);
    e += a.hi * b.lo + a.lo * b.hi;
    quick_two_sum(p, e, p, e);
    return {p, e};
}

dd_real dd_div(dd_real a, dd_real b) {
    double q1 = a.hi / b.hi;
    dd_real r = dd_add(a, dd_mul({-q1, 0.0}, b));
    double q2 = r.hi / b.hi;
    r = dd_add(r, dd_mul({-q2, 0.0}, b));
    double q3 = r.hi / b.hi;
    quick_two_sum(q1, q2, q1, q2);
    return dd_add({q1, q2}, {q3, 0.0});
}

// 10^n, n >= 0; exact up to n ~ 45
dd_real dd_pow10(int n) {
    dd_real r = {1.0, 0.0}, p = {10.0, 0.0};
    while (n > 0) {
        if (n & 1)
            r = dd_mul(r, p);
        p = dd_mul(p, p);
        n >>= 1;
    }
    return r;
}

// x * 10^n, done in steps so that the power of 10 does not overflow
dd_real dd_scale10(dd_real x, int n) {
    while (n != 0) {
        int step = std::max(-280, std::min(280, n));
        if (step > 0)
            x = dd_mul(x, dd_pow10(step));
        else
            x = dd_div(x, dd_pow10(-step));
        n -= step;
    }
    return x;
}

} // namespace

/**
 * @brief Decimal representation of the value in scientific notation, e.g.
 * -6.38876966070759031234567890e-14
 * @details The value is first renormalized exactly to double-double (~32 significant
 * digits), which is then converted to decimal.  Digits beyond ~31 are not meaningful.
 * @param digits - number of significant digits, 1 ... 34
 */
std::string ExtendedPrecision::to_string(int digits) const {
    digits = std::max(1, std::min(34, digits));

    if (!std::isfinite(to_double())) {
        std::ostringstream ss;
        ss << to_double();
        return ss.str();
    }

    // renormalize exactly: hi + mid (+ lo, below double-double precision)
    double t, e1, e2, hi, mid, lo;
    two_sum(compensation, compensation2, t, e2);
    two_sum(value, t, hi, e1);
    two_sum(e1, e2, mid, lo);
    quick_two_sum(hi, mid, hi, mid);

    std::string s;
    if (hi < 0 || (hi == 0 && std::signbit(value))) {
        s = "-";
        hi = -hi;
        mid = -mid;
    }

    std::vector<int> d(digits, 0);
    int exp10 = 0;

    if (hi != 0) {
        // scale to [1, 10)
        exp10 = static_cast<int>(std::floor(std::log10(hi)));
        dd_real x = dd_scale10({hi, mid}, -exp10);
        // x.hi may be exactly 1 or 10 while x.lo has the opposite sign
        while (x.hi > 10.0 || (x.hi == 10.0 && x.lo >= 0)) {
            x = dd_div(x, {10.0, 0.0});
            exp10++;
        }
        while (x.hi < 1.0 || (x.hi == 1.0 && x.lo < 0)) {
            x = dd_mul(x, {10.0, 0.0});
            exp10--;
        }

        for (int i = 0; i < digits; i++) {
            double di = std::floor(x.hi);
            x = dd_add(x, {-di, 0.0});
            // x.hi may be rounded to an integer while x.lo has the opposite sign
            if (x.hi < 0) {
                di -= 1;
                x = dd_add(x, {1.0, 0.0});
            } else if (x.hi >= 1.0) {
                di += 1;
                x = dd_add(x, {-1.0, 0.0});
            }
            d[i] = std::max(0, std::min(9, static_cast<int>(di)));
            x = dd_mul(x, {10.0, 0.0});
        }

        // round to nearest using the remainder
        if (x.hi >= 5.0) {
            int i = digits - 1;
            while (i >= 0 && d[i] == 9) {
                d[i] = 0;
                i--;
            }
            if (i >= 0) {
                d[i]++;
            } else {
                d[0] = 1;
                exp10++;
            }
        }
    }

    s += static_cast<char>('0' + d[0]);
    if (digits > 1) {
        s += '.';
        for (int i = 1; i < digits; i++)
            s += static_cast<char>('0' + d[i]);
    }

    char ebuf[16];
    std::snprintf(ebuf, sizeof(ebuf), "e%+03d", exp10);
    s += ebuf;
    return s;
}
