#include "datatypes/extended.h"


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
