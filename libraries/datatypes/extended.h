/**
 * @file extended.h
 * @brief This files containts definitions for the extended precision class that allows for high precision reductions.
 */

#ifndef EXTENDED_H
#define EXTENDED_H

#include "plumbing/defs.h"

template <typename T>
class Extended;


template <typename T>
Extended<T> fast_two_sum(T a, T b);

template <typename T>
class Extended {
    public:
        using base_type = hila::arithmetic_type<T>;
        using argument_type = T;
        T value, compensation;
        Extended() : value(0), compensation(0) {}
        Extended(const T &v) : value(v), compensation(0) {}
        Extended(const T &v, const T &c) : value(v), compensation(c) {}
        Extended(const Extended &rhs) : value(rhs.value), compensation(rhs.compensation) {}
        Extended &operator=(const Extended &rhs) {
            if (this != &rhs) {
                value = rhs.value;
                compensation = rhs.compensation;
            }
            return *this;
        }
        Extended &operator=(const T &v) {
            value = v;
            compensation = 0;
            return *this;
        }
        
        Extended operator+() const {
            return *this;
        }

        Extended operator-() const {
            return Extended(-value, -compensation);
        }
        /**
         * @brief += addition assignment operator
         * @details Impmlements addition assignment with second order Neumaier algorithm. 
         * 
         * @param rhs 
         * @return Extended& 
         */
        Extended &operator+=(const Extended &rhs) {

            Extended<T> temp = fast_two_sum((*this).value, rhs.value);
            Extended<T> temp2 = fast_two_sum((*this).compensation, rhs.compensation);
            (*this).compensation = temp.compensation + temp2.value;
            (*this).value = temp.value;

            return *this;

        }
        // double double sum
        // Extended &operator+=(const Extended &rhs) {

        //     Extended<T> temp = fast_two_sum((*this).value, rhs.value);
        //     (*this) = fast_two_sum(temp.value, temp.compensation + (*this).compensation);

        //     return *this;

        // }

        /**
         * @brief Conversion back to arithmetic internal type Extended<T> 
         * @details Converts the Extended<T> object back to the base arithmetic type T by summing the compensation and returning.
         * @return T 
         */
        operator T() const {
            return value + compensation;
        }

        T return_extended() const {
            return value + compensation;
        }
        
};

template <typename T>
Extended<T> fast_two_sum(T a, T b) {
    Extended<T> result;
    result.value = a + b;
    if (abs(a) >= abs(b)) {
        result.compensation = b - (result.value - a);
    } else {
        result.compensation = a - (result.value - b);
    }
    return result;
}

template <typename T>
std::ostream &operator<<(std::ostream &strm, const Extended<T> &var) {
    return strm << var.value + var.compensation;
}


#endif // EXTENDED_H