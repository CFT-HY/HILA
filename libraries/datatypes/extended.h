/**
 * @file extended.h
 * @brief This files containts definitions for the extended precision class that allows for high precision reductions.
 */

#ifndef EXTENDED_H
#define EXTENDED_H

#include "plumbing/defs.h"

class Extended {
        
    public:
        double value;
        double compensation;
        Extended() = default;
        ~Extended() = default;
        Extended(double v) : value(v), compensation(0) {}
        Extended(double v, double c) : value(v), compensation(c) {}
        Extended(const Extended &rhs) : value(rhs.value), compensation(rhs.compensation) {}
        
        #pragma hila loop_function
        Extended &operator=(const Extended &rhs) {
            if (this != &rhs) {
                value = rhs.value;
                compensation = rhs.compensation;
            }
            return *this;
        }
        #pragma hila loop_function
        template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
        Extended &operator=(const T &v) {
            value = v;
            compensation = 0;
            return *this;
        }

        
        #pragma hila loop_function
        Extended operator+() const {
            return *this;
        }

        #pragma hila loop_function
        template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
        Extended operator+(const T& rhs) const {
            Extended temp(rhs);
            Extended result = *this;
            result += temp;
            return result;
        }

        #pragma hila loop_function
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
        #pragma hila loop_function
        Extended &operator+=(const Extended &rhs);

        #pragma hila loop_function
        template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
        inline Extended& operator+=(const T& rhs) {
            return *this += Extended(rhs);
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
        #pragma hila loop_function  
        operator double() const {
            return value + compensation;
        }

        #pragma hila loop_function
        double return_extended() const {
            return value + compensation;
        }
        
};

#pragma hila loop_function
inline Extended fast_two_sum(double a, double b) {
    double sum;
    double compensation;
    sum = a + b;
    if (abs(a) >= abs(b)) {
        compensation = b - (sum - a);
    } else {
        compensation = a - (sum - b);
    }
    return Extended(static_cast<double>(sum), compensation);
}

#pragma hila loop_function
inline Extended &Extended::operator+=(const Extended &rhs) {
    Extended temp = fast_two_sum(this->value, rhs.value);
    Extended temp2 = fast_two_sum(this->compensation, rhs.compensation);
    this->compensation = temp.compensation + temp2.value;
    this->value = temp.value;
    return *this;
}

inline std::ostream &operator<<(std::ostream &strm, const Extended &var) {
    return strm << var.value + var.compensation;
}


#endif // EXTENDED_H