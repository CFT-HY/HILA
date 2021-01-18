#ifndef RANDOM_H_
#define RANDOM_H_

#include <cmath>
#include "plumbing/mersenne.h"
#include "plumbing/defs.h"

double gaussian_ran(double variance=0.5);
double gaussian_ran2(double & out2);
float gaussian_ran2(float & out2);

// this defines gaussian_random(T & arg) for all double -> T convertible types
template <typename T, std::enable_if_t<std::is_convertible<double,T>::value, int> = 0 >
#pragma hila loop_function
inline void gaussian_random(T & d) { d = static_cast<T>(gaussian_ran()); }

// this does the same for std. random
template <typename T, std::enable_if_t<std::is_convertible<double,T>::value, int> = 0 >
#pragma hila loop_function
inline void random(T & d) { d = static_cast<T>(hila_random()); }

#endif