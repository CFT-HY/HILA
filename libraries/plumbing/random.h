#ifndef RANDOM_H_
#define RANDOM_H_

#include <cmath>
#include "plumbing/mersenne.h"
#include "plumbing/defs.h"

/// It is important that the random number generators gaussian_ran() and gaussian_ran2()
/// are marked as "loop function" and "contains rng", because hilapp does not have a view
/// inside them from all compilation units

#pragma hila contains_rng loop_function
double gaussian_ran();

#pragma hila contains_rng loop_function
double gaussian_ran2(double &out2);

// this defines gaussian_random(T & arg) for all double -> T convertible types
template <typename T, std::enable_if_t<std::is_convertible<double, T>::value, int> = 0>

inline void gaussian_random(T &d) {
    d = static_cast<T>(gaussian_ran());
}

// this does the same for std. random
template <typename T, std::enable_if_t<std::is_convertible<double, T>::value, int> = 0>

inline void random(T &d) {
    d = static_cast<T>(hila_random());
}

#endif