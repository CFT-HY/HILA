#ifndef RANDOM_H_
#define RANDOM_H_

#include <cmath>
#include "plumbing/mersenne.h"
#include "plumbing/defs.h"

/// It is important that the random number generators hila::gaussrand() and gaussrand2()
/// are marked as "loop function" and "contains rng", because hilapp does not have a view
/// inside them from all compilation units

namespace hila {
#pragma hila contains_rng loop_function
double gaussrand();

#pragma hila contains_rng loop_function
double gaussrand2(double &out2);
} // namespace hila

// this defines gaussian_random(T & arg) for all double -> T convertible types
template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
inline void gaussian_random(T &d) {
    d = static_cast<T>(hila::gaussrand());
}

// this does the same for std. random
template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>

inline void random(T &d) {
    d = static_cast<T>(hila::random());
}

#endif