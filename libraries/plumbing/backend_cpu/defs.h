#ifndef CPU_DEFS_H
#define CPU_DEFS_H

#include "plumbing/defs.h"

#define VANILLA

// Define random number generator
namespace hila {
    inline double random() { return mersenne(); }
    inline void seed_random(int seed) {

        /* First seed the generator */
        seed_mersenne(seed);

        /* "Warm up" to create a full state */
        for (int i = 0; i < 543210; i++)
            mersenne();
    }
} // namespace hila

// Trivial synchronization
inline void synchronize_threads() {}

/// Implements test for basic in types, similar to
/// std::is_arithmetic, but allows the backend to add
/// it's own basic tyes (such as AVX vectors)
template <class T>
struct is_arithmetic : std::integral_constant<bool, std::is_arithmetic<T>::value> {};

template <class T, class U>
struct is_assignable : std::integral_constant<bool, std::is_assignable<T, U>::value> {};

#endif
