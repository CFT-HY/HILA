#ifndef CPU_DEFS_H
#define CPU_DEFS_H

#include "../../plumbing/defs.h"

#define VANILLA


// Define random number generator
#define seed_random(seed) seed_mersenne(seed)
inline double hila_random(){ return mersenne(); }

// Trivial synchronization
inline void synchronize_threads(){}


/// Implements test for basic in types, similar to 
/// std::is_arithmetic, but allows the backend to add
/// it's own basic tyes (such as AVX vectors)
template< class T >
struct is_arithmetic : std::integral_constant<
  bool,
  std::is_arithmetic<T>::value
> {};



#endif
