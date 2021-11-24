#ifndef RANDOM_H_
#define RANDOM_H_

#include <cmath>
#include "plumbing/mersenne.h"

/// It is important that the random number generators hila::gaussrand() and gaussrand2()
/// are marked as "loop function" and "contains rng", because hilapp does not have a view
/// inside them from all compilation units

namespace hila {

#pragma hila contains_rng loop_function
double gaussrand();

#pragma hila contains_rng loop_function
double gaussrand2(double &out2);

} // namespace hila


#endif