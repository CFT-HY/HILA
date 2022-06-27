#ifndef RANDOM_H_
#define RANDOM_H_

namespace hila {

/// Seed random generators with 64-bit unsigned
void seed_random(uint64_t seed);

/// It is important that the random number generators random(), gaussrand() and gaussrand2()
/// are marked as "loop function" and "contains rng", because hilapp does not have a view
/// inside them from all compilation units - although hila::random() has special treatment


#pragma hila contains_rng loop_function
double random();

// alias for host (non-loop) random number, used in GPU code
double host_random();

#pragma hila contains_rng loop_function
double gaussrand();

#pragma hila contains_rng loop_function
double gaussrand2(double &out2);

// do what the name says
void check_that_rng_is_initialized();

} // namespace hila


#endif