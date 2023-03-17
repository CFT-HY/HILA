#ifndef RANDOM_H_
#define RANDOM_H_

namespace hila {

constexpr bool device_rng_off = false;
constexpr bool device_rng_on = true;

/// Seed random generators with 64-bit unsigned
void seed_random(uint64_t seed, bool device_rng = true);

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

void seed_device_rng(unsigned long long seed);

/// free gpu rng generator state, random does not work after this (unless seeded again)
void free_gpu_rng();

bool is_device_rng_on();

} // namespace hila


#endif