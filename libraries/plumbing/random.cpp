#include <cmath>

#define IN_HILA_RANDOM 

#include "hila.h"

/////////////////////////////////////////////////////////////////////////

static bool rng_is_initialized = false;


#ifdef USE_PHILOX_RNG
// define global vars for philox

// global index for this node/rank for non-onsites code
static uint64_t hila_philox_my_index;

#else
// Now philox is not defined

#include <random>

// static variable which holds the random state
// Use 64-bit mersenne twister
static std::mt19937_64 mersenne_twister_gen;

// random numbers are in interval [0,1)
static std::uniform_real_distribution<double> real_rnd_dist(0.0, 1.0);

// #endif


// In GPU code hila::random() defined in hila_gpu.cpp
#if !defined(CUDA) && !defined(HIP)
double hila::random() {
    return real_rnd_dist(mersenne_twister_gen);
}

#endif

// Generate random number in non-kernel (non-loop) code.  Not meant to
// be used in "user code"
double hila::host_random() {
    return real_rnd_dist(mersenne_twister_gen);
}


/////////////////////////////////////////////////////////////////////////


namespace hila {

/**
 * @brief Random shuffling of rng seed for MPI nodes.
 * @details Do it in a manner makes it difficult to give the same seed by mistake
 * and also avoids giving the same seed for 2 nodes
 * For single MPI node seed remains unchanged
 */
uint64_t shuffle_rng_seed(uint64_t seed) {

    uint64_t n = hila::myrank();
    if (hila::partitions.number() > 1)
        n += hila::partitions.mylattice() * hila::number_of_nodes();

    return (seed + n) ^ (n << 31);
}


/**
 *@param seed unsigned int_64  
 *@note On MPI, this function shuffles different seed values for all MPI ranks.  
 */
void initialize_host_rng(uint64_t seed) {

    seed = hila::shuffle_rng_seed(seed);

    // #if !defined(OPENMP)
    mersenne_twister_gen.seed(seed);
    // warm it up
    for (int i = 0; i < 9000; i++)
        mersenne_twister_gen();
    // #endif
}

} // namespace hila

// philox defns
#endif

/**
 * @details The optional 2nd argument indicates whether to initialize the RNG on GPU device:
 * `hila::device_rng::on` (default) or `hila::device_rng::off`.  This argument does nothing if no GPU
 * platform.  If `hila::device_rng::off` is used, `onsites()` -loops cannot contain random number calls
 * (Runtime error will be flagged and program exits).
 * 
 * 
 * Seed is shuffled so that different nodes
 * get different rng seeds.  If `seed == 0`,
 * seed is generated through using the `time()` -function.
 */  
void hila::seed_random(uint64_t seed, bool device_init) {

    rng_is_initialized = true;

    if (!lattice.is_initialized()) {
        hila::error("lattice.setup() must be called before hila::seed_random()");
    }


    uint64_t n = hila::myrank();
    if (hila::partitions.number() > 1)
        n += hila::partitions.mylattice() * hila::number_of_nodes();

    if (seed == 0) {
        // get seed from time
        if_rank0() {
            struct timespec tp;

            clock_gettime(CLOCK_MONOTONIC, &tp);
            seed = tp.tv_sec;
            seed = (seed << 30) ^ tp.tv_nsec;
            hila::out0 << "Random seed from time: " << seed << '\n';
        }
        hila::broadcast(seed);
    }

    if (hila::partitions.number() > 1)
        seed = seed ^ ((static_cast<uint64_t>(hila::partitions.mylattice())) << 28);

#ifndef USE_PHILOX_RNG

    hila::out0 << "Using old (deprecated) random number generators\n";
    
    hila::initialize_host_rng(seed);

#if defined(CUDA) || defined(HIP)

    // we can use the same seed, the generator is different
    if (device_init) {
        hila::initialize_device_rng(seed);
    } else {
        hila::out0 << "Not initializing GPU random numbers\n";
    }

#endif

#else

    // Now use PHILOX
    philox_seed0 = static_cast<uint32_t>(seed);
    philox_seed1 = static_cast<uint32_t>(seed >> 32);

    hila::out0 << " SEED IS " << philox_seed0() << philox_seed1() << '\n';

    hila_philox_loop_counter = 0;
    hila_philox_my_index = std::numeric_limits<uint64_t>::max() - 1 - hila::myrank();

    hila::philox_setup(hila_philox_loop_counter, hila_philox_my_index);

    hila::out0 << "Using Philox 4x32-10 random number generator\n";

#endif
    hila::out0 << "RNG seed for node 0: " << seed << std::endl;

}


#ifdef USE_PHILOX_RNG

void hila::philox_after_onsites() {
    hila::philox_setup(++hila_philox_loop_counter, hila_philox_my_index);
}

#endif

////////////////////////////////////////////////////////////////////
// Def here gpu rng functions for non-gpu
////////////////////////////////////////////////////////////////////

#if !(defined(CUDA) || defined(HIP)) && !defined(USE_PHILOX_RNG)

/**
 *@details `hila::random()` does not work inside `onsites()` after this, 
 *unless seeded again using `initialize_device_rng()`. Frees the memory RNG takes on the device. 
 */  
void hila::free_device_rng() {}


/**
 *@details Returns `true` on non-GPU archs.
 */
bool hila::is_device_rng_on() {
    return rng_is_initialized;
}

/**
 *@details This function shuffles the seed for different MPI ranks on MPI.  Called by `seed_random()` unless its 2nd
 *argument is `hila::device_rng_off`. This can reinitialize device RNG free'd by `free_device_rng()`.
 */
void hila::initialize_device_rng(uint64_t seed) {}

#endif


#define VARIANCE 1.0

double hila::gaussrand2(out_only double &out2) {

    double urnd;
    double phi = 2.0 * M_PI * hila::random2(urnd);

    // this should not really trigger
    while (urnd <= 0.0 || urnd > 1.0) {
        urnd = hila::random();
    }

    double r = sqrt(-::log(urnd) * (2.0 * VARIANCE));
    out2 = r * cos(phi);
    return r * sin(phi);
}

//#if !defined(CUDA) && !defined(HIP)
#if 0

/**
 * @details By default these gives random numbers with variance \f$1.0\f$ and expectation value \f$0.0\f$, i.e.
 * \f[
 *    e^{-(\frac{x^{2}}{2})}
 * \f]
 * with variance
 * \f[
 *    < x^{2}-0.0> = 1
 * \f] 
 *
 * If you want random numbers with variance \f$ \sigma^{2} \f$, multiply the
 * result by \f$ \sqrt{\sigma^{2}} \f$ i.e.,
 * \code {.cpp}
 *       sqrt(sigma * sigma) * gaussrand();
 * \endcode
 *
 * @return double
 */  
double hila::gaussrand() {
    static double second;
    static bool draw_new = true;
    if (draw_new) {
        draw_new = false;
        return hila::gaussrand2(second);
    }
    draw_new = true;
    return second;
}

#else

// Cuda and other stuff which does not accept
// static variables - just throw away another gaussian number.

// #pragma hila loop function contains rng
double hila::gaussrand() {
    double second;
    return hila::gaussrand2(second);
}

#endif

/**
 *@returns bool
 */
bool hila::is_rng_seeded() {
    return rng_is_initialized;
}


///////////////////////////////////////////////////////////////
// RNG initialization check - emitted on loops
///////////////////////////////////////////////////////////////

/**
 *@details program quits with error message if RNG is not initialized.
 *It also quit with error messages if the device RNG is not initialized.
 */
void hila::check_that_rng_is_initialized() {

    if (!rng_is_initialized) {
        hila::out0 << "ERROR: trying to use random numbers without initialization"
                   << std::endl;
        hila::terminate(1);
    }
#if !defined(USE_PHILOX_RNG) && (defined(CUDA) || defined(HIP))
    if (!hila::is_device_rng_on()) {
        hila::out0 << "ERROR: GPU random number generator is not initialized and onsites()-loop is "
                      "using random numbers"
                   << std::endl;
        hila::terminate(1);
    }
#endif
}






