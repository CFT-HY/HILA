#include <cmath>

#define IN_HILA_RANDOM 

#include "hila.h"

/////////////////////////////////////////////////////////////////////////

static bool rng_is_initialized = false;

// define global vars for philox

// global index for this node/rank for non-onsites code
static uint64_t hila_philox_my_index;

/**
 * philox_after_onsites() must be called after each onsites()-block where RNGs are used
 * This is installed by hilapp, not by user
 */

void hila::philox_after_onsites() {
    hila::philox_setup(++hila_philox_loop_counter, hila_philox_my_index);
}

/**
 * @brief Seed the RNG, must be called before RNG is used
 *
 * @details Philox random numbers are identical independent of the number of MPI ranks or
 * GPU/CPU (inside onsites()-loops).  If `seed == 0`, seed is generated using the clock.
 * The 2nd argument is not used, it is kept for compatibility.
 */
void hila::seed_random(uint64_t seed, bool device_init) {

    rng_is_initialized = true;

    if (!lattice.is_initialized()) {
        hila::error("lattice.setup() must be called before hila::seed_random()");
    }

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

    // Now use PHILOX
    philox_seed0 = static_cast<uint32_t>(seed);
    philox_seed1 = static_cast<uint32_t>(seed >> 32);

    hila_philox_loop_counter = 0;
    hila_philox_my_index = std::numeric_limits<uint64_t>::max() - hila::myrank();

    hila::philox_setup(hila_philox_loop_counter, hila_philox_my_index);

    hila::out0 << "Using Philox 4x32-10 random number generator\n";

    hila::out0 << "Seed is " << seed << std::endl;
}



#define VARIANCE 1.0

/**
 * @brief gaussrand() and gaussrand(val2) return gaussian distribution double prec numbers
 * 
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

/**
 * @brief First, 2 gaussian random numbbers
 */

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
#if 0 // comment out 

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
 *@brief check if rng initialization has been done
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
}






