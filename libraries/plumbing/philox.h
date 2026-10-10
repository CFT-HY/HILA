#ifndef HILA_PHILOX_H_
#define HILA_PHILOX_H_

#include <array>
#include <cstdint>
#include <iostream>
#include <limits>

#include "plumbing/globals.h"

/**
 * @file philox.h
 * @brief This contains definition of the philox 4x32-10 pseudorandom number generator
 * 
 * @details 
 * Philox is a stateless RNG, based on "cryptographic" function P_s:
 * R = P_s(C), where input C is a 128-bit value, counter, and R 128-bit encryption of C,
 * which is now taken as a 128-bit random number.  The function P_s also depends on 
 * a 64-bit seed s, different seeds producing a different sequence.
 * 
 * The function P_s is implemented below in philox4x32_10.  The 128-bit values are
 * implemented as 4 32-bit values (of type uint32_t).
 * (see https://www.thesalmons.org/john/random123/papers/random123sc11.pdf)
 * 
 * From the 128-bit output we can construct e.g. 1 or 2 double precision values.
 * 
 * In hila the counter C is calculated as follows:
 * 
 * a: in onsites(), on each site:
 *    - 64 bits for SiteIndex (x + nx*(y + ny*(z + nz*t)))
 *    - 32 bits bits for "global" counter, which is incremented by one for each 
 *      onsites-loop where random numbers are used 
 *    - 32 bits for "local" counter for each site, which is set to 0 at the beginning
 *      of onsites-block and incremented each time random number generator is called.
 *    This scheme guarantees unique C for each call of the RNG on each site.
 * 
 * b: outside onsites():
 *    - 64-bit SiteIndex is substituted by (max uint64_t - hila::myrank())
 *    - the same 32-bit "global" counter as above
 *    - 32-bit counter which is incremented for each RNG call between onsites-blocks.
 * 
 * The above scheme guarantees identical random numbers for each lattice size independent of the 
 * node division or computing platform.
 * 
 * [with the exception of random values generated outside onsites in nodes !=0, which naturally
 *  depend on the existence (number) of nodes.  It is recommended that random values outside onsites
 *  are only generated on node 0 if independence on node division is the goal.]
 * 
 * On GPUs this uses (3*N_threads + 1)*sizeof(uint32_t) bytes __shared__ memory as a scratchpad
 * in onsites-loops which use RNG. If RNG is not used no __shared__ memory is used.
 */



// Philox4x32 constant multipliers
#define PHILOX_M0 0xD2511F53U
#define PHILOX_M1 0xCD9E8D57U

// Philox4x32 key Weyl-sequence
#define PHILOX_W0 0x9E3779B9U
#define PHILOX_W1 0xBB67AE85U

#ifdef IN_HILA_RANDOM
#define HILA_PHILOX_EXTERN /* nothing */
#else
#define HILA_PHILOX_EXTERN extern
#endif

// store seed in constant memory
HILA_PHILOX_EXTERN hila::global<uint32_t> philox_seed0;
HILA_PHILOX_EXTERN hila::global<uint32_t> philox_seed1;

// global loop counter for philox, increased after onsites() which contains random
HILA_PHILOX_EXTERN uint32_t hila_philox_loop_counter;


namespace hila {

/**
 * @brief philox 4x32-10 -function - generates 4 random uint32_t:s from input "counters"
 * Both counters and output are handled by the src-array
 *
 * @note This is a modified version of Philox generator from Random123 library.
 * This uses 4x 32-bit counter, 2x 32-bit key along with 10 rounds.
 *
 * @internal this is not meant to be directly used by the user
 */

#pragma hila loop_function novector
inline void philox4x32_10(uint32_t (&src)[4]) {

    // let hilapp see this function, it has to modify global var references
    uint32_t k0 = philox_seed0();
    uint32_t k1 = philox_seed1();

    // Unroll as much as you can
    #pragma unroll
    for (int round = 0; round < 10; ++round) {
        // 1. From 32x32 -> 64-bits multiplications
        uint64_t prod0 = static_cast<uint64_t>(src[0]) * PHILOX_M0;
        uint64_t prod1 = static_cast<uint64_t>(src[2]) * PHILOX_M1;

        // 2. Divide into hi and lo
        uint32_t hi0 = static_cast<uint32_t>(prod0 >> 32);
        uint32_t lo0 = static_cast<uint32_t>(prod0);

        uint32_t hi1 = static_cast<uint32_t>(prod1 >> 32);
        uint32_t lo1 = static_cast<uint32_t>(prod1);

        // 3. Philox-permutation and XOR-mixing
        uint32_t next_x0 = hi1 ^ src[1] ^ k0;
        uint32_t next_x1 = lo1;
        uint32_t next_x2 = hi0 ^ src[3] ^ k1;
        uint32_t next_x3 = lo0;

        // Next round, and return value here too
        src[0] = next_x0;
        src[1] = next_x1;
        src[2] = next_x2;
        src[3] = next_x3;

        // 4. Update key for the next round
        if (round < 9) {
            k0 += PHILOX_W0;
            k1 += PHILOX_W1;
        }
    }
}

/**
 * @brief philox main function - needs to be a single function to setup the philox and generate, in
 * order for the __shared__ memory to be accessible inside random generator calls
 *
 * in input, counter[0] != 0 means setup. Then
 *   counter[1]: low bits of siteindex
 *   counter[2]: high bits of siteindex
 *   counter[3]: rng "outer" loop (incremented in each onsites())
 *   No output
 * if counter[0] == 0, then in generator phase. In output all counter[i]'s contain random uint32_t
 * values.
 *
 * @internal
 */


#pragma hila loop_function
inline void philox_generate(uint32_t (&counter)[4]) {

#ifndef HILAPP

    // hilapp does not like the function, exclude contents from analysis

#ifdef _GPU_DEVICE_COMPILE_
    // This branch below is used only inside onsites() (kernels) in GPU code

    // Store siteindex and loop index in __shared__ memory for the duration of the kernel,
    // increment local loopcount as random numbers per thread are being used inside the kernel.
    __shared__ uint32_t philox_index0[N_threads];
    __shared__ uint32_t philox_index1[N_threads];
    __shared__ uint32_t philox_loopcount[N_threads];
    __shared__ uint32_t philox_global_loop;

    if (counter[0] == 0) {
        // generate here the random number
        counter[0] = philox_index0[threadIdx.x];
        counter[1] = philox_index1[threadIdx.x];
        counter[2] = philox_global_loop;
        counter[3] = philox_loopcount[threadIdx.x]++; // increment local loop count

        philox4x32_10(counter);

    } else if (counter[0] == 1) {
        // initializing branch
        philox_loopcount[threadIdx.x] = 0;
        philox_index0[threadIdx.x] = counter[1];
        philox_index1[threadIdx.x] = counter[2];
    
    } else {
        // separate initialization for global loop index
        // Needs to be separate because __syncthreads() may cause 
        // trouble otherwise!
        if (threadIdx.x == 0) 
            philox_global_loop = counter[3];
        __syncthreads();
    }

#else

    // This is for non-gpu code and for host rng in GPU code

    static uint32_t philox_index0;
    static uint32_t philox_index1;
    static uint32_t philox_loopcount;
    static uint32_t philox_global_loop;

    if (counter[0] == 0) {
        // generate here the random number
        counter[0] = philox_index0;
        counter[1] = philox_index1;
        counter[2] = philox_global_loop;
        counter[3] = philox_loopcount++;

        philox4x32_10(counter);

    } else {
        // initializing branch
        philox_loopcount = 0;
        philox_index0 = counter[1];
        philox_index1 = counter[2];
        philox_global_loop = counter[3];
    }
#endif

#else
    // this is HILAPP - insert call to philox to have chain of functions
    uint32_t cntr[4];
    philox4x32_10(cntr);
#endif

    // End of function philox_generate
}

/**
 * @brief philox_setup must be called at the beginning of each onsites-block if random numbers
 * are used there - note, on GPU loops we have 2 setup functions!
 */
#pragma hila loop_function
inline void philox_setup_index(uint64_t index) {
    uint32_t counter[4];
    counter[0] = 1;
    counter[1] = static_cast<uint32_t>(index);
    counter[2] = static_cast<uint32_t>(index >> 32);
    philox_generate(counter);
}

/**
 * @brief philox_setup_loopcount must be called at the very beginning of each onsites-block if random numbers
 * are used there
 */
#pragma hila loop_function
inline void philox_setup_loopcount(uint32_t loop) {
    uint32_t counter[4];
    counter[0] = 2;
    counter[3] = loop;
    philox_generate(counter);
}

inline void philox_setup(uint32_t loop, uint64_t index) {
    uint32_t counter[4];
    counter[0] = 1;
    counter[1] = static_cast<uint32_t>(index);
    counter[2] = static_cast<uint32_t>(index >> 32);
    counter[3] = loop;
    philox_generate(counter);
}

/// set up philox for host use, defined in random.cpp
/// to be called after onsites(), inserted by hilapp
void philox_after_onsites();


void seed_random(uint64_t seed, bool device_init = true);
/**
 *@brief Check if RNG is seeded already
 */
bool is_rng_seeded();

/// Empty stub is fine for this
inline void free_device_rng() {}
inline bool is_device_rng_on() {
    return is_rng_seeded();
}


//////////////////////////////////////////////////////////////////////////////////////////
/// generating functions
//////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief return double precision random number in the interval [0,1)
 */

#pragma hila loop_function contains_rng
inline double random() {
    uint32_t counter[4];
    counter[0] = 0;
    philox_generate(counter);
    // make one 53 bit quantity for the double.
    uint64_t x = (static_cast<uint64_t>(counter[0]) << 32) | counter[1];
    return (x >> 11) * 0x1.0p-53; // 2^-53
}

/**
 * @brief return 2 double precision random numbers in the interval [0,1)
 */

#pragma hila loop_function contains_rng
inline double random2(out_only double &d2) {
    uint32_t counter[4];
    counter[0] = 0;
    philox_generate(counter);
    // make one 53 bit quantity for the double.
    uint64_t x = (static_cast<uint64_t>(counter[0]) << 32) | counter[1];
    d2 = (x >> 11) * 0x1.0p-53; // 2^-53
    x = (static_cast<uint64_t>(counter[2]) << 32) | counter[3];
    return (x >> 11) * 0x1.0p-53; // 2^-53
}

/**
 * @brief return 64-bit unsigned int (random bits)
 */

#pragma hila loop_function contains_rng
inline double random_uint64() {
    uint32_t counter[4];
    counter[0] = 0;
    philox_generate(counter);
    return (static_cast<uint64_t>(counter[0]) << 32) | counter[1];
}


/**
 * @brief Gaussian random generation routine
 */
#pragma hila contains_rng loop_function
double gaussrand();

/**
 *@brief `hila::gaussrand2` returns 2 Gaussian distributed random numbers with variance \f$1.0\f$.
 *@details Useful because Box-Muller algorithm computes 2 values at the same time.
 */
#pragma hila contains_rng loop_function
double gaussrand2(double &out2);

/**
 *@brief Check if RNG is initialized, do what the name says.
 */
void check_that_rng_is_initialized();

/**
 *@brief Template function `const T & hila::random(T & var)`
 *       sets the argument to a random value, and return a constant reference to it.
 *@details For example
 *\code{.cpp}
 *   Complex<double> c;
 *   auto n = hila::random(c).abs();
 *\endcode
 *sets the variable `c` to complex random value and calculates its absolute value.
 *`c.real()` and `c.imag()` will be \f$\in [0,1)\f$.
 *
 *For hila classes relies on the existence of method `T::random()` (i.e. `var.random()`),
 *this function typically sets the argument real numbers to interval \f$[0,1)\f$ if `type T` is
 *arithmatic. if T is more commplicated classes such as `SU<N,T>`-matrix, this function sets the
 *argument to valid random `SU<N,T>`.
 *
 * Advantage of this function over class function `T::random()` is that the argument can be
 * of elementary arithmetic type.
 */
template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
T random(out_only T &val) {
    val = hila::random();
    return val;
}

template <typename T, std::enable_if_t<!std::is_arithmetic<T>::value, int> = 0>
T &random(out_only T &val) {
    val.random();
    return val;
}

/**
 *@brief Template function `T hila::random<T>()` without argument.
 *@details This is used to generate random value for `type T` without defined variable.
 *Example:
 *\code{.cpp}
 *    auto n = hila::random<Complex<double>>().abs();
 *\endcode
 *  calculates the norm of a random complex value.
 *`hila::random<double>()` is functionally equivalent to `hila::random()`
 */
template <typename T>
T random() {
    T val;
    hila::random(val);
    return val;
}

/**
 * @brief Template function
 *        `const T & hila::gaussian_random(T & variable,double width=1)`
 * @details Sets the argument to a gaussian random value, and return a constant reference to it.
 * Optional second argument width sets the \f$variance=width^{2}\f$ (\f$default==1\f$)
 *
 * For example:
 * \code {.cpp}
 * Complex<double> c;
 * auto n = sqr(hila::gaussian_random(c));
 * \endcode
 * sets the variable `c` to complex gaussian random value and stores its square in `n`.
 *
 * This function is for hila classes relies on the existence of method `T::gaussian_random()`.
 * The advantage for this function over class function `T::random()` is that the argument can be
 * of elementary arithmetic type.
 * @return T by reference
 */
template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
T gaussian_random(out_only T &val, double w = 1.0) {
    val = hila::gaussrand() * w;
    return val;
}


template <typename T, std::enable_if_t<!std::is_arithmetic<T>::value, int> = 0>
T &gaussian_random(out_only T &val, double w = 1.0) {
    val.gaussian_random(w);
    return val;
}


/**
 *@brief Template function
 *       `T hila::gaussian_random<T>()`,generates gaussian random value of `type T`, with variance
 *\f$1\f$.
 *@details For example,
 *\code{.cpp}
 *  auto n = hila::gaussian_random<Complex<double>>().abs();
 *\endcode
 *calculates the norm of a gaussian random complex value.
 *
 *@note there is no width/variance parameter, because of danger of confusion
 *with above `hila::gaussian_random(value)`
 */
template <typename T>
T gaussian_random() {
    T val;
    hila::gaussian_random(val);
    return val;
}


// just trivial stubs here

// inline double philox_double() {
//     return 0.0;
// }
// inline double philox_double2(double &d2) {
//     return 0.0;
// }
// inline uint64_t philox_uint64() {
//     return 0;
// }
// inline void philox_setup(uint32_t loop, uint64_t index) {}

} // namespace hila


#endif // HILA_PHILOX_H_
