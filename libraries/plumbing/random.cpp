#include <cmath>
#include "hila.h"
#include "plumbing/random.h"

/////////////////////////////////////////////////////////////////////////

#include <random>

//#ifndef OPENMP

// static variable which holds the random state
// Use 64-bit mersenne twister
static std::mt19937_64 mersenne_twister_gen;

// random numbers are in interval [0,1)
static std::uniform_real_distribution<double> real_rnd_dist(0.0,1.0);

//#endif


// In GPU code hila::random() defined in hila_gpu.cpp
#if !defined(CUDA) && !defined(HIP)

double hila::random() {
    return real_rnd_dist( mersenne_twister_gen );
}

#endif

// Generate random number in non-kernel (non-loop) code.  Not meant to
// be used in "user code"
double hila::host_random() {
    return real_rnd_dist( mersenne_twister_gen );
}

/////////////////////////////////////////////////////////////////////////

/// Seed random number generators
/// Seed is shuffled so that different nodes
/// get different rng seeds.  If seed == 0,
/// generate seed using the time() -function.

void hila::seed_random(uint64_t seed) {

#ifndef SITERAND

    uint64_t n = hila::myrank();
    if (hila::partitions.number() > 1)
        n += hila::partitions.mylattice() * hila::number_of_nodes();

    if (seed == 0) {
        // get seed from time
        if (hila::myrank() == 0) {
            struct timespec tp;

            clock_gettime(CLOCK_MONOTONIC, &tp);
            seed = tp.tv_sec;
            seed = (seed << 30) ^ tp.tv_nsec;
            output0 << "Random seed from time: " << seed << '\n';
        }
        hila::broadcast(seed);
    }
    // do node shuffling for seed
    // do it in a manner makes it difficult to give the same seed by mistake
    // and also avoids giving the same seed for 2 nodes
    // n=0 remains unchanged

    seed = seed ^ (n ^ ((7*n) << 25));

    output0 << "Using node random numbers, seed for node 0: " << seed << std::endl;

// #if !defined(OPENMP)
    mersenne_twister_gen.seed( seed );
    // warm it up
    for (int i = 0; i < 9000; i++)
        mersenne_twister_gen();
// #endif


#if defined(CUDA) || defined(HIP)
    // we can use the same seed, the generator is different
    hila::seed_device_rng(seed);
#endif

    // taus_initialize();

#else

    // TODO: clean SITERAND!
    // Now SITERAND is defined
    // This is usually used only for occasional benchmarking, where identical output
    // independent of the node number is desired

    output0 << "*** SITERAND is in use!\n";

    random_seed_arr =
        (unsigned short(*)[3])memalloc(3 * node.sites * sizeof(unsigned short));
    forallsites(i) {
        random_seed_arr[i][0] = (unsigned short)(seed + site[i].index);
        random_seed_arr[i][1] = (unsigned short)(seed + 2 * site[i].index);
        random_seed_arr[i][2] = (unsigned short)(seed + 3 * site[i].index);
    }

    random_seed_ptr = random_seed_arr[0];

#endif
}




///////////////////////////////////////////////////////////////////
/// gaussian rng generation routines
/// By default these give random numbers with variance 1.0, i.e.
/// probability distribution is
///           exp( -x*x/2 ), so < x^2 > = 1
/// If you want random numbers with variance sigma, multiply the
/// result by  sqrt(sigma):
///       sqrt(sigma) * gaussrand();
////////////////////////////////////////////////////////////////////

#define VARIANCE 1.0

double hila::gaussrand2(double &out2) {

    double phi, urnd, r;
    phi = 2.0 * M_PI * hila::random();

    // this should not really trigger
    do {
        urnd = 1.0 - hila::random();
    } while (urnd == 0.0);

    r = sqrt(-::log(urnd) * (2.0 * VARIANCE));
    out2 = r * cos(phi);
    return r * sin(phi);
}

#if !defined(CUDA) && !defined(HIP)

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

//#pragma hila loop function contains rng
double hila::gaussrand() {
    double second;
    return hila::gaussrand2(second);
}

#endif
