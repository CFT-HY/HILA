#ifndef RANDOM_H_
#define RANDOM_H_

namespace hila {

constexpr bool device_rng_off = false;
constexpr bool device_rng_on = true;

/// Seed random generators with 64-bit unsigned value. On MPI shuffles the seed so that different
/// MPI ranks are seeded with different values.
///
/// The optional 2nd argument indicates whether to initialize the RNG on GPU device:
/// hila::device_rng_on (default) or hila::device_rng_off.  This argument does nothing if not GPU
/// platform.  If hila::device_rng_off is used, onsites() -loops cannot contain random number calls
/// (Runtime error will be flagged and program exits).

void seed_random(uint64_t seed, bool device_rng = true);

/// @brief Check if RNG is seeded already
bool is_rng_seeded();

///
/// Initialize host (CPU) random number generator separately, done implicitly by seed_random()
/// Parameter random number seed.  On MPI shuffles different seed values for all MPI ranks.

void initialize_host_rng(uint64_t seed);

///
/// Initialize device random number generator on GPUs, if on GPU platform.  No effect on other
/// archs. On MPI shuffles the seed for different MPI ranks.  Called by seed_random() unless its 2nd
/// argument is hila::device_rng_off.  This can reinitialize device RNG free'd by free_device_rng().

void initialize_device_rng(uint64_t seed);

///
/// Free GPU RNG state, hila::random() does not work inside onsites() after this
/// (unless seeded again using initialize_device_rng()).  Frees the memory RNG takes on the device.
/// Does nothing on non-GPU archs.

void free_device_rng();

///
/// Check if the RNG on GPU is allocated and ready to use.  Returns true on non-GPU archs.

bool is_device_rng_on();


/////////////////////////////////////////////////////////////////////////////////////////////////

// It is important that the random number generators random(), gaussrand() and gaussrand2()
// are marked as "loop function" and "contains rng", because hilapp does not have a view
// inside them from all compilation units - although hila::random() has special treatment


/**
 * @brief Real valued uniform random number generator
 * @details Returns a uniform double precision random number in interval [0,1).  Can be
 * called from outside or inside site loops (on GPU if the device rng is initialized).
 *
 *  Uses
 * [std::uniform_real_distribution](https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution)
 *
 * @return double
 */
#pragma hila contains_rng loop_function
double random();

// alias for host (CPU) rng, used in GPU code generation. Must not be used inside onsites() {}. This
// routine is not normally needed in user code, instead use standard hila::random().
double host_random();

///
/// hila::gaussrand() returns a Gaussian distributed random number with variance 1.0, i.e.
/// the probability distribution is
///           exp( -x*x/2 ), so < x^2 > = 1
/// If you want random numbers with variance sigma, multiply the
/// result by  sqrt(sigma):
///       sqrt(sigma) * gaussrand();

#pragma hila contains_rng loop_function
double gaussrand();

/// hila::gaussrand2 returns 2 Gaussian distributed random numbers with variance 1.0.
/// Useful because Box-Muller algorithm computes 2 values at the same time.

#pragma hila contains_rng loop_function
double gaussrand2(double &out2);

///
/// do what the name says - program quits with error message if RNG is not initialized and on GPU
/// archs the device RNG is not initialized.

void check_that_rng_is_initialized();


///
/// Template function
///     const T & hila::random(T & var)
/// Sets the argument to a random value, and return a constant reference to it.
/// Example:
///     Complex<double> c;
///     auto n = hila::random(c).abs();
/// sets the variable c to complex random value and calculates its absolute value.
/// c.real() and c.imag() will be \in [0,1)
///
/// For hila classes relies on the existence of method T::random() (i.e. var.random())
///
/// Typically sets the argument real numbers to interval [0,1), but not
/// always: for example, if T is SU<N,T> -matrix sets the argument to
/// valid random SU<N,T>.
///
/// Advantage over class function T::random() is that the argument can be
/// of elementary arithmetic type.

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

///
/// Template function
///     T hila::random<T>();
/// without argument.  This is used to generate random value for T without defined variable.
/// Example:
///     auto n = hila::random<Complex<double>>().abs();
/// calculates the norm of a random complex value.
///
/// hila::random<double>() is functionally equivalent to hila::random()

template <typename T>
T random() {
    T val;
    hila::random(val);
    return val;
}


///
/// Template function
///     const T & hila::gaussian_random(T & variable,double width=1)
/// Sets the argument to a gaussian random value, and return a constant reference to it.
/// Optional second argument width sets the variance=width^2 (default=1)
///
/// Example:
///     Complex<double> c;
///     auto n = sqr(hila::gaussian_random(c));
/// sets the variable c to complex gaussian random value and stores its square in n.
///
/// For hila classes relies on the existence of method T::gaussian_random().
///
/// Advantage over class function T::random() is that the argument can be
/// of elementary arithmetic type.

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

///
/// Template function
///     T hila::gaussian_random<T>();
/// generates gaussian random value of type T, with variance 1.
/// Example:
///     auto n = hila::gaussian_random<Complex<double>>().abs();
/// calculates the norm of a gaussian random complex value.
///
/// Note: there is no width/variance parameter, because of danger of confusion
/// with above hila::gaussian_random(value)

template <typename T>
T gaussian_random() {
    T val;
    hila::gaussian_random(val);
    return val;
}


} // namespace hila


#endif