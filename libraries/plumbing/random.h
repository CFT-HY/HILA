#ifndef HILA_RANDOM_H_
#define HILA_RANDOM_H_

namespace hila {

constexpr bool device_rng_off = false;
constexpr bool device_rng_on = true;


/**
 *@brief Seed random generators with 64-bit unsigned value. 
 *       On MPI shuffles the seed so that different MPI ranks are seeded with different values.
 */  
void seed_random(uint64_t seed, bool device_rng = true);

/**
 *@brief Check if RNG is seeded already
 */
bool is_rng_seeded();

/**
 *@brief Initialize host (CPU) random number generator separately, done implicitly by `seed_random()`
 */
void initialize_host_rng(uint64_t seed);

/**
 *@brief Initialize device random number generator on GPUs, if application run on GPU platform. 
 *       No effect on other archs.
 */
void initialize_device_rng(uint64_t seed);

/**
 *@brief Free GPU RNG state, does nothing on non-GPU archs.
 */  
void free_device_rng();

/**
 *@brief Check if the RNG on GPU is allocated and ready to use. 
 */  
bool is_device_rng_on();


/////////////////////////////////////////////////////////////////////////////////////////////////

// It is important that the random number generators random(), gaussrand() and gaussrand2()
// are marked as "loop function" and "contains rng", because hilapp does not have a view
// inside them from all compilation units - although hila::random() has special treatment


/**
 * @brief Real valued uniform random number generator.
 * @details Returns an uniform double precision random number in interval \f$[0,1)\f$.  
 * This function can be called from outside or inside site loops (on GPU if the device rng is initialized).
 *
 * Uses
 * [std::uniform_real_distribution](https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution)
 *
 * @return double
 */
#pragma hila contains_rng loop_function
double random();

// alias for host (CPU) rng, used in GPU code generation. Must not be used inside onsites() {}. This
// routine is not normally needed in user code, instead use standard hila::random().
double host_random();


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
 *this function typically sets the argument real numbers to interval \f$[0,1)\f$ if `type T` is arithmatic.
 *if T is more commplicated classes such as `SU<N,T>`-matrix, this function sets the argument to
 * valid random `SU<N,T>`.
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
 *       `T hila::gaussian_random<T>()`,generates gaussian random value of `type T`, with variance \f$1\f$.
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


} // namespace hila


#endif
