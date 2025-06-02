#ifndef PARAMS_H_
#define PARAMS_H_
/**
 * @file params.h
 * @brief This file contains #defined constants
 * @details
 * These can be overruled in application Makefile, with APP_OPTS := -DPARAMETER=value.
 *
 * There are two types of #define variables, True/False switches or parameter variables.
 *
 * True/False statements can be set with either 0 (False) or 1 (True) as -DPARAMETER=0.
 *
 * Parameter variables are set similary with -DPARAMETER=var where var is the chosen variable
 */

#ifdef RELEASE
/**
 * @brief Turn off asserts which are on by default.
 * @details By defining either RELEASE or NDEBUG (No debug) asserts will be turned off.
 * Static asserts naturally remain active
 */
#ifndef NDEBUG
#define NDEBUG
#endif
#endif

#ifndef NDIM
/**
 * @brief HILA system dimensionality
 * @details Set's HILA dimensionality for which 4 is default. Options are 2,3,4
 */
#define NDIM 4
#endif

#ifndef DEFAULT_OUTPUT_NAME
/**
 * @def DEFAULT_OUTPUT_NAME
 * @brief Default output file name
 */
#define DEFAULT_OUTPUT_NAME "output"
#endif

#ifndef EVEN_SITES_FIRST
/**
 * @brief  EVEN_SITES_FIRST is default. To traverse sites on natural (parity invariant) order use
 * -DEVEN_SITES_FIRST=0
 */
#define EVEN_SITES_FIRST
#elif EVEN_SITES_FIRST == 0
#undef EVEN_SITES_FIRST
#endif

/// NODE_LAYOUT_TRIVIAL or NODE_LAYOUT_BLOCK determine how MPI ranks are laid out on logical
/// lattice.  TRIVIAL lays out the lattice on logical order where x-direction runs fastest etc.
/// if NODE_LAYOUT_BLOCK is defined, NODE_LAYOUT_BLOCK consecutive MPI ranks are laid out so that
/// these form a compact "block" of ranks logically close togeter.
/// Define NODE_LAYOUT_BLOCK to be the number of
/// MPI processes within one compute node - tries to maximize the use of fast local communications.
/// Either one of these must be defined.

#ifndef NODE_LAYOUT_TRIVIAL
#ifndef NODE_LAYOUT_BLOCK
#define NODE_LAYOUT_BLOCK 4
#endif
#endif

/// WRITE_BUFFER SIZE
/// Size of the write buffer in field writes, in bytes
/// Larger buffer -> less MPI calls in writing, but more memory
#ifndef WRITE_BUFFER_SIZE
#define WRITE_BUFFER_SIZE 2000000
#endif


// boundary conditions are "off" by default -- no need to do anything here
// #ifndef SPECIAL_BOUNDARY_CONDITIONS

///////////////////////////////////////////////////////////////////////////
// Special defines for GPU targets
#if defined(CUDA) || defined(HIP)

/// Use gpu memory pool by default
/// turn off by using -DGPU_MEMORY_POOL=0 in Makefile
#ifndef GPU_MEMORY_POOL
#define GPU_MEMORY_POOL
#elif GPU_MEMORY_POOL == 0
#undef GPU_MEMORY_POOL
#endif

/// GPU_AWARE_MPI
/// By default GPU aware MPI is on. Turn it off in Makefile with -DGPU_AWARE_MPI=0
#ifndef GPU_AWARE_MPI
#define GPU_AWARE_MPI 1
#elif GPU_AWARE_MPI == 0
#undef GPU_AWARE_MPI
#endif

/// GPU_RNG_THREAD_BLOCKS
/// Number of thread blocks (of N_threads threads) to use in onsites()-loops containing random
/// numbers. GPU_RNG_THREAD_BLOCKS=0 or undefined means use one RNG on each lattice site, and the
/// thread block number is not restricted.  RNG takes about 48 B/generator (with XORWOW). When
/// GPU_RNG_THREAD_BLOCKS > 0 only (N_threads * GPU_RNG_THREAD_BLOCKS) generators are in use, which
/// reduces the memory footprint substantially (and bandwidth demand) Too small number slows down
/// onsites()-loops containing RNGs, because less threads are active. Example:
///     Field<Vector<4,double>> vfield;
///     onsites(ALL) {
///        vfield[X].gaussian_random();      // there's RNG here, so this onsites() is handled by
///                                          // GPU_RNG_THREAD_BLOCKS thread blocks
///     }
/// GPU_RNG_THREAD_BLOCKS<0 disables GPU random numbers entirely, and loops like above will crash if
/// executed. hilapp will emit a warning, but program is compiled
///
///  Default: 32 seems to be OK compromise. Can be set to 0 if memory is not a problem.

#ifndef GPU_RNG_THREAD_BLOCKS
#define GPU_RNG_THREAD_BLOCKS 32
#endif

/// GPU_VECTOR_REDUCTION_THREAD_BLOCKS
/// # of thread blocks (of N_threads threads) used in ReductionVector (weighted histogram) ops.
/// A value > 0 for GPU_VECTOR_REDUCTION_THREAD_BLOCKS means that the onsites-loop where the
/// reduction is done is handled by GPU_VECTOR_REDUCTION_THREAD_BLOCKS thread blocks of N_threads
/// threads.  Each thread handles its own histogram, thus there are
/// (GPU_VECTOR_REDUCTION_THREAD_BLOCKS*N_threads) working copies of the histogram which are then
/// combined. Too small value slows the loop where this happens computation, too large uses
/// (temporarily) more memory. Example:
///      ReductionVector<double> rv(100);
///      Field<int> index;
///      ... (set index to values 0 .. 99)
///      onsites(ALL) {
///           rv[index[X]] += ..
///           ..
///      }
///
/// GPU_VECTOR_REDUCTION_THREAD_BLOCKS = 0 or undefined means that the thread block number is not
/// restricted and only a single histogram is used with atomic operations (atomicAdd).  This
/// can be slower, but the performance is GPU hardware/driver dependent.  In some
/// cases GPU_VECTOR_REDUCTION_THREAD_BLOCKS = 0 turns out to be faster.
///
/// Default: 32 is currently OK compromise (32 thread blocks)

#ifndef GPU_VECTOR_REDUCTION_THREAD_BLOCKS
#define GPU_VECTOR_REDUCTION_THREAD_BLOCKS 32
#endif


/// GPU_VECTOR_REDUCTION_THREADS defines max threads per block for block reduction
/// ReductionVector uses cub::blockreduce with this many threads per block. Too large
/// value can exhaust the resources on GPUs, which will give runtime error. 

#ifndef GPU_BLOCK_REDUCTION_THREADS
#define GPU_BLOCK_REDUCTION_THREADS 128
#endif


/// GPU_VECTOR_REDUCTION_SIZE_THRESHOLD is an optimization parameter. If reduction size
/// is large than threshold, we use "single pass" kernel; if smaller, hierarchial
/// kernel launch. Both use the same amount of memory. The algorithms become
/// correspondingly better at extreme ends, but around 500-1000 there is
/// a slow changeover, depending on computing hardware. 
/// NOT USED IN PRESENT ReductionVector implementation

// #ifndef GPU_VECTOR_REDUCTION_SIZE_THRESHOLD
// #define GPU_VECTOR_REDUCTION_SIZE_THRESHOLD 700
// #endif


/// GPUFFT_BATCH_SIZE
/// How many complex fft's in parallel - large value can be faster, small uses less memory.
/// Performance is reduced if the value is too small, but levels to a ~constant
/// when sufficiently large.
#ifndef GPUFFT_BATCH_SIZE
#define GPUFFT_BATCH_SIZE 256
#endif

/** @brief GPU_SYNCHRONIZE_TIMERS : if set and !=0 synchronize GPU on timer calls, in order to
 * obtain meaningful timer values
 *
 * @details Because GPU kernel launch is asynchronous process, the timers by default may not measure
 * the actual time used in GPU kernel execution. Defining GPU_SYNCHRONIZE_TIMERS inserts GPU
 * synchronization calls to timers.  This is off by default, because this may slow down GPU code.
 * Turn on in order to measure more accurately the time spent in different parts of the code.
 */

#ifdef GPU_SYNCHRONIZE_TIMERS
#if GPU_SYNCHRONIZE_TIMERS == 0
#undef GPU_SYNCHRNONIZE_TIMERS
#endif
#endif


/** @brief GPU_GLOBAL_ARG_MAX_SIZE : in some __global__functions gives the max size of variable
 * passed directly as an argument of the function call. Larger value sizes are passed with
 * gpuMemcopy() and a pointer. CUDA < 12.1 limits the total parameter size to 4K, >= 12.1 it is 32K.
 * We set the default to 2K. in HIP/rocm I have not found the size. Passing as an arg is faster, but
 * size limit is relevant only for big "matrices"
 */

#ifndef GPU_GLOBAL_ARG_MAX_SIZE
#define GPU_GLOBAL_ARG_MAX_SIZE 2048
#endif

#endif // CUDA || HIP

///////////////////////////////////////////////////////////////////////////
// Special defines for CUDA target

#if defined(CUDA)

/// General number of threads in a thread block
#ifndef N_threads
#define N_threads 256
#endif


#ifndef GPU_MEMORY_POOL

// CUDA_MALLOC_ASYNC
#ifndef CUDA_MALLOC_ASYNC
// Use async malloc only if version is large enough
// NOTE: does not seem to work with OpenMPI, disable
#if 0 && CUDART_VERSION >= 11020
#define CUDA_MALLOC_ASYNC
#endif

#elif CUDA_MALLOC_ASYNC == 0
#undef CUDA_MALLOC_ASYNC
#endif

#endif // if not GPU_MEMORY_POOL

#endif // CUDA

///////////////////////////////////////////////////////////////////////////
// Same for HIP

#if defined(HIP)

// General number of threads in a thread block
#ifndef N_threads
#define N_threads 256
#endif


// End of GPU defines
#endif // HIP

#endif
