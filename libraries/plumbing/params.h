#ifndef PARAMS_H_
#define PARAMS_H_

///////////////////////////////////////////////////////////////////////////
///  This file contains #defined constants
///  These can be overruled in app Makefile, with
///  APP_OPTS := -DPARAMETER=value
///  On switches which are by default on, "-DPARAMETER=0" undefines them
///
///  These can be set on the make command line with
///  make OPTS="-DPARAMETER=value -DPARAMETER2=value2"
///////////////////////////////////////////////////////////////////////////

/// Assertions on by default? Turn them off by defining NDEBUG or RELEASE
#ifdef RELEASE
#ifndef NDEBUG
#define NDEBUG
#endif
#endif

/// Dimensionality
#ifndef NDIM
#define NDIM 4
#endif

/// output file name
#ifndef DEFAULT_OUTPUT_NAME
#define DEFAULT_OUTPUT_NAME "output"
#endif

// EVEN_SITES_FIRST is the default
#ifndef EVEN_SITES_FIRST
#define EVEN_SITES_FIRST
#elif EVEN_SITES_FIRST == 0
#undef EVEN_SITES_FIRST
#endif

// NODE_LAYOUT_TRIVIAL or NODE_LAYOUT_BLOCK must be defined
// Define NODE_LAYOUT_BLOCK to be the number of
// MPI processes within a compute node - tries to maximize
// locality somewhat
#ifndef NODE_LAYOUT_TRIVIAL
#ifndef NODE_LAYOUT_BLOCK
#define NODE_LAYOUT_BLOCK 4
#endif
#endif

// Size of the write buffer in field writes, in bytes
// Larger buffer -> less MPI calls in writing, but more memory
#ifndef WRITE_BUFFER_SIZE
#define WRITE_BUFFER_SIZE 2000000
#endif


// boundary conditions are "off" by default -- no need to do anything here
// #ifndef SPECIAL_BOUNDARY_CONDITIONS

///////////////////////////////////////////////////////////////////////////
// Special defines for GPU targets
#if defined(CUDA) || defined(HIP)

// Use gpu memory pool by default
// set off by using -DGPU_MEMORY_POOL=0 in Makefile
#ifndef GPU_MEMORY_POOL
#define GPU_MEMORY_POOL
#elif GPU_MEMORY_POOL == 0
#undef GPU_MEMORY_POOL
#endif

// Undef cuda/hip -aware mpi at makefile with -DGPU_AWARE_MPI=0
#ifndef GPU_AWARE_MPI
#define GPU_AWARE_MPI 1
#elif GPU_AWARE_MPI == 0
#undef GPU_AWARE_MPI
#endif

// GPU_RNG_THREAD_BLOCKS
// Number of thread blocks (of N_threads threads) to use in onsites()-loops containing random numbers.
// GPU_RNG_THREAD_BLOCKS=0 or undefined means use one RNG on each lattice site, and the thread block
// number is not restricted.  RNG takes about 48 B/generator (with XORWOW).
// When GPU_RNG_THREAD_BLOCKS > 0 only (N_threads * GPU_RNG_THREAD_BLOCKS)
// generators are in use, which reduces the memory footprint substantially (and bandwidth demand)
// Too small number slows down onsites()-loops containing RNGs, because less threads are active.
// Example:
//     Field<Vector<4,double>> vfield;
//     onsites(ALL) {
//        vfield[X].gaussian_random();      // there's RNG here, so this onsites() is handled by
//                                          // GPU_RNG_THREAD_BLOCKS thread blocks
//     }
// GPU_RNG_THREAD_BLOCKS<0 disables GPU random numbers entirely, and loops like above will crash if executed.
// hilapp will emit a warning, but program is compiled

#ifndef GPU_RNG_THREAD_BLOCKS
#define GPU_RNG_THREAD_BLOCKS 32
#endif

// GPU_VECTOR_REDUCTION_THREAD_BLOCKS
// # of thread blocks (of N_threads threads) used in ReductionVector (weighted histogram) ops.
// A value > 0 for GPU_VECTOR_REDUCTION_THREAD_BLOCKS means that the onsites-loop where the
// reduction is done is handled by GPU_VECTOR_REDUCTION_THREAD_BLOCKS thread blocks of N_threads
// threads.  Each thread handles its own histogram, thus there are 
// (GPU_VECTOR_REDUCTION_THREAD_BLOCKS*N_threads) working copies of the histogram which are then
// combined. Too small value slows the loop where this happens computation, too large uses (temporarily)
// more memory.
// Example:
//      ReductionVector<double> rv(100);
//      Field<int> index;
//      ...
//      onsites(ALL) {
//           rv[index[X]] += ..
//           ..
//      }
//
// GPU_VECTOR_REDUCTION_THREAD_BLOCKS = 0 or undefined means that the thread block number is not
// restricted and only a single histogram is used with atomic operations (atomicAdd).  This 
// can slow down tight loops, but in many cases this turns out to be actually faster.
// 
// Default: leave it off.  Otherwise 32 is currently OK compromise (32 thread blocks)

// #define GPU_VECTOR_REDUCTION_THREAD_BLOCKS 32


// GPUFFT_BATCH_SIZE:
// How many complex fft's in parallel - large value faster, small less memory.
// Performance is reduced if the value is too small, but levels to a ~constant
// when sufficiently large.
#ifndef GPUFFT_BATCH_SIZE
#define GPUFFT_BATCH_SIZE 256
#endif


#endif // CUDA || HIP

///////////////////////////////////////////////////////////////////////////
// Special defines for CUDA target

#if defined(CUDA)

// General # of threads
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

// General # of threads
#ifndef N_threads
#define N_threads 256
#endif


// End of GPU defines
#endif // HIP

#endif
