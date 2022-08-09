#ifndef PARAMS_H_
#define PARAMS_H_

///////////////////////////////////////////////////////////////////////////
///  This file contains #defined constants
///  These can be overruled in app Makefile, with 
///  APP_OPTS := -DPARAMETER=value
///  On switches which are by default on, "-DPARAMETER=0" undefines them
///////////////////////////////////////////////////////////////////////////


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
// locality
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

// If SLOW_GPU_REDUCTION is defined, use slow but memory stingy
// reduction.  Probably should not be used.
// #define SLOW_GPU_REDUCTION

#endif  // CUDA || HIP

///////////////////////////////////////////////////////////////////////////
// Special defines for CUDA target


#if defined(CUDA)

// General # of threads
#ifndef N_threads
#define N_threads 256
#endif

// # of threads in reduction - as large as HW supports
#ifndef N_GPU_REDUCE_THREADS
#define N_GPU_REDUCE_THREADS 512
#endif

// which reduction kernel to use - see gpu_reduction.h
#ifndef GPU_REDUCE_KERNEL
#define GPU_REDUCE_KERNEL 6
#endif

// How many fft's in parallel - large value faster, small less memory.
#ifndef GPUFFT_BATCH_SIZE
#define GPUFFT_BATCH_SIZE 256
#endif

// Undef cuda-aware mpi at makefile with -DGPU_AWARE_MPI=0
#ifndef GPU_AWARE_MPI
#define GPU_AWARE_MPI 1
#elif GPU_AWARE_MPI == 0
#undef GPU_AWARE_MPI
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

// # of threads in reduction - as large as HW supports
#ifndef N_GPU_REDUCE_THREADS
#define N_GPU_REDUCE_THREADS 512
#endif

// which reduction kernel to use - see gpu_reduction.h
#ifndef GPU_REDUCE_KERNEL
#define GPU_REDUCE_KERNEL 6
#endif

// How many fft's in parallel - large value faster, small less memory.
#ifndef GPUFFT_BATCH_SIZE
#define GPUFFT_BATCH_SIZE 256
#endif


// End of GPU defines
#endif   // HIP

#endif
