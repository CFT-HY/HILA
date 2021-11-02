#ifndef PARAMS_H_
#define PARAMS_H_

///////////////////////////////////////////////////////////////////////////
//  This file contains #defined constants
//  These can be overruled in Makefile, with "-DPARAM=value"
//  On switches which are by default on, "-DPARAM=0" undefs' them
///////////////////////////////////////////////////////////////////////////


// Dimensionality
#ifndef NDIM
#define NDIM 4
#endif

// output file name
#ifndef DEFAULT_OUTPUT_NAME
#define DEFAULT_OUTPUT_NAME "output"
#endif


///////////////////////////////////////////////////////////////////////////
// Special defines for GPU targets
#if defined(CUDA) || defined(HIP)

// Use gpu memory pool by default 
// set off by using -DGPU_MEMORY_POOL=0 in Makefile
#if defined(GPU_MEMORY_POOL)
#if GPU_MEMORY_POOL == 0
#undef GPU_MEMORY_POOL
#endif
#else 
#define GPU_MEMORY_POOL
#endif

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


#ifndef GPU_MEMORY_POOL

// CUDA_MALLOC_ASYNC 
#if defined(CUDA_MALLOC_ASYNC)
#if CUDA_MALLOC_ASYNC == 0
#undef CUDA_MALLOC_ASYNC
#endif
// Now CUDA_MALLOC_ASYNC is not defined
// Use async malloc only if version is large enough
// NOTE: does not seem to work with OpenMPI, disable
#elif 0 && CUDART_VERSION >= 11020
#define CUDA_MALLOC_ASYNC
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
