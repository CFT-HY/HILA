#ifndef PARAMS_H_
#define PARAMS_H_

///////////////////////////////////////////////////////////////////////////
///  This file contains #defined constants
///  These can be overruled in app Makefile, with 
///  APP_OPTS := -DPARAMETER=value
///  On switches which are by default on, "-DPARAMETER=0" undefines them
///////////////////////////////////////////////////////////////////////////


// Dimensionality
#ifndef NDIM
#define NDIM 4
#endif

// output file name
#ifndef DEFAULT_OUTPUT_NAME
#define DEFAULT_OUTPUT_NAME "output"
#endif

// EVEN_SITES_FIRST is the default
#ifndef EVEN_SITES_FIRST
#define EVEN_SITES_FIRST
#elif EVEN_SITES_FIRST == 0
#undef EVEN_SITES_FIRST
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

// Undef cuda-aware mpi at makefile with -DCUDA_AWARE_MPI=0
#ifndef CUDA_AWARE_MPI
#define CUDA_AWARE_MPI 1
#elif CUDA_AWARE_MPI == 0
#undef CUDA_AWARE_MPI
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