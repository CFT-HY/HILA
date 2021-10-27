#ifndef PARAMS_H_
#define PARAMS_H_

///////////////////////////////////////////////////////////////////////////
//  This file contains #defined constants
//  These can be overruled in Makefile, with "-DPARAM=value"
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

///////////////////////////////////////////////////////////////////////////
// Same for HIP

#elif defined(HIP)

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
#endif     

#endif
