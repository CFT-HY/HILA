#ifndef CUDA_INCLUDE_H
#define CUDA_INCLUDE_H


#ifdef __CUDACC__

#include <cuda_runtime.h>
#define random() curand_uniform_double(&state);

#define loop_callable __host__ __device__

#else

#define loop_callable 

// define cuda functions in order to avoid compilation errors
// in transformer
#define random() 0.5
#define cudaMalloc(a,b) while(0)
#define cudaFree(a) while(0)


#endif

#endif
