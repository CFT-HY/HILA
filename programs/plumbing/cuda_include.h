#ifndef CUDA_INCLUDE_H
#define CUDA_INCLUDE_H


#ifdef __CUDACC__

#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>

#define N_threads 128

#define loop_callable __host__ __device__

/* Random number generator */
curandState * curandstate;
__device__ curandState * d_curandstate;

/* Set seed */
__global__ void seed_random_kernel( curandState * state, unsigned long seed )
{
  int id = threadIdx.x;
  d_curandstate = state;
  curand_init ( seed, id, 0, &d_curandstate[id] );
}

static inline void seed_random(unsigned long seed){
  cudaMalloc( &curandstate, N_threads*sizeof( curandState ) );
  seed_random_kernel<<< 1, N_threads >>>( curandstate, seed );
  seed_mersenne(seed+N_threads);
}

loop_callable double hila_random(){
  #ifdef __CUDA_ARCH__
  return curand_uniform( &d_curandstate[threadIdx.x] );
  #else
  return mersenne();
  #endif
}


#else

#define loop_callable 

// define cuda functions in order to avoid compilation errors
// in transformer
#define hila_random() 0.5
#define seed_random(seed) while(0)
#define cudaMalloc(a,b) while(0)
#define cudaFree(a) while(0)


#endif

#endif
