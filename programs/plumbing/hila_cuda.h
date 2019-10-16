#ifndef CUDA_INCLUDE_H
#define CUDA_INCLUDE_H


#ifdef __CUDACC__

#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>

#define N_threads 128

#define loop_callable __host__ __device__

extern __device__ unsigned * d_neighb[NDIRS];
void set_neighbour_pointers( unsigned * neighb, int d );

/* Random number generator */
extern curandState * curandstate;
__device__ extern curandState * d_curandstate;

__global__ void seed_random_kernel( curandState * state, unsigned long seed );
void seed_random(unsigned long seed);
loop_callable double hila_random();

static inline void check_cuda_error(std::string message){
  cudaError code = cudaGetLastError();
  if( cudaSuccess != code ){
    std::cout << message << ": "
              << cudaGetErrorString(code) << "\n";
    exit(1);
  }
}

/* Reduction */
template<typename T>
__global__ void cuda_reduce_sum_kernel( T * vector, int N)
{
  for( int i=1; i<N; i++ ) {
    vector[0] += vector[i];
  }
}

template<typename T>
T cuda_reduce_sum(  T * vector, int N ){
  T sum;
  cuda_reduce_sum_kernel<<< 1, 1 >>>( vector, N );
  check_cuda_error("cuda_reduce_sum kernel");
  cudaMemcpy( &sum, vector, sizeof(T), cudaMemcpyDeviceToHost );
  check_cuda_error("Memcpy in reduction");
  return sum;
}

#else

#define loop_callable 

// define cuda functions in order to avoid compilation errors
// in transformer
#define hila_random() 0.5
#define seed_random(seed) while(0)

#define cudaMalloc(a,b) 0
#define cudaFree(a) while(0)
#define check_cuda_error(a) while(0)

#endif

#endif
