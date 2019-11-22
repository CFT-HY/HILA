#ifndef CUDA_INCLUDE_H
#define CUDA_INCLUDE_H


#ifdef __CUDACC__

#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <cub/cub.cuh>

#define N_threads 64

extern __device__ unsigned * d_neighb[NDIRS];

/* Random number generator */
extern curandState * curandstate;
__device__ extern curandState * d_curandstate;
__global__ void seed_random_kernel( curandState * state, unsigned long seed );
void seed_random(unsigned long seed);
__host__ __device__ double hila_random();


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
T cuda_reduce_sum(  T * vector, int N ){
  static bool initialized = false;
  static T * d_sum;
  static void *d_temp_storage;
  static size_t temp_storage_size;

  T sum;

  if( !initialized ) {
    cudaMalloc( (void **)&d_sum, sizeof(T) );
    d_temp_storage = NULL;
    temp_storage_size = 0;
    cub::DeviceReduce::Sum(d_temp_storage, temp_storage_size, vector, d_sum, N);
    initialized = true;
  }
  
  // Allocate temporary storage
  cudaMalloc(&d_temp_storage, temp_storage_size);

  // Run sum-reduction
  cub::DeviceReduce::Sum(d_temp_storage, temp_storage_size, vector, d_sum, N);

  cudaMemcpy( &sum, d_sum, sizeof(T), cudaMemcpyDeviceToHost );

  return sum;
}


// A simple hand-written reduction that does not require a library
//
//template<typename T>
//__global__ void cuda_reduce_sum_kernel( T * vector, int vector_size, int new_size, int elems)
//{
//  int Index = threadIdx.x + blockIdx.x * blockDim.x;
//  if( Index < new_size ){
//    for(int i=1; i<elems; i++){
//      int ind = Index + i*new_size;
//      if( ind < vector_size ){
//        vector[Index] += vector[ ind ];
//      }
//    }
//  }
//}
//
//
//template<typename T>
//T cuda_reduce_sum(  T * vector, int N ){
//  const int reduce_step = 32;
//  T sum;
//  T * host_vector = (T *)malloc(N*sizeof(T));
//  int vector_size = N;
//
//  while( vector_size > reduce_step ){
//    int new_size = vector_size/reduce_step + 1;
//    int blocks = new_size/N_threads + 1;
//    cuda_reduce_sum_kernel<<< blocks, N_threads >>>( vector, vector_size, new_size, reduce_step );
//    vector_size = new_size;
//    cudaDeviceSynchronize();
//  }
//  
//
//  check_cuda_error("cuda_reduce_sum kernel");
//  cudaMemcpy( host_vector, vector, vector_size*sizeof(T), cudaMemcpyDeviceToHost );
//  check_cuda_error("Memcpy in reduction");
//
//  for(int i=0; i<vector_size; i++){
//    sum += host_vector[i];
//  }
//
//  free(host_vector);
//  
//  return sum;
//}


inline void synchronize_threads(){
  cudaDeviceSynchronize();
}

void initialize_cuda(int rank);

#else

// define cuda functions in order to avoid compilation errors
// in transformer
#define hila_random() 0.5
#define seed_random(seed) while(0)

#define cudaMalloc(a,b) 0
#define cudaFree(a) while(0)
#define check_cuda_error(a) while(0)

inline void synchronize_threads(){}
void initialize_cuda(int rank){};

#endif

#endif
