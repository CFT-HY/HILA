#ifndef CUDA_INCLUDE_H
#define CUDA_INCLUDE_H


// On Puhti, use UCX_MEMTYPE_CACHE=n with
// CUDA_AWARE_MPI
#define CUDA_AWARE_MPI


#ifdef __CUDACC__

#include <sstream>
#include<iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
//#include <cub/cub.cuh>

#define N_threads 128

/* Random number generator */
extern curandState * curandstate;
__device__ extern curandState * d_curandstate;
__global__ void seed_random_kernel( curandState * state, unsigned long seed, unsigned int stride );
void seed_random(unsigned long seed);
__host__ __device__ double hila_random();

#define check_cuda_error(msg) cuda_exit_on_error(msg,__FILE__,__LINE__)
#define check_cuda_error_code(code,msg) cuda_exit_on_error(code,msg,__FILE__,__LINE__)
void cuda_exit_on_error(const char *msg, const char *file, int line);
void cuda_exit_on_error(cudaError code, const char *msg, const char *file, int line);


/* Reduction */
/*
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
*/


// A simple hand-written reduction that does not require a library
template<typename T>
__global__ void cuda_reduce_sum_kernel( T * vector, int vector_size, int new_size, int elems)
{
  int Index = threadIdx.x + blockIdx.x * blockDim.x;
  if( Index < new_size ){
    for(int i=1; i<elems; i++){
      int ind = Index + i*new_size;
      vector[Index] += vector[ ind ];
    }
  }
}

template<typename T>
T cuda_reduce_sum(  T * vector, int N ){
  const int reduce_step = 32;
  T sum = 0;
  T * host_vector = (T *)malloc(N*sizeof(T));
  int vector_size = N;

  while( vector_size > reduce_step ){
    // Take the last n elements that are divisible by reduce_step
    int first = vector_size%reduce_step;
    // Calculate the size of the reduced list
    int new_size = (vector_size-first)/reduce_step;
    // Find number of blocks and launch the kernel
    int blocks = (new_size-1)/N_threads + 1;
    cuda_reduce_sum_kernel<<< blocks, N_threads >>>( vector+first, vector_size, new_size, reduce_step );
    check_cuda_error("cuda_reduce_sum kernel");
    // Find the full size of the resulting array
    vector_size = new_size + first;
    cudaDeviceSynchronize();
  }

  cudaMemcpy( host_vector, vector, vector_size*sizeof(T), cudaMemcpyDeviceToHost );
  check_cuda_error("Memcpy in reduction");

  for(int i=0; i<vector_size; i++){
    sum += host_vector[i];
  }

  free(host_vector);

  return sum;
}


template<typename T>
__global__ void cuda_reduce_prod_kernel( T * vector, int vector_size, int new_size, int elems)
{
  int Index = threadIdx.x + blockIdx.x * blockDim.x;
  if( Index < new_size ){
    for(int i=1; i<elems; i++){
      int ind = Index + i*new_size;
      vector[Index] *= vector[ ind ];
    }
  }
}

template<typename T>
T cuda_reduce_prod(  T * vector, int N ){
  const int reduce_step = 32;
  T prod=1;
  T * host_vector = (T *)malloc(N*sizeof(T));
  int vector_size = N;

  while( vector_size > reduce_step ){
    // Take the last n elements that are divisible by reduce_step
    int first = vector_size%reduce_step;
    // Calculate the size of the reduced list
    int new_size = (vector_size-first)/reduce_step;
    // Find number of blocks and launch the kernel
    int blocks = new_size/N_threads + 1;
    cuda_reduce_prod_kernel<<< blocks, N_threads >>>( vector+first, vector_size, new_size, reduce_step );
    // Find the full size of the resulting array
    vector_size = new_size + first;
    cudaDeviceSynchronize();
  }

  check_cuda_error("cuda_reduce_prod kernel");
  cudaMemcpy( host_vector, vector, vector_size*sizeof(T), cudaMemcpyDeviceToHost );
  check_cuda_error("Memcpy in reduction");

  for(int i=0; i<vector_size; i++){
    prod *= host_vector[i];
  }

  free(host_vector);

  return prod;
}



template<typename T>
void cuda_multireduce_sum( std::vector<T> &vector, T * d_array, int N ){
  for( int v=0; v<vector.size(); v++ ){
    vector[v] += cuda_reduce_sum( d_array + v*N, N );
  }
}


template<typename T>
void cuda_multireduce_product( std::vector<T> vector, T * d_array, int N ){
  for( int v=0; v<vector.size(); v++ ){
    vector[v] += cuda_reduce_product( d_array + v*N, N );
  }
}


template<typename T>
__global__ void cuda_set_zero_kernel( T * vector, int elems)
{
  int Index = threadIdx.x + blockIdx.x * blockDim.x;
  if( Index < elems ){
    vector[Index] = 0;
  }
}

#if 0
template<typename T>
__global__ void cuda_set_one_kernel( T * vector, int elems)
{
  int Index = threadIdx.x + blockIdx.x * blockDim.x;
  if( Index < elems ){
    vector[Index] = 0;
  }
}

template<typename T>
void cuda_set_one_kernel( T * vec, size_t N ){
  int blocks = N/N_threads + 1;
  cuda_set_one_kernel<<<blocks, N_threads>>>(vec, N);
}

#endif // 0

template<typename T>
void cuda_set_zero( T * vec, size_t N ){
  int blocks = N/N_threads + 1;
  cuda_set_zero_kernel<<<blocks, N_threads>>>(vec, N);
}



inline void synchronize_threads(){
  cudaDeviceSynchronize();
}


void initialize_cuda(int rank);
void cuda_device_info();


#else
// This is not the CUDA compiler
// Maybe hilapp?

// define cuda functions in order to avoid compilation errors
// in hilapp
#define hila_random() 0.5
#define seed_random(seed) while(0)

#define cudaMalloc(a,b) 0
#define cudaFree(a) 0

#define check_cuda_error(a)
#define check_cuda_error_code(c,a)



inline void synchronize_threads(){};
void initialize_cuda(int rank){};
void cuda_device_info();



#endif





/// Implements test for basic in types, similar to
/// std::is_arithmetic, but allows the backend to add
/// it's own basic tyes (such as AVX vectors)
template< class T >
struct is_arithmetic : std::integral_constant<
  bool,
  std::is_arithmetic<T>::value
> {};

template< class T, class U >
struct is_assignable : std::integral_constant<
  bool,
  std::is_assignable<T,U>::value
> {};


#endif
