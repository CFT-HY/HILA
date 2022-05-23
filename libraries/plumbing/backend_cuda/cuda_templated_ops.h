#ifndef CUDA_TEMPLATED_OPS_H
#define CUDA_TEMPLATED_OPS_H

#include "plumbing/defs.h"
#include "plumbing/backend_cuda/defs.h"

// this include has to be after the backend defs, because those define hila::random()
#include "plumbing/random.h"

#include "plumbing/type_tools.h"

#if defined(__CUDACC__) || defined(__HIPCC__)

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
template <typename T>
__global__ void gpu_reduce_sum_kernel(T *vector, int vector_size, int new_size,
                                      int elems) {
    unsigned Index = threadIdx.x + blockIdx.x * blockDim.x;
    if (Index < new_size) {
        for (int i = 1; i < elems; i++) {
            int ind = Index + i * new_size;
            vector[Index] += vector[ind];
        }
    }
}

template <typename T>
T gpu_reduce_sum(T *vector, int N) {
    const int reduce_step = 32;
    T sum = 0;
    T *host_vector = (T *)malloc(N * sizeof(T));
    int vector_size = N;

    while (vector_size > reduce_step) {
        // Take the last n elements that are divisible by reduce_step
        int first = vector_size % reduce_step;
        // Calculate the size of the reduced list
        int new_size = (vector_size - first) / reduce_step;
        // Find number of blocks and launch the kernel
        int blocks = (new_size - 1) / N_threads + 1;
        gpu_reduce_sum_kernel<<<blocks, N_threads>>>(vector + first, vector_size,
                                                     new_size, reduce_step);
        check_device_error("gpu_reduce_sum kernel");
        // Find the full size of the resulting array
        vector_size = new_size + first;
        gpuDeviceSynchronize();
    }

    gpuMemcpy(host_vector, vector, vector_size * sizeof(T), gpuMemcpyDeviceToHost);

    for (int i = 0; i < vector_size; i++) {
        sum += host_vector[i];
    }

    free(host_vector);

    return sum;
}

template <typename T>
__global__ void gpu_reduce_product_kernel(T *vector, int vector_size, int new_size,
                                          int elems) {
    unsigned Index = threadIdx.x + blockIdx.x * blockDim.x;
    if (Index < new_size) {
        for (int i = 1; i < elems; i++) {
            int ind = Index + i * new_size;
            vector[Index] *= vector[ind];
        }
    }
}

template <typename T>
T gpu_reduce_product(T *vector, int N) {
    const int reduce_step = 32;
    T prod;
    prod = 1;
    T *host_vector = (T *)malloc(N * sizeof(T));
    int vector_size = N;
    while (vector_size > reduce_step) {
        // Take the last n elements that are divisible by reduce_step
        int first = vector_size % reduce_step;
        // Calculate the size of the reduced list
        int new_size = (vector_size - first) / reduce_step;
        // Find number of blocks and launch the kernel
        int blocks = new_size / N_threads + 1;
        gpu_reduce_product_kernel<<<blocks, N_threads>>>(vector + first, vector_size,
                                                         new_size, reduce_step);
        // Find the full size of the resulting array
        vector_size = new_size + first;
        gpuDeviceSynchronize();
    }

    check_device_error("gpu_reduce_product kernel");

    gpuMemcpy(host_vector, vector, vector_size * sizeof(T), gpuMemcpyDeviceToHost);

    for (int i = 0; i < vector_size; i++) {
        prod *= host_vector[i];
    }

    free(host_vector);
    return prod;
}

#if 0

template <typename T>
void gpu_multireduce_sum(std::vector<T> &vector, T *d_array, int N) {
    for (int v = 0; v < vector.size(); v++) {
        vector[v] += gpu_reduce_sum(d_array + v * N, N);
    }
}

template <typename T>
void gpu_multireduce_product(std::vector<T> vector, T *d_array, int N) {
    for (int v = 0; v < vector.size(); v++) {
        vector[v] += gpu_reduce_product(d_array + v * N, N);
    }
}

#endif

/// Implement here atomicAdd for double precision for less than 6.0 capability
/// motivated by the cuda toolkit documentation.
/// atomicCAS(long long *a, long long b, long long v)
/// does the operation *a = (*a == b ? v : *a), i.e. assigns v to *a if *a == b,
/// atomically.  Returns *a.

#if __CUDA_ARCH__ < 600
__device__ inline double atomic_Add(double *dp, double v) {

    unsigned long long int *dp_ull = (unsigned long long int *)dp;
    unsigned long long int old = *dp_ull;
    unsigned long long int av;

    do {
        av = old;
        old = atomicCAS(dp_ull, av, __double_as_longlong(v + __longlong_as_double(av)));

        // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (av != old);

    return __longlong_as_double(old);
}

#else
__device__ inline double atomic_Add(double *dp, double v) {
    return atomicAdd(dp, v);
}

#endif

/// Do addition for "long long" -sized int type, which in cuda is 64 bits
/// Cuda includes a ULL atomicAdd, which is used here.  Signed
/// version can be obtained just by casting to ULL values, and
/// using ULL atomicAdd - this gives the correct signed addition
/// with "wraparound" overflow behaviour

template <typename T, typename B,
          std::enable_if_t<std::is_integral<T>::value &&
                               sizeof(T) == sizeof(unsigned long long) &&
                               std::is_convertible<B, T>::value,
                           int> = 0>
__device__ inline T atomicAdd(T *dp, B v) {

    T tv = v;
    return atomicAdd((unsigned long long int *)dp, (unsigned long long int)tv);
}

/// Atomic add for composed datatypes - do element-by-element
/// requires that hila::number_type is defined
template <
    typename T, typename B,
    std::enable_if_t<!std::is_arithmetic<T>::value && std::is_convertible<B, T>::value,
                     int> = 0>
__device__ inline void atomicAdd(T *d, const B &bv) {

    T v = bv;
    hila::number_type<T> *dp;
    const hila::number_type<T> *dv;
    constexpr int N = sizeof(T) / sizeof(hila::number_type<T>);

    dp = (hila::number_type<T> *)(void *)d;
    dv = (hila::number_type<T> *)(void *)&v;

    for (int i = 0; i < N; ++i) {
        atomic_Add(dp + i, dv[i]);
    }
}

/// Multiply 2 doubles atomically
__device__ inline double atomicMultiply(double *dp, double v) {

    unsigned long long int *dp_ull = (unsigned long long int *)dp;
    unsigned long long int old = *dp_ull;
    unsigned long long int av;

    do {
        av = old;
        old = atomicCAS(dp_ull, av, __double_as_longlong(v * __longlong_as_double(av)));

        // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (av != old);

    return __longlong_as_double(old);
}

/// Multiply 2 floats atomically
__device__ inline float atomicMultiply(float *dp, float v) {
    unsigned int *dp_ui = (unsigned int *)dp;
    unsigned int old = *dp_ui;
    unsigned int av;

    do {
        av = old;
        old = atomicCAS(dp_ui, av, __float_as_int(v * __int_as_float(av)));

        // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (av != old);

    return __int_as_float(old);
}

///////////////////////

template <typename T>
__global__ void gpu_set_one_kernel(T *vector, int elems) {
    unsigned Index = threadIdx.x + blockIdx.x * blockDim.x;
    if (Index < elems) {
        vector[Index] = 1;
    }
}

template <typename T>
__global__ void gpu_set_zero_kernel(T *vector, int elems) {
    unsigned Index = threadIdx.x + blockIdx.x * blockDim.x;
    if (Index < elems) {
        vector[Index] = 0;
    }
}

template <typename T>
void gpu_set_one(T *vec, size_t N) {
    int blocks = N / N_threads + 1;
    gpu_set_one_kernel<<<blocks, N_threads>>>(vec, N);
}

template <typename T>
void gpu_set_zero(T *vec, size_t N) {
    int blocks = N / N_threads + 1;
    gpu_set_zero_kernel<<<blocks, N_threads>>>(vec, N);
}

#endif // __CUDACC__

#endif
