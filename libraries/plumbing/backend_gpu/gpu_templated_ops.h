#ifndef GPU_TEMPLATED_OPS_H
#define GPU_TEMPLATED_OPS_H

#include "plumbing/defs.h"
#include "plumbing/backend_gpu/defs.h"

// this include has to be after the backend defs, because those define hila::random()
#include "plumbing/random.h"

#include "plumbing/type_tools.h"

// #if defined(__CUDACC__) || defined(__HIPCC__)

#if !defined(HILAPP)
#if defined(CUDA)
#include <cub/cub.cuh>
namespace gpucub = cub;
using gpuStream_t = cudaStream_t;

#endif

// Ensure no #pragma directives are embedded within macro arguments in the code.

#if defined(HIP)
#include <hipcub/hipcub.hpp>
namespace gpucub = hipcub;
using gpuStream_t = hipStream_t;
#endif

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
    cudaMalloc( &d_sum, sizeof(T) );
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

// Functions used by hilapp in reductions
// (cannot guarantee operator= and operator+= are marked __device__)

template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
__device__ static void _hila_kernel_set_zero(T &v) {
    v = 0;
}

template <typename T, std::enable_if_t<!hila::is_arithmetic<T>::value, int> = 0>
__device__ static void _hila_kernel_set_zero(T &v) {
    using ntype = hila::arithmetic_type<T>;
    constexpr int N = sizeof(T) / sizeof(ntype);

    ntype *arr = reinterpret_cast<ntype *>(&v);
    for (int i = 0; i < N; i++)
        arr[i] = 0;
}

template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
__device__ static void _hila_kernel_copy_var(T &a, const T &b) {
    a = b;
}

template <typename T, std::enable_if_t<!hila::is_arithmetic<T>::value, int> = 0>
__device__ static void _hila_kernel_copy_var(T &a, const T &b) {
    using ntype = hila::arithmetic_type<T>;
    constexpr int N = sizeof(T) / sizeof(ntype);

    ntype *ar = reinterpret_cast<ntype *>(&a);
    const ntype *br = reinterpret_cast<const ntype *>(&b);

    for (int i = 0; i < N; i++)
        ar[i] = br[i];
}

template <typename T, std::enable_if_t<hila::is_arithmetic<T>::value, int> = 0>
__device__ static void _hila_kernel_add_var(T &a, const T &b) {
    a += b;
}

template <typename T, std::enable_if_t<!hila::is_arithmetic<T>::value, int> = 0>
__device__ static void _hila_kernel_add_var(T &a, const T &b) {
    using ntype = hila::arithmetic_type<T>;
    constexpr int N = sizeof(T) / sizeof(ntype);

    ntype *ar = reinterpret_cast<ntype *>(&a);
    const ntype *br = reinterpret_cast<const ntype *>(&b);

    for (int i = 0; i < N; i++)
        ar[i] += br[i];
}


////////////////////////////////////////////////////////////

// A simple hand-written reduction that does not require a library
template <typename T>
__global__ void gpu_reduce_sum_kernel(T *vector, int vector_size, int new_size, int elems) {
    unsigned Index = threadIdx.x + blockIdx.x * blockDim.x;
    if (Index < new_size) {
        for (int i = 1; i < elems; i++) {
            int ind = Index + i * new_size;
            // vector[Index] += vector[ind];
            _hila_kernel_add_var(vector[Index], vector[ind]);
        }
    }
}

template <typename T>
T gpu_reduce_sum(T *vector, int N) {
    const int reduce_step = 32;
    T sum = 0;
    T *host_vector = (T *)memalloc(N * sizeof(T));
    int vector_size = N;

    while (vector_size > reduce_step) {
        // Take the last n elements that are divisible by reduce_step
        int first = vector_size % reduce_step;
        // Calculate the size of the reduced list
        int new_size = (vector_size - first) / reduce_step;
        // Find number of blocks and launch the kernel
        int blocks = (new_size - 1) / N_threads + 1;
        gpu_reduce_sum_kernel<<<blocks, N_threads>>>(vector + first, vector_size, new_size,
                                                     reduce_step);
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
__global__ void gpu_reduce_product_kernel(T *vector, int vector_size, int new_size, int elems) {
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
        gpu_reduce_product_kernel<<<blocks, N_threads>>>(vector + first, vector_size, new_size,
                                                         reduce_step);

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

/// Float atomicAdd should exist

__device__ inline float atomic_Add(float *dp, float v) {
    return atomicAdd(dp, v);
}


/// Do addition for "long long" -sized int type, which in cuda is 64 bits
/// Cuda includes a ULL atomicAdd, which is used here.  Signed
/// version can be obtained just by casting to ULL values, and
/// using ULL atomicAdd - this gives the correct signed addition
/// with "wraparound" overflow behaviour

template <typename T, typename B,
          std::enable_if_t<std::is_integral<T>::value && sizeof(T) == sizeof(unsigned long long) &&
                               std::is_convertible<B, T>::value,
                           int> = 0>
__device__ inline T atomicAdd(T *dp, B v) {

    T tv = v;
    return atomicAdd((unsigned long long int *)dp, (unsigned long long int)tv);
}

/// Atomic add for composed datatypes - do element-by-element
/// requires that hila::arithmetic_type is defined
template <
    typename T, typename B,
    std::enable_if_t<!std::is_arithmetic<T>::value && std::is_assignable<T &, T>::value, int> = 0>
__device__ inline void atomicAdd(T *d, const B &bv) {

    T v;
    v = bv;
    hila::arithmetic_type<T> *dp;
    const hila::arithmetic_type<T> *dv;
    constexpr int N = sizeof(T) / sizeof(hila::arithmetic_type<T>);

    dp = (hila::arithmetic_type<T> *)(void *)d;
    dv = (hila::arithmetic_type<T> *)(void *)&v;

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

///////////////////////////////////////////////////////////////////////////////

template <typename T>
__global__ void sum_blocked_vectorreduction_kernel(T *D, const int reduction_size,
                                                   const int threads) {
    int id = threadIdx.x + blockIdx.x * blockDim.x;

    T sum;

    if (id < reduction_size) {
        // id is now the reduction coordinate
        sum = D[id];
        for (int i = 1; i < threads; i++) {
            // add everything to zero
            sum += D[id + i * reduction_size];
        }
        D[id] = sum;
    }

}


template <typename T>
void sum_blocked_vectorreduction(T *data, const int reduction_size, const int threads) {

    // straightforward implementation, use as many threads as elements in reduction vector

    int blocks = (reduction_size + N_threads - 1) / N_threads;
    T* host_data = (T*)malloc(reduction_size*threads * sizeof(T));
    gpuMemcpy(host_data, data, reduction_size*threads * sizeof(T), gpuMemcpyDeviceToHost);

    std::ofstream output_file("reduction_output.txt");
    for (int i = 0; i < reduction_size * threads; ++i) {
        output_file << host_data[i] << " ";
        if ((i + 1) % reduction_size == 0) {
            output_file << " ; ";
        }
        if ((i + 1) % threads == 0) {
            output_file << std::endl ;
        }
    }
    output_file << std::endl;
    output_file.close();

    free(host_data);
    sum_blocked_vectorreduction_kernel<<<blocks, N_threads>>>(data, reduction_size, threads);

    check_device_error("sum_blocked_vectorreduction");
}

template <typename T>
void sum_blocked_vectorreduction_cub(T *data, const int reduction_size, const int threads) {

    // straightforward implementation, use as many threads as elements in reduction vector

    int blocks = (reduction_size + N_threads - 1) / N_threads;
    T* host_data = (T*)malloc(reduction_size*threads * sizeof(T));
    gpuMemcpy(host_data, data, reduction_size*threads * sizeof(T), gpuMemcpyDeviceToHost);

    std::ofstream output_file("reduction_output.txt");
    for (int i = 0; i < reduction_size * threads; ++i) {
        output_file << host_data[i] << " ";
        if ((i + 1) % reduction_size == 0) {
            output_file << " ; ";
        }
        if ((i + 1) % threads == 0) {
            output_file << std::endl ;
        }
    }
    output_file << std::endl;
    output_file.close();

    free(host_data);

    size_t temp_storage_bytes = 0;
    //sum_blocked_vectorreduction_kernel<<<blocks, N_threads>>>(data, reduction_size, threads);
    T* temp_storage = nullptr;
    int stream_count = 16;
    std::vector<gpuStream_t> streams(stream_count);
    for (int i = 0; i < stream_count; i++) {
        gpuStreamCreate(&streams[i]);
    }
    for (int z = 0; z < reduction_size; z++) {
        T* z_index_ptr = data + z * threads;
        int stream_id = z % stream_count;
        gpucub::DeviceReduce::Sum(temp_storage, temp_storage_bytes, z_index_ptr, z_index_ptr, threads,streams[stream_id]);
    }
    
    for (int i = 0; i < stream_count; i++) {
        gpuStreamSynchronize(streams[i]);
        gpuStreamDestroy(streams[i]);
    }
    for (int z = 0; z < reduction_size; z++) {
        int z_threads_stride = z * threads;
        data[z] = data[z_threads_stride];
    }
    check_device_error("sum_blocked_vectorreduction");
}


///////////////////////

template <typename T>
__global__ void gpu_set_value_kernel(T *vector, T value, int elems) {
    int Index = threadIdx.x + blockIdx.x * blockDim.x;
    if (Index < elems) {
        vector[Index] = value;
    }
}

// passing a ptr instead of value directly
template <typename T>
__global__ void gpu_set_value_kernel_ptr(T *vector, const T* valptr, int elems) {
    int Index = threadIdx.x + blockIdx.x * blockDim.x;
    if (Index < elems) {
        vector[Index] = *valptr;
    }
}


// template <typename T>
// __global__ void gpu_set_zero_kernel(T *vector, int elems) {
//     unsigned Index = threadIdx.x + blockIdx.x * blockDim.x;
//     if (Index < elems) {
//         vector[Index] = 0;
//     }
// }


/// use memset to set value to zero - not useful for other values
template <typename T>
inline void gpu_set_zero(T *vec, size_t N) {
    gpuMemset(vec, 0, N * sizeof(T));
}

/// Setting to some value needs to be done elem-by-elem
template <typename T>
void gpu_set_value(T *vec, const T &val, size_t N) {
    int blocks = N / N_threads + 1;
    if constexpr (sizeof(T) <= GPU_GLOBAL_ARG_MAX_SIZE)
        // small T size, pass as arg to __global__
        gpu_set_value_kernel<<<blocks, N_threads>>>(vec, val, N);
    else {
        // bigger size, memcopy
        T *buf;
        gpuMalloc(&buf, sizeof(T));
        gpuMemcpy(buf, (char *)&val, sizeof(T), gpuMemcpyHostToDevice);

        // call the kernel to set correct indexes
        gpu_set_value_kernel_ptr<<<blocks, N_threads>>>(vec, buf, N);
        gpuFree(buf);
    }
}

// template <typename T>
// void gpu_set_zero(T *vec, size_t N) {
//     int blocks = N / N_threads + 1;
//     gpu_set_zero_kernel<<<blocks, N_threads>>>(vec, N);
// }


#endif // !HILAPP

#endif
