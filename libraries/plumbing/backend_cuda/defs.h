#ifndef CUDA_INCLUDE_H
#define CUDA_INCLUDE_H

// On Puhti, use UCX_MEMTYPE_CACHE=n with
// CUDA_AWARE_MPI
#define CUDA_AWARE_MPI

#ifdef __CUDACC__

#include <sstream>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
//#include <cub/cub.cuh>


#define N_threads 256   // Threads per block for CUDA

/* Random number generator */
extern curandState *curandstate;
__device__ extern curandState *d_curandstate;
__global__ void seed_random_kernel(curandState *state, unsigned long seed,
                                   unsigned int stride);
namespace hila {
    __host__ __device__ double random();
    void seed_device_rng(unsigned long seed);
} 

#define check_cuda_error(msg) cuda_exit_on_error(msg, __FILE__, __LINE__)
#define check_cuda_error_code(code, msg) cuda_exit_on_error(code, msg, __FILE__, __LINE__)
void cuda_exit_on_error(const char *msg, const char *file, int line);
void cuda_exit_on_error(cudaError code, const char *msg, const char *file, int line);



inline void synchronize_threads() { cudaDeviceSynchronize(); }

void initialize_cuda(int rank);
void cuda_device_info();

#else
// This is not the CUDA compiler
// Maybe hilapp?

// define cuda functions in order to avoid compilation errors
// in hilapp
namespace hila {
    double random() { return 0.5; }
    void seed_device_rng(unsigned long seed);
}

#define cudaMalloc(a, b) 0
#define cudaFree(a) 0

#define check_cuda_error(a)
#define check_cuda_error_code(c, a)

inline void synchronize_threads(){};
void initialize_cuda(int rank){};
void cuda_device_info();

#endif

namespace hila {

// Implements test for basic in types, similar to
/// std::is_arithmetic, but allows the backend to add
/// it's own basic tyes (such as AVX vectors)
template <class T>
struct is_arithmetic : std::integral_constant<bool, std::is_arithmetic<T>::value> {};

template <class T, class U>
struct is_assignable : std::integral_constant<bool, std::is_assignable<T, U>::value> {};

} // namespace hila

#endif
