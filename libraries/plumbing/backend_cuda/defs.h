#ifndef CUDA_INCLUDE_H
#define CUDA_INCLUDE_H

// On Puhti, use UCX_MEMTYPE_CACHE=n with
// CUDA_AWARE_MPI

#define CUDA_AWARE_MPI

#include <sstream>
#include <iostream>

// Prototypes for memory pool ops 
void gpu_memory_pool_alloc(void **p, size_t req_size);
void gpu_memory_pool_free(void * ptr);
void gpu_memory_pool_purge();
void gpu_memory_pool_report();


#ifndef HILAPP

// GPU specific definitions

////////////////////////////////////////////////////////////////////////////////////
// Some cuda-specific definitions
////////////////////////////////////////////////////////////////////////////////////
#if defined(CUDA)

#include <cuda.h>
#include <cuda_runtime.h>
//#include <cub/cub.cuh>

// THis set in params.h now
// #define N_threads 256 // Threads per block for CUDA   TODO: make configurable?

using gpuError = cudaError;
#define gpuSuccess cudaSuccess

/////////////////////////////////////////////
// If gpu memory pool in use, the interface to memory
#ifdef GPU_MEMORY_POOL
#define gpuMalloc(a,b) gpu_memory_pool_alloc(a,b)
#define gpuFree(a) gpu_memory_pool_free(a)
#define gpuMemPoolPurge() gpu_memory_pool_purge()
#define gpuMemPoolReport() gpu_memory_pool_report()

#else 
// here std interfaces

#define gpuMemPoolPurge() do {} while(0)
#define gpuMemPoolReport() do {} while(0)

#ifdef CUDA_MALLOC_ASYNC
#define gpuMalloc(a, b) GPU_CHECK(cudaMallocAsync(a, b, 0))
#define gpuFree(a) GPU_CHECK(cudaFreeAsync(a, 0))

#else
#define gpuMalloc(a, b) GPU_CHECK(cudaMalloc(a, b))
#define gpuFree(a) GPU_CHECK(cudaFree(a))

#endif

#endif  // gpu memory pool
/////////////////////////////////////////////


#define gpuGetLastError cudaGetLastError
#define gpuMemcpy(a, b, c, d) GPU_CHECK(cudaMemcpy(a, b, c, d))
#define gpuMemcpyHostToDevice cudaMemcpyHostToDevice
#define gpuMemcpyDeviceToHost cudaMemcpyDeviceToHost
#define gpuDeviceSynchronize() GPU_CHECK(cudaDeviceSynchronize())

#define GPUTYPESTR "CUDA"

#ifdef __CUDA_ARCH__
#define __GPU_DEVICE_COMPILE__ __CUDA_ARCH__
#endif

////////////////////////////////////////////////////////////////////////////////////
// Same for HIP
////////////////////////////////////////////////////////////////////////////////////
#elif defined(HIP)

#include <hip/hip_runtime.h>
#include <hiprand.h>

//#include <hipcub/hipcub.hpp>*

// Set in params.h now
// #define N_threads 256 // Threads per block for CUDAs

using gpuError = hipError_t;
#define gpuSuccess hipSuccess

/////////////////////////////////////////////
// If gpu memory pool in use, the interface to memory
#ifdef GPU_MEMORY_POOL
#define gpuMalloc(a,b) gpu_memory_pool_alloc(a,b)
#define gpuFree(a) gpu_memory_pool_free(a)
#define gpuMemPoolPurge() gpu_memory_pool_purge()

#else
// here std interfaces

#define gpuMemPoolPurge() do {} while(0)

#define gpuMalloc(a, b) GPU_CHECK(hipMalloc(a, b))
#define gpuFree(a) GPU_CHECK(hipFree(a))

#endif // ifdef memory pool

#define gpuGetLastError hipGetLastError
#define gpuMemcpy(a, b, siz, d) GPU_CHECK(hipMemcpy(a, b, siz, d))
#define gpuMemcpyHostToDevice hipMemcpyHostToDevice
#define gpuMemcpyDeviceToHost hipMemcpyDeviceToHost
#define gpuDeviceSynchronize() GPU_CHECK(hipDeviceSynchronize())

#define GPUTYPESTR "HIP"

#ifdef __HIP_DEVICE_COMPILE__
#define __GPU_DEVICE_COMPILE__ __HIP_DEVICE_COMPILE__
#endif

#endif
////////////////////////////////////////////////////////////////////////////////////
// General GPU (cuda/hip) definitions
////////////////////////////////////////////////////////////////////////////////////

namespace hila {
__device__ __host__ double random();
void seed_device_rng(unsigned long seed);
} // namespace hila

#define GPU_CHECK(cmd)                                                           \
    do {                                                                               \
        auto code = cmd;                                                               \
        gpu_exit_on_error(code, #cmd, __FILE__, __LINE__);                             \
    } while (0)

#define check_device_error(msg) gpu_exit_on_error(msg, __FILE__, __LINE__)
#define check_device_error_code(code, msg)                                             \
    gpu_exit_on_error(code, msg, __FILE__, __LINE__)
void gpu_exit_on_error(const char *msg, const char *file, int line);
void gpu_exit_on_error(gpuError code, const char *msg, const char *file, int line);

inline void synchronize_threads() {
    gpuDeviceSynchronize();
 } 
 

#else  // NOW HILAPP

////////////////////////////////////////////////////////////////////////////////////
// Now not cuda or hip - hilapp stage scans this section
///////////////////////////////////////////////////////////////////////////////////

namespace hila {
double random();
void seed_device_rng(unsigned long seed);
} // namespace hila

using gpuError = int;

// Define empty stubs - return 1 (true)
#define gpuMalloc(a, b) do {} while(0)
#define gpuFree(a) do {} while(0)
#define gpuMemcpy(a, b, siz, d) do {} while(0)
#define gpuMemPoolPurge() do {} while(0)
#define gpuMemPoolReport() do {} while(0)

#define check_device_error(msg) do {} while(0)
#define check_device_error_code(code, msg) do {} while(0)

#define GPUTYPESTR "NONE"

inline void synchronize_threads() {}

#endif
////////////////////////////////////////////////////////////////////////////////////

void initialize_cuda(int rank);
void gpu_device_info();

// This is not the CUDA compiler
// Maybe hilapp?

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
