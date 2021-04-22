#ifndef CUDA_INCLUDE_H
#define CUDA_INCLUDE_H

// On Puhti, use UCX_MEMTYPE_CACHE=n with
// CUDA_AWARE_MPI
#define CUDA_AWARE_MPI

#include <sstream>
#include <iostream>

#ifndef HILAPP

// GPU specific definitions

////////////////////////////////////////////////////////////////////////////////////
// Some cuda-specific definitions
////////////////////////////////////////////////////////////////////////////////////
#if defined(CUDA)

#include <cuda.h>
#include <cuda_runtime.h>
//#include <cub/cub.cuh>

#define N_threads 256 // Threads per block for CUDA   TODO: make configurable?

using gpuError = cudaError;
#define gpuSuccess cudaSuccess
#define gpuMalloc(a, b) check_gpu_error(cudaMalloc(a, b))
#define gpuFree(a) check_gpu_error(cudaFree(a))
#define gpuGetLastError cudaGetLastError
#define gpuMemcpy(a, b, c, d) check_gpu_error(cudaMemcpy(a, b, c, d))
#define gpuMemcpyHostToDevice cudaMemcpyHostToDevice
#define gpuMemcpyDeviceToHost cudaMemcpyDeviceToHost
#define gpuDeviceSynchronize() check_gpu_error(cudaDeviceSynchronize())

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

//#include <hipcub/hipcub.hpp>

#define N_threads 256 // Threads per block for CUDA

using gpuError = hipError_t;
#define gpuSuccess hipSuccess
#define gpuMalloc(a, b) check_gpu_error(hipMalloc(a, b))
#define gpuFree(a) check_gpu_error(hipFree(a))
#define gpuGetLastError hipGetLastError
#define gpuMemcpy(a, b, siz, d) check_gpu_error(hipMemcpy(a, b, siz, d))
#define gpuMemcpyHostToDevice hipMemcpyHostToDevice
#define gpuMemcpyDeviceToHost hipMemcpyDeviceToHost
#define gpuDeviceSynchronize() check_gpu_error(hipDeviceSynchronize())

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

#define check_gpu_error(cmd)                                                           \
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
 

#else

////////////////////////////////////////////////////////////////////////////////////
// Now not cuda or hip - hilapp stage scans this section
///////////////////////////////////////////////////////////////////////////////////

namespace hila {
double random();
void seed_device_rng(unsigned long seed);
} // namespace hila

using gpuError = int;

// Define empty stubs - return 1 (true)
#define gpuMalloc(a, b) 1
#define gpuFree(a) 1
#define gpuMemcpy(a, b, siz, d) 1
#define check_device_error(msg) 1
#define check_device_error_code(code, msg) 1

#define GPUTYPESTR "NONE"

inline void synchronize_threads() {}

#endif
////////////////////////////////////////////////////////////////////////////////////

void initialize_cuda(int rank);
void cuda_device_info();

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
