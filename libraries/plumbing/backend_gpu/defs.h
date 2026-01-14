#ifndef HILA_GPU_DEFS_H
#define HILA_GPU_DEFS_H

// On Puhti, use UCX_MEMTYPE_CACHE=n with
// GPU_AWARE_MPI

#include <sstream>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <memory>
#include <utility>
#include <queue>

// Prototypes for memory pool ops
void gpu_memory_pool_alloc(void **p, size_t req_size);
void gpu_memory_pool_free(void *ptr);
void gpu_memory_pool_purge();
void gpu_memory_pool_report();

////////////////////////////////////////////////////////////////////////////////////
// some device rng headers
////////////////////////////////////////////////////////////////////////////////////
namespace hila {
// double random();  // defined in random.h
void seed_device_rng(unsigned long long seed);
} // namespace hila

namespace hila {
void free_device_rng();
} // namespace hila


#ifndef HILAPP

// GPU specific definitions

////////////////////////////////////////////////////////////////////////////////////
// Some cuda-specific definitions
////////////////////////////////////////////////////////////////////////////////////
#if defined(CUDA)

#include <cuda.h>
#include <cuda_runtime.h>
#include <cub/cub.cuh>

using gpuError = cudaError;
#define gpuSuccess cudaSuccess

/////////////////////////////////////////////
// If gpu memory pool in use, the interface to memory
#ifdef GPU_MEMORY_POOL
#define gpuMalloc(a, b) gpu_memory_pool_alloc((void **)a, b)
#define gpuFree(a) gpu_memory_pool_free(a)
#define gpuMemPoolPurge() gpu_memory_pool_purge()
#define gpuMemPoolReport() gpu_memory_pool_report()

#else
// here std interfaces

// clang-format off
#define gpuMemPoolPurge()  do { } while (0)
#define gpuMemPoolReport() do { } while (0)
// clang-format on

#ifdef CUDA_MALLOC_ASYNC
#define gpuMalloc(a, b) GPU_CHECK(cudaMallocAsync(a, b, 0))
#define gpuFree(a) GPU_CHECK(cudaFreeAsync(a, 0))

#else
#define gpuMalloc(a, b) GPU_CHECK(cudaMalloc((void **)a, b))
#define gpuFree(a) GPU_CHECK(cudaFree(a))

#endif

#endif // gpu memory pool
/////////////////////////////////////////////


#define gpuGetLastError cudaGetLastError
#define gpuMemcpy(a, b, c, d) GPU_CHECK(cudaMemcpy(a, b, c, d))
#define gpuMemcpyHostToDevice cudaMemcpyHostToDevice
#define gpuMemcpyDeviceToHost cudaMemcpyDeviceToHost
#define gpuMemcpyDeviceToDevice cudaMemcpyDeviceToDevice
#define gpuDeviceSynchronize() GPU_CHECK(cudaDeviceSynchronize())
#define gpuStreamSynchronize(a) GPU_CHECK(cudaStreamSynchronize(a))
#define gpuStreamCreate(a) GPU_CHECK(cudaStreamCreate(a))
#define gpuStreamCreateWithFlags(a,b) GPU_CHECK(cudaStreamCreateWithFlags(a,b))
#define gpuStreamQuery(a) cudaStreamQuery(a)
#define gpuStreamDestroy(a) GPU_CHECK(cudaStreamDestroy(a))
#define gpuStreamWaitEvent(a,b,c) GPU_CHECK(cudaStreamWaitEvent(a,b,c))
#define gpuEventCreate(a) GPU_CHECK(cudaEventCreate(a))
#define gpuEventCreateWithFlags(a,b) GPU_CHECK(cudaEventCreateWithFlags(a,b))
#define gpuEventRecord(a,b) GPU_CHECK(cudaEventRecord(a,b))
#define gpuEventQuery(a) cudaEventQuery(a)
#define gpuEventSynchronize(a) GPU_CHECK(cudaEventSynchronize(a))
#define gpuEventDestroy(a) GPU_CHECK(cudaEventDestroy(a))
#define gpuMemset(a, b, c) GPU_CHECK(cudaMemset(a, b, c))
#define gpuMemsetAsync(a, b, c, d) GPU_CHECK(cudaMemsetAsync(a, b, c, d))
#define gpuMemcpyToSymbol(a, b, size, c, dir) GPU_CHECK(cudaMemcpyToSymbol(a, b, size, c, dir))
#define gpuFuncAttributes cudaFuncAttributes
#define gpuFuncGetAttributes cudaFuncGetAttributes
#define gpuSuccess cudaSuccess
#define gpuStream_t cudaStream_t
#define gpuEvent_t cudaEvent_t
#define gpuStreamNonBlocking cudaStreamNonBlocking
#define gpuSetDevice(dev) GPU_CHECK(cudaSetDevice(dev))

#define GPUTYPESTR "CUDA"

#ifdef __CUDA_ARCH__
#define _GPU_DEVICE_COMPILE_ __CUDA_ARCH__
#endif

////////////////////////////////////////////////////////////////////////////////////
// Same for HIP
////////////////////////////////////////////////////////////////////////////////////
#elif defined(HIP)

#include <hip/hip_runtime.h>
#include <hiprand/hiprand.h>

// #include <hipcub/hipcub.hpp>*

using gpuError = hipError_t;
#define gpuSuccess hipSuccess

/////////////////////////////////////////////
// If gpu memory pool in use, the interface to memory
#ifdef GPU_MEMORY_POOL
#define gpuMalloc(a, b) gpu_memory_pool_alloc((void **)a, b)
#define gpuFree(a) gpu_memory_pool_free(a)
#define gpuMemPoolPurge() gpu_memory_pool_purge()
#define gpuMemPoolReport() gpu_memory_pool_report()


#else
// here std interfaces

// clang-format off
#define gpuMemPoolPurge() do {} while (0)
#define gpuMemPoolReport() do {} while (0)
// clang-format on

#define gpuMalloc(a, b) GPU_CHECK(hipMalloc((void **)a, b))
#define gpuFree(a) GPU_CHECK(hipFree(a))

#endif // ifdef memory pool

#define gpuGetLastError hipGetLastError
#define gpuMemcpy(a, b, siz, d) GPU_CHECK(hipMemcpy(a, b, siz, d))
#define gpuMemcpyHostToDevice hipMemcpyHostToDevice
#define gpuMemcpyDeviceToHost hipMemcpyDeviceToHost
#define gpuMemcpyDeviceToDevice hipMemcpyDeviceToDevice
#define gpuDeviceSynchronize() GPU_CHECK(hipDeviceSynchronize())
#define gpuStreamSynchronize(a) GPU_CHECK(hipStreamSynchronize(a))
#define gpuStreamCreate(a) GPU_CHECK(hipStreamCreate(a))
#define gpuStreamCreateWithFlags(a,b) GPU_CHECK(hipStreamCreateWithFlags(a,b))
#define gpuStreamQuery(a) hipStreamQuery(a)
#define gpuStreamDestroy(a) GPU_CHECK(hipStreamDestroy(a))
#define gpuStreamWaitEvent(a,b,c) GPU_CHECK(hipStreamWaitEvent(a,b,c))
#define gpuEventCreate(a) GPU_CHECK(hipEventCreate(a))
#define gpuEventCreateWithFlags(a,b) GPU_CHECK(hipEventCreateWithFlags(a,b))
#define gpuEventRecord(a,b) GPU_CHECK(hipEventRecord(a,b))
#define gpuEventQuery(a) hipEventQuery(a)
#define gpuEventSynchronize(a) GPU_CHECK(hipEventSynchronize(a))
#define gpuEventDestroy(a) GPU_CHECK(hipEventDestroy(a))
#define gpuMemset(a, b, c) GPU_CHECK(hipMemset(a, b, c))
#define gpuMemsetAsync(a, b, c, d) GPU_CHECK(hipMemsetAsync(a, b, c, d))
#define gpuMemcpyToSymbol(a, b, size, c, dir)                                                      \
    GPU_CHECK(hipMemcpyToSymbol(HIP_SYMBOL(a), b, size, c, dir))
#define gpuFuncAttributes hipFuncAttributes
#define gpuFuncGetAttributes hipFuncGetAttributes
#define gpuSuccess hipSuccess
#define gpuStream_t hipStream_t
#define gpuEvent_t hipEvent_t
#define gpuStreamNonBlocking hipStreamNonBlocking
#define gpuSetDevice(dev) GPU_CHECK(hipSetDevice(dev))

#define GPUTYPESTR "HIP"

#ifdef __HIP_DEVICE_COMPILE__
#define _GPU_DEVICE_COMPILE_ __HIP_DEVICE_COMPILE__
#endif
#endif // CUDA or HIP
#ifdef GPU_CCL // if GPU_CCL and CUDA is defined use NCCL

#if defined(CUDA)
#include <nccl.h>
#elif defined(HIP)
#include <rccl/rccl.h> 
#endif

namespace hila {
    void initialize_gccl_communications();
}

template<class Scalar_type>
struct gccl_type;

template<> struct gccl_type<int8_t>   { static constexpr ncclDataType_t value = ncclInt8;  };
template<> struct gccl_type<uint8_t>  { static constexpr ncclDataType_t value = ncclUint8; };
template<> struct gccl_type<int32_t>  { static constexpr ncclDataType_t value = ncclInt32; };
template<> struct gccl_type<uint32_t> { static constexpr ncclDataType_t value = ncclUint32; };
template<> struct gccl_type<int64_t>  { static constexpr ncclDataType_t value = ncclInt64; };
template<> struct gccl_type<uint64_t> { static constexpr ncclDataType_t value = ncclUint64; };
template<> struct gccl_type<float>    { static constexpr ncclDataType_t value = ncclFloat32; };
template<> struct gccl_type<double>   { static constexpr ncclDataType_t value = ncclFloat64; };

#define gcclComm_t ncclComm_t 
#define gcclUniqueId ncclUniqueId
#define gcclGetUniqueId(a) NCCL_CHECK(ncclGetUniqueId(a))
#define gcclCommInitRank(a,b,c,d) NCCL_CHECK(ncclCommInitRank(a,b,c,d))
#define gcclCommDestroy(a) NCCL_CHECK(ncclCommDestroy(a))
#define gcclBroadcast(a,b,c,d,e,f,g) NCCL_CHECK(ncclBroadcast(a,b,c,d,e,f,g))
#define gcclAllReduce(a,b,c,d,e,f,g) NCCL_CHECK(ncclAllReduce(a,b,c,d,e,f,g))
#define gcclReduce(a,b,c,d,e,f,g,h) NCCL_CHECK(ncclReduce(a,b,c,d,e,f,g,h))
#define gcclAllGather(a,b,c,d,e,f) NCCL_CHECK(ncclAllGather(a,b,c,d,e,f))
#define gcclReduceScatter(a,b,c,d,e,f,g) NCCL_CHECK(ncclReduceScatter(a,b,c,d,e,f,g))
#define gcclSend(a,b,c,d,e,f) NCCL_CHECK(ncclSend(a,b,c,d,e,f))
#define gcclRecv(a,b,c,d,e,f) NCCL_CHECK(ncclRecv(a,b,c,d,e,f))
#define gcclGroupStart() NCCL_CHECK(ncclGroupStart())
#define gcclGroupEnd() NCCL_CHECK(ncclGroupEnd())

#define NCCL_CHECK(cmd) \
    { \
        ncclResult_t _result = cmd; \
        if (_result != ncclSuccess) { \
            fprintf(stderr, "NCCL error: %s\n", ncclGetErrorString(_result)); \
            exit(EXIT_FAILURE); \
        } \
    }

#endif // GPU_CCL

////////////////////////////////////////////////////////////////////////////////////
// General GPU (cuda/hip) definitions
////////////////////////////////////////////////////////////////////////////////////


#define GPU_CHECK(cmd)                                                                             \
    do {                                                                                           \
        auto code = cmd;                                                                           \
        gpu_exit_on_error(code, #cmd, __FILE__, __LINE__);                                         \
    } while (0)

#define check_device_error(msg) gpu_exit_on_error(msg, __FILE__, __LINE__)
#define check_device_error_code(code, msg) gpu_exit_on_error(code, msg, __FILE__, __LINE__)
void gpu_exit_on_error(const char *msg, const char *file, int line);
void gpu_exit_on_error(gpuError code, const char *msg, const char *file, int line);

namespace hila {
inline void synchronize_threads() {
    gpuDeviceSynchronize();
}
} // namespace hila

class gpuStreamPool {
public:

    explicit gpuStreamPool() : stream_queue(), active_streams() {
        gpuStream_t stream;
        gpuStreamCreateWithFlags(&stream, gpuStreamNonBlocking);
        stream_queue.push(stream);
    }

    ~gpuStreamPool() {
        while(!stream_queue.empty()) {
            gpuStream_t s = stream_queue.front();
            gpuStreamDestroy(s);
            stream_queue.pop();
        }

        while(!active_streams.empty()) {
            for (auto it = active_streams.begin(); it != active_streams.end(); ) {
                gpuStreamDestroy(*it);
                it = active_streams.erase(it);
            }
        }
        
    }

    gpuStream_t next_stream() {
        if (stream_queue.empty()) {
            gpuStream_t stream;
            gpuStreamCreateWithFlags(&stream, gpuStreamNonBlocking);
            active_streams.push_back(stream);
            return stream;
        } else {
            gpuStream_t stream = stream_queue.front();
            active_streams.push_back(stream);
            stream_queue.pop();
            return stream ;
        }
    }

    void wait_all() {
        while (!active_streams.empty()) {
            for (auto it = active_streams.begin(); it != active_streams.end(); ) {
                if (gpuStreamQuery(*it) == gpuSuccess) {
                    stream_queue.push(*it);
                    it = active_streams.erase(it);
                } else {
                    ++it;
                }
            }
        }
    }

    bool has_active_streams() const {
        return !active_streams.empty();
    }

private:
    static constexpr int DEFAULT_STREAM_POOL_SIZE = 1;
    std::queue<gpuStream_t> stream_queue;
    std::vector<gpuStream_t> active_streams;
};

namespace hila {
    gpuStreamPool& stream_pool();
    gpuStream_t& halo_stream(); 
    gpuEvent_t& halo_event();
    gpuStream_t& bulk_stream();
    gpuEvent_t& bulk_event();
}


#else // NOW HILAPP


///////////////////////////////////////////////////////////////////////////////////
// Now not cuda or hip - hilapp stage scans this section
///////////////////////////////////////////////////////////////////////////////////


using gpuError = int;

// Define empty stubs - return 1 (true)
// clang-format off
#define gpuMalloc(a, b) do {} while(0)
#define gpuFree(a) do {} while(0)
#define gpuMemcpy(a, b, siz, d) do {} while(0)
#define gpuMemcpyHostToDevice 1
#define gpuMemcpyDeviceToHost 2
#define gpuMemset(a,b,c) do {} while(0)
#define gpuMemcpyToSymbol(a, b, size, c, dir) do {} while(0)

#define gpuMemPoolPurge() do {} while(0)
#define gpuMemPoolReport() do {} while(0)

#define check_device_error(msg) do {} while(0)
#define check_device_error_code(code, msg) do {} while(0)

#define gpuStreamSynchronize(a) do {} while(0)
#define gpuDeviceSynchronize() do {} while(0)

#define gpuEventRecord(a,b) do {} while(0)
#define gpuEventSynchronize(a) do {} while(0)


#define gpuGetLastError cudaGetLastError

// Placeholder for hilapp
typedef void* gpuStream_t;
typedef void* gpuEvent_t;
typedef void* gcclComm_t;

template<class Scalar_type>
struct gccl_type;

namespace hila {

struct gpuStreamPool {
    void wait_all() {} 
    gpuStream_t next_stream() {
        return nullptr;
    }
};

inline gpuStreamPool& stream_pool() {
    static gpuStreamPool dummy;
    return dummy;
}

inline gpuStream_t& halo_stream() {
    static gpuStream_t dummy = nullptr;
    return dummy;
}

inline gpuStream_t& bulk_stream() {
    static gpuStream_t dummy = nullptr;
    return dummy;
}

inline gpuStream_t& halo_event() {
    static gpuStream_t dummy = nullptr;
    return dummy;
}

inline gpuStream_t& bulk_event() {
    static gpuStream_t dummy = nullptr;
    return dummy;
}


void initialize_gccl_communications();

} 
// clang-format on


#define GPUTYPESTR "NONE"

namespace hila {
inline void synchronize_threads() {}
} // namespace hila

#endif
////////////////////////////////////////////////////////////////////////////////////

void initialize_gpu(int rank, int device);
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

template <class T>
struct is_floating_point : std::integral_constant<bool, std::is_floating_point<T>::value> {};

} // namespace hila

#endif