#ifndef GPU_REDUCTION_H_
#define GPU_REDUCTION_H_

// Slightly modified from the NVIDIA reduction code
//   - Make kernel number to be a fixed constant
//   - Make number of threads to be a template parameter
//   - change the main interface

/* Copyright (c) 2021, NVIDIA CORPORATION. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of NVIDIA CORPORATION nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
    Parallel reduction kernels
*/

//#include "hila.h"

#if !defined(HILAPP) && !defined(SLOW_GPU_REDUCTION)

#define _CG_ABI_EXPERIMENTAL
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
static constexpr int whichKernel = GPU_REDUCE_KERNEL;
static constexpr int numThreads = N_GPU_REDUCE_THREADS;

// Define what reduction kernel to use - a local variable
// Number from 0 to 9, can benchmark ...
// TODO: make this makefile-define!

// Utility class used to avoid linker errors with extern
// unsized shared memory arrays with templated type
template <class T>
struct SharedMemory {
    __device__ inline operator T *() {
        extern __shared__ int __smem[];
        return (T *)__smem;
    }

    __device__ inline operator const T *() const {
        extern __shared__ int __smem[];
        return (T *)__smem;
    }
};

// specialize for double to avoid unaligned memory
// access compile errors
template <>
struct SharedMemory<double> {
    __device__ inline operator double *() {
        extern __shared__ double __smem_d[];
        return (double *)__smem_d;
    }

    __device__ inline operator const double *() const {
        extern __shared__ double __smem_d[];
        return (double *)__smem_d;
    }
};

template <class T>
__device__ __forceinline__ T warpReduceSum(unsigned int mask, T mySum) {
    for (int offset = warpSize / 2; offset > 0; offset /= 2) {
        mySum += __shfl_down_sync(mask, mySum, offset);
    }
    return mySum;
}

#if __CUDA_ARCH__ >= 800
// Specialize warpReduceFunc for int inputs to use __reduce_add_sync intrinsic
// when on SM 8.0 or higher
template <>
__device__ __forceinline__ int warpReduceSum<int>(unsigned int mask, int mySum) {
    mySum = __reduce_add_sync(mask, mySum);
    return mySum;
}
#endif

/*
    Parallel sum reduction using shared memory
    - takes log(n) steps for n input elements
    - uses n threads
    - only works for power-of-2 arrays
*/

/* This reduction interleaves which threads are active by using the modulo
   operator.  This operator is very expensive on GPUs, and the interleaved
   inactivity means that no whole warps are active, which is also very
   inefficient */
template <class T>
__global__ void reduce0(T *g_idata, T *g_odata, unsigned int n) {
    // Handle to thread block group
    cooperative_groups::thread_block cta = cooperative_groups::this_thread_block();
    T *sdata = SharedMemory<T>();

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

    sdata[tid] = (i < n) ? g_idata[i] : 0;

    cooperative_groups::sync(cta);

    // do reduction in shared mem
    for (unsigned int s = 1; s < blockDim.x; s *= 2) {
        // modulo arithmetic is slow!
        if ((tid % (2 * s)) == 0) {
            sdata[tid] += sdata[tid + s];
        }

        cooperative_groups::sync(cta);
    }

    // write result for this block to global mem
    if (tid == 0)
        g_odata[blockIdx.x] = sdata[0];
}

/* This version uses contiguous threads, but its interleaved
   addressing results in many shared memory bank conflicts.
*/
template <class T>
__global__ void reduce1(T *g_idata, T *g_odata, unsigned int n) {
    // Handle to thread block group
    cooperative_groups::thread_block cta = cooperative_groups::this_thread_block();
    T *sdata = SharedMemory<T>();

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

    sdata[tid] = (i < n) ? g_idata[i] : 0;

    cooperative_groups::sync(cta);

    // do reduction in shared mem
    for (unsigned int s = 1; s < blockDim.x; s *= 2) {
        int index = 2 * s * tid;

        if (index < blockDim.x) {
            sdata[index] += sdata[index + s];
        }

        cooperative_groups::sync(cta);
    }

    // write result for this block to global mem
    if (tid == 0)
        g_odata[blockIdx.x] = sdata[0];
}

/*
    This version uses sequential addressing -- no divergence or bank conflicts.
*/
template <class T>
__global__ void reduce2(T *g_idata, T *g_odata, unsigned int n) {
    // Handle to thread block group
    cooperative_groups::thread_block cta = cooperative_groups::this_thread_block();
    T *sdata = SharedMemory<T>();

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

    sdata[tid] = (i < n) ? g_idata[i] : 0;

    cooperative_groups::sync(cta);

    // do reduction in shared mem
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] += sdata[tid + s];
        }

        cooperative_groups::sync(cta);
    }

    // write result for this block to global mem
    if (tid == 0)
        g_odata[blockIdx.x] = sdata[0];
}

/*
    This version uses n/2 threads --
    it performs the first level of reduction when reading from global memory.
*/
template <class T>
__global__ void reduce3(T *g_idata, T *g_odata, unsigned int n) {
    // Handle to thread block group
    cooperative_groups::thread_block cta = cooperative_groups::this_thread_block();
    T *sdata = SharedMemory<T>();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * (blockDim.x * 2) + threadIdx.x;

    T mySum = (i < n) ? g_idata[i] : 0;

    if (i + blockDim.x < n)
        mySum += g_idata[i + blockDim.x];

    sdata[tid] = mySum;
    cooperative_groups::sync(cta);

    // do reduction in shared mem
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] = mySum = mySum + sdata[tid + s];
        }

        cooperative_groups::sync(cta);
    }

    // write result for this block to global mem
    if (tid == 0)
        g_odata[blockIdx.x] = mySum;
}

/*
    This version uses the warp shuffle operation if available to reduce
    warp synchronization. When shuffle is not available the final warp's
    worth of work is unrolled to reduce looping overhead.

    See
   http://devblogs.nvidia.com/parallelforall/faster-parallel-reductions-kepler/
    for additional information about using shuffle to perform a reduction
    within a warp.

    Note, this kernel needs a minimum of 64*sizeof(T) bytes of shared memory.
    In other words if blockSize <= 32, allocate 64*sizeof(T) bytes.
    If blockSize > 32, allocate blockSize*sizeof(T) bytes.
*/
template <class T, unsigned int blockSize>
__global__ void reduce4(T *g_idata, T *g_odata, unsigned int n) {
    // Handle to thread block group
    cooperative_groups::thread_block cta = cooperative_groups::this_thread_block();
    T *sdata = SharedMemory<T>();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * (blockDim.x * 2) + threadIdx.x;

    T mySum = (i < n) ? g_idata[i] : 0;

    if (i + blockSize < n)
        mySum += g_idata[i + blockSize];

    sdata[tid] = mySum;
    cooperative_groups::sync(cta);

    // do reduction in shared mem
    for (unsigned int s = blockDim.x / 2; s > 32; s >>= 1) {
        if (tid < s) {
            sdata[tid] = mySum = mySum + sdata[tid + s];
        }

        cooperative_groups::sync(cta);
    }

    cooperative_groups::thread_block_tile<32> tile32 =
        cooperative_groups::tiled_partition<32>(cta);

    if (cta.thread_rank() < 32) {
        // Fetch final intermediate sum from 2nd warp
        if (blockSize >= 64)
            mySum += sdata[tid + 32];
        // Reduce final warp using shuffle
        for (int offset = tile32.size() / 2; offset > 0; offset /= 2) {
            mySum += tile32.shfl_down(mySum, offset);
        }
    }

    // write result for this block to global mem
    if (cta.thread_rank() == 0)
        g_odata[blockIdx.x] = mySum;
}

/*
    This version is completely unrolled, unless warp shuffle is available, then
    shuffle is used within a loop.  It uses a template parameter to achieve
    optimal code for any (power of 2) number of threads.  This requires a switch
    statement in the host code to handle all the different thread block sizes at
    compile time. When shuffle is available, it is used to reduce warp
   synchronization.

    Note, this kernel needs a minimum of 64*sizeof(T) bytes of shared memory.
    In other words if blockSize <= 32, allocate 64*sizeof(T) bytes.
    If blockSize > 32, allocate blockSize*sizeof(T) bytes.
*/
template <class T, unsigned int blockSize>
__global__ void reduce5(T *g_idata, T *g_odata, unsigned int n) {
    // Handle to thread block group
    cooperative_groups::thread_block cta = cooperative_groups::this_thread_block();
    T *sdata = SharedMemory<T>();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * (blockSize * 2) + threadIdx.x;

    T mySum = (i < n) ? g_idata[i] : 0;

    if (i + blockSize < n)
        mySum += g_idata[i + blockSize];

    sdata[tid] = mySum;
    cooperative_groups::sync(cta);

    // do reduction in shared mem
    if ((blockSize >= 512) && (tid < 256)) {
        sdata[tid] = mySum = mySum + sdata[tid + 256];
    }

    cooperative_groups::sync(cta);

    if ((blockSize >= 256) && (tid < 128)) {
        sdata[tid] = mySum = mySum + sdata[tid + 128];
    }

    cooperative_groups::sync(cta);

    if ((blockSize >= 128) && (tid < 64)) {
        sdata[tid] = mySum = mySum + sdata[tid + 64];
    }

    cooperative_groups::sync(cta);

    cooperative_groups::thread_block_tile<32> tile32 =
        cooperative_groups::tiled_partition<32>(cta);

    if (cta.thread_rank() < 32) {
        // Fetch final intermediate sum from 2nd warp
        if (blockSize >= 64)
            mySum += sdata[tid + 32];
        // Reduce final warp using shuffle
        for (int offset = tile32.size() / 2; offset > 0; offset /= 2) {
            mySum += tile32.shfl_down(mySum, offset);
        }
    }

    // write result for this block to global mem
    if (cta.thread_rank() == 0)
        g_odata[blockIdx.x] = mySum;
}

/*
    This version adds multiple elements per thread sequentially.  This reduces
   the overall cost of the algorithm while keeping the work complexity O(n) and
   the step complexity O(log n). (Brent's Theorem optimization)

    Note, this kernel needs a minimum of 64*sizeof(T) bytes of shared memory.
    In other words if blockSize <= 32, allocate 64*sizeof(T) bytes.
    If blockSize > 32, allocate blockSize*sizeof(T) bytes.
*/
template <class T, unsigned int blockSize, bool nIsPow2>
__global__ void reduce6(T *g_idata, T *g_odata, unsigned int n) {
    // Handle to thread block group
    cooperative_groups::thread_block cta = cooperative_groups::this_thread_block();
    T *sdata = SharedMemory<T>();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int gridSize = blockSize * gridDim.x;

    T mySum = 0;

    // we reduce multiple elements per thread.  The number is determined by the
    // number of active thread blocks (via gridDim).  More blocks will result
    // in a larger gridSize and therefore fewer elements per thread
    if (nIsPow2) {
        unsigned int i = blockIdx.x * blockSize * 2 + threadIdx.x;
        gridSize = gridSize << 1;

        while (i < n) {
            mySum += g_idata[i];
            // ensure we don't read out of bounds -- this is optimized away for
            // powerOf2 sized arrays
            if ((i + blockSize) < n) {
                mySum += g_idata[i + blockSize];
            }
            i += gridSize;
        }
    } else {
        unsigned int i = blockIdx.x * blockSize + threadIdx.x;
        while (i < n) {
            mySum += g_idata[i];
            i += gridSize;
        }
    }

    // each thread puts its local sum into shared memory
    sdata[tid] = mySum;
    cooperative_groups::sync(cta);

    // do reduction in shared mem
    if ((blockSize >= 512) && (tid < 256)) {
        sdata[tid] = mySum = mySum + sdata[tid + 256];
    }

    cooperative_groups::sync(cta);

    if ((blockSize >= 256) && (tid < 128)) {
        sdata[tid] = mySum = mySum + sdata[tid + 128];
    }

    cooperative_groups::sync(cta);

    if ((blockSize >= 128) && (tid < 64)) {
        sdata[tid] = mySum = mySum + sdata[tid + 64];
    }

    cooperative_groups::sync(cta);

    cooperative_groups::thread_block_tile<32> tile32 =
        cooperative_groups::tiled_partition<32>(cta);

    if (cta.thread_rank() < 32) {
        // Fetch final intermediate sum from 2nd warp
        if (blockSize >= 64)
            mySum += sdata[tid + 32];
        // Reduce final warp using shuffle
        for (int offset = tile32.size() / 2; offset > 0; offset /= 2) {
            mySum += tile32.shfl_down(mySum, offset);
        }
    }

    // write result for this block to global mem
    if (cta.thread_rank() == 0)
        g_odata[blockIdx.x] = mySum;
}

/*
     Kernel 7
 */

template <typename T, unsigned int blockSize, bool nIsPow2>
__global__ void reduce7(const T *__restrict__ g_idata, T *__restrict__ g_odata,
                        unsigned int n) {
    T *sdata = SharedMemory<T>();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int gridSize = blockSize * gridDim.x;
    unsigned int maskLength = (blockSize & 31); // 31 = warpSize-1
    maskLength = (maskLength > 0) ? (32 - maskLength) : maskLength;
    const unsigned int mask = (0xffffffff) >> maskLength;

    T mySum = 0;

    // we reduce multiple elements per thread.  The number is determined by the
    // number of active thread blocks (via gridDim).  More blocks will result
    // in a larger gridSize and therefore fewer elements per thread
    if (nIsPow2) {
        unsigned int i = blockIdx.x * blockSize * 2 + threadIdx.x;
        gridSize = gridSize << 1;

        while (i < n) {
            mySum += g_idata[i];
            // ensure we don't read out of bounds -- this is optimized away for
            // powerOf2 sized arrays
            if ((i + blockSize) < n) {
                mySum += g_idata[i + blockSize];
            }
            i += gridSize;
        }
    } else {
        unsigned int i = blockIdx.x * blockSize + threadIdx.x;
        while (i < n) {
            mySum += g_idata[i];
            i += gridSize;
        }
    }

    // Reduce within warp using shuffle or reduce_add if T==int & CUDA_ARCH ==
    // SM 8.0
    mySum = warpReduceSum<T>(mask, mySum);

    // each thread puts its local sum into shared memory
    if ((tid % warpSize) == 0) {
        sdata[tid / warpSize] = mySum;
    }

    __syncthreads();

    const unsigned int shmem_extent =
        (blockSize / warpSize) > 0 ? (blockSize / warpSize) : 1;
    const unsigned int ballot_result = __ballot_sync(mask, tid < shmem_extent);
    if (tid < shmem_extent) {
        mySum = sdata[tid];
        // Reduce final warp using shuffle or reduce_add if T==int & CUDA_ARCH ==
        // SM 8.0
        mySum = warpReduceSum<T>(ballot_result, mySum);
    }

    // write result for this block to global mem
    if (tid == 0) {
        g_odata[blockIdx.x] = mySum;
    }
}

/*
    Kernel 8  gc_reduce
 */

// Performs a reduction step and updates numTotal with how many are remaining
template <typename T, typename Group>
__device__ T cg_reduce_n(T in, Group &threads) {
    return cooperative_groups::reduce(threads, in, cooperative_groups::plus<T>());
}

template <class T>
__global__ void cg_reduce(T *g_idata, T *g_odata, unsigned int n) {
    // Shared memory for intermediate steps
    T *sdata = SharedMemory<T>();
    // Handle to thread block group
    cooperative_groups::thread_block cta = cooperative_groups::this_thread_block();
    // Handle to tile in thread block
    cooperative_groups::thread_block_tile<32> tile =
        cooperative_groups::tiled_partition<32>(cta);

    unsigned int ctaSize = cta.size();
    unsigned int numCtas = gridDim.x;
    unsigned int threadRank = cta.thread_rank();
    unsigned int threadIndex = (blockIdx.x * ctaSize) + threadRank;

    T threadVal = 0;
    {
        unsigned int i = threadIndex;
        unsigned int indexStride = (numCtas * ctaSize);
        while (i < n) {
            threadVal += g_idata[i];
            i += indexStride;
        }
        sdata[threadRank] = threadVal;
    }

    // Wait for all tiles to finish and reduce within CTA
    {
        unsigned int ctaSteps = tile.meta_group_size();
        unsigned int ctaIndex = ctaSize >> 1;
        while (ctaIndex >= 32) {
            cta.sync();
            if (threadRank < ctaIndex) {
                threadVal += sdata[threadRank + ctaIndex];
                sdata[threadRank] = threadVal;
            }
            ctaSteps >>= 1;
            ctaIndex >>= 1;
        }
    }

    // Shuffle redux instead of smem redux
    {
        cta.sync();
        if (tile.meta_group_rank() == 0) {
            threadVal = cg_reduce_n(threadVal, tile);
        }
    }

    if (threadRank == 0)
        g_odata[blockIdx.x] = threadVal;
}

/*
    Kernel 9
 */

template <class T, size_t BlockSize, size_t MultiWarpGroupSize>
__global__ void multi_warp_cg_reduce(T *g_idata, T *g_odata, unsigned int n) {
    // Shared memory for intermediate steps
    T *sdata = SharedMemory<T>();
    __shared__ cooperative_groups::experimental::block_tile_memory<sizeof(T), BlockSize>
        scratch;

    // Handle to thread block group
    auto cta = cooperative_groups::experimental::this_thread_block(scratch);
    // Handle to multiWarpTile in thread block
    auto multiWarpTile =
        cooperative_groups::experimental::tiled_partition<MultiWarpGroupSize>(cta);

    unsigned int gridSize = BlockSize * gridDim.x;
    T threadVal = 0;

    // we reduce multiple elements per thread.  The number is determined by the
    // number of active thread blocks (via gridDim).  More blocks will result
    // in a larger gridSize and therefore fewer elements per thread
    int nIsPow2 = !(n & n - 1);
    if (nIsPow2) {
        unsigned int i = blockIdx.x * BlockSize * 2 + threadIdx.x;
        gridSize = gridSize << 1;

        while (i < n) {
            threadVal += g_idata[i];
            // ensure we don't read out of bounds -- this is optimized away for
            // powerOf2 sized arrays
            if ((i + BlockSize) < n) {
                threadVal += g_idata[i + blockDim.x];
            }
            i += gridSize;
        }
    } else {
        unsigned int i = blockIdx.x * BlockSize + threadIdx.x;
        while (i < n) {
            threadVal += g_idata[i];
            i += gridSize;
        }
    }

    threadVal = cg_reduce_n(threadVal, multiWarpTile);

    if (multiWarpTile.thread_rank() == 0) {
        sdata[multiWarpTile.meta_group_rank()] = threadVal;
    }
    cooperative_groups::sync(cta);

    if (threadIdx.x == 0) {
        threadVal = 0;
        for (int i = 0; i < multiWarpTile.meta_group_size(); i++) {
            threadVal += sdata[i];
        }
        g_odata[blockIdx.x] = threadVal;
    }
}

////////////////////////////////////////////////////////////////////////////////

inline bool isPow2(unsigned int x) {
    return ((x & (x - 1)) == 0);
}

////////////////////////////////////////////////////////////////////////////////
// Wrapper function for kernel launch
// Now make kernel number a template parameter, also number of threads
////////////////////////////////////////////////////////////////////////////////

template <class T, int threads>
void reduce_kernel(int size, int blocks, T *d_idata, T *d_odata) {
    dim3 dimBlock(threads, 1, 1);
    dim3 dimGrid(blocks, 1, 1);

    // when there is only one warp per block, we need to allocate two warps
    // worth of shared memory so that we don't index shared memory out of bounds
    int smemSize = (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);

    // choose which of the optimized versions of reduction to launch
    switch (whichKernel) {
    case 0:
        reduce0<T><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;

    case 1:
        reduce1<T><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;

    case 2:
        reduce2<T><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;

    case 3:
        reduce3<T><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;

    case 4:
        reduce4<T, threads><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;

    case 5:
        reduce5<T, threads><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;

    case 6:
        if (isPow2(size)) {
            reduce6<T, threads, true>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        } else {
            reduce6<T, threads, false>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        }

        break;

    case 7:
        // For reduce7 kernel we require only blockSize/warpSize
        // number of elements in shared memory
        smemSize = ((threads / 32) + 1) * sizeof(T);
        if (isPow2(size)) {
            reduce7<T, threads, true>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        } else {
            reduce7<T, threads, false>
                <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        }
        break;

    case 8:
        cg_reduce<T><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;

    case 9:
        constexpr int numOfMultiWarpGroups = 2;
        smemSize = numOfMultiWarpGroups * sizeof(T);

        static_assert(threads >= 64,
                      "thread block size of < 64 is not supported for this kernel");
        multi_warp_cg_reduce<T, threads, threads / numOfMultiWarpGroups>
            <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
    }
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Compute the number of (threads and) blocks to use for the given reduction
// kernel For the kernels >= 3, we set threads / block to the minimum of
// maxThreads and n/2. For kernels < 3, we set to the minimum of maxThreads and
// n.  For kernel 6, we observe the maximum specified number of blocks, because
// each thread in that kernel can process a variable number of elements.
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// This function performs a reduction of the input data
// d_odata is the
////////////////////////////////////////////////////////////////////////////////
template <class T>
T gpu_reduce(int size, T *d_idata, bool keep_buffers) {

    // Reasonable defaults?  TODO: make adjustable
    int maxBlocks = 64; // only for kernels >= 6

    int numBlocks = 0;

    if (whichKernel < 3) {
        numBlocks = (size + numThreads - 1) / numThreads;
    } else {
        numBlocks = (size + (numThreads * 2 - 1)) / (numThreads * 2);
        if (whichKernel >= 6) {
            if (numBlocks < maxBlocks)
                numBlocks = maxBlocks;
        }
    }

    static T *d_odata = nullptr;
    static T *h_odata = nullptr;

    // allocate buffers if not done before
    if (d_odata == nullptr) {
        gpuMalloc((void **)&d_odata, numBlocks * sizeof(T));
    }
    if (h_odata == nullptr) {
        h_odata = (T *)memalloc(numBlocks * sizeof(T));
    }

    // execute the kernel
    reduce_kernel<T, numThreads>(size, numBlocks, d_idata, d_odata);
    check_device_error("reduction");

    // sum partial sums from each block on CPU
    // copy result from device to host
    gpuMemcpy(h_odata, d_odata, numBlocks * sizeof(T), gpuMemcpyDeviceToHost);

    T gpu_result = 0;
    for (int i = 0; i < numBlocks; i++) {
        gpu_result += h_odata[i];
    }

    // release buffers if not needed again
    if (!keep_buffers) {
        free(h_odata);
        gpuFree(d_odata);
        h_odata = nullptr;
        d_odata = nullptr;
    }

    return gpu_result;
}

/////////////////////////////////////////////////////////////////////////////////////////

#else // if !defined(HILAPP) && !defined(SLOW_GPU_REDUCTION)

// just declare the name
template <class T>
T gpu_reduce(int size, T *d_idata, bool keep_buffers);

#endif // ifndef HILAPP

#ifdef USE_MPI
#include "com_mpi.h"
#endif

/////////////////////////////////////////////////////////////////////////////////////////
// Reduce field var over the lattice

template <typename T>
T Field<T>::gpu_reduce_sum(bool allreduce, Parity par, bool do_mpi) const {

#ifndef EVEN_SITES_FIRST
    assert(par == Parity::all &&
           "EVEN_SITES_FIRST neede for gpu reduction with parity");
#endif

    using base_t = hila::number_type<T>;
    constexpr int n = sizeof(T) / sizeof(base_t);

    // fields are stored 1 after another on gpu fields
    // -- this really relies on the data layout!
    T *fb = this->field_buffer();
    // address data with bptr
    base_t *bptr = (base_t *)(fb);
    const lattice_struct *lat = this->fs->lattice;

    // use this union to set the value element by element
    union {
        T value;
        base_t element[n];
    } result;

    unsigned fsize = lat->field_alloc_size();
    unsigned size = lat->mynode.volume();

    if (par == Parity::even) {
        size = lat->mynode.evensites;
    } else if (par == Parity::odd) {
        size = lat->mynode.oddsites;
        bptr += lat->mynode.evensites;
    }

    for (int i = 0; i < n; i++) {
        // keep buffers until last element
        result.element[i] = gpu_reduce<base_t>(size, bptr, i < (n - 1));
        bptr += fsize;
    }

#ifdef USE_MPI
    // we don't always do MPI - not in generated loops
    if (do_mpi) {
        MPI_Datatype dtype;
        int s;

        dtype = get_MPI_number_type<T>(s);
        assert(dtype != MPI_BYTE && "Unknown number_type in gpu_reduce_sum");

        if (allreduce) {
            MPI_Allreduce(MPI_IN_PLACE, &result.value, size, dtype, MPI_SUM,
                          this->fs->lattice->mpi_comm_lat);
        } else {
            MPI_Reduce(MPI_IN_PLACE, &result.value, size, dtype, MPI_SUM, 0,
                       this->fs->lattice->mpi_comm_lat);
        }
    }
#endif

    return result.value;
}

template <class T>
__global__ T minmax_kernel(T *i_data, T min_or_max_out, bool min_or_max) {
    int thIdx = threadIdx.x;
    int gthIdx = thIdx + blockIdx.x*N_threads;
    const int gridSize = N_threads*gridDim.x;
    min_or_max_out = 1;
}

template <typename T>
T Field<T>::gpu_minmax(bool min_or_max) const {

    T *field_data = this->field_buffer();
    T *return_value_d;
    T return_value_h;
    cudaMalloc(&return_value_d, sizeof(T));

    const lattice_struct *lat = this->fs->lattice;
    unsigned const node_system_size = lat->mynode.volume();
    int const gridSize = (node_system_size + N_threads - 1) / N_threads;
    int const blockSize = N_threads;

    minmax_kernel<<<gridSize, blockSize>>>(field_data, return_value_d, min_or_max);
    cudaMemcpy(&return_value_h, return_value_d, sizeof(T), cudaMemcpyDeviceToHost);\
    cudaFree(return_value_d);
    output0 << "test " << return_value_h << '\n';
    return 0;
}

#endif
