///////////////////////////////////////////
/// simple list-based alloc program for cuda/hip

#include "plumbing/memory_pool.h"

// Compile with make .. OPTS="-DPOOL_DEBUG"
// #define POOL_DEBUG

#if !defined(HILAPP)

#if defined(GPU_MEMORY_POOL)

#if defined(HIP)
#define gpuMallocDirect(a, b) GPU_CHECK(hipMalloc(a, b))
#define gpuFreeDirect(a) GPU_CHECK(hipFree(a))
#elif defined(CUDA)
#define gpuMallocDirect(a, b) GPU_CHECK(cudaMalloc(a, b))
#define gpuFreeDirect(a) GPU_CHECK(cudaFree(a))
#else
static_assert(0 && "HIP or CUDA must be defined");
#endif

#if defined(GPU_SHMEM)
#define gpuMallocSharedDirect(a, b)                                                                \
    do {                                                                                           \
        *(a) = nvshmem_malloc(b);                                                                  \
    } while (0)
#define gpuFreeSharedDirect(a) nvshmem_free(a)
#endif
#endif


// keep relatively large min allocation
#define MIN_ALLOC_SIZE 128


void *hila::memory_pool::alloc(size_t req_size) {

    if (req_size < MIN_ALLOC_SIZE) {
        req_size = MIN_ALLOC_SIZE;
    }

    p.n_allocs++;

    // do we have free stuff?  Simple linear search - list should not be too large
    bool found_match = false;
    int steps = 0;
    int first_free_block = -1;
    auto ptr = blocklist.begin();
    for (auto it = blocklist.begin() + p.first_free_block; it != blocklist.end(); it++) {
        steps++;
        if (it->in_use)
            continue;

        if (first_free_block < 0)
            first_free_block = it - blocklist.begin();

        if (it->size == req_size) {
            ptr = it;
            found_match = true;
            break; // perfect match, that's it and stop
        }

        // allow allocated blocks at most twice larger - find smallest
        if (it->size > req_size && it->size < 2 * req_size) {
            if (!found_match || ptr->size > it->size) {
                ptr = it;
            }
            found_match = true;
        }
    }

    p.blocklist_avg_search += steps;

    if (first_free_block < 0) {
        if (found_match)
            first_free_block = ptr - blocklist.begin();
        else
            first_free_block = blocklist.size();
    }

    if (first_free_block > p.first_free_block)
        p.first_free_block = first_free_block;

    if (found_match) {
        // got it
        ptr->in_use = true;

#ifdef POOL_DEBUG
        hila::out << "GPU MEMORY: request " << req_size << " gave block " << ptr->size
                  << " current total " << total_size << '\n';
#endif
        return (ptr->ptr);

    } else {

        // did not find free memory - allocate
        // alloc failure caught by gpuMalloc
        allocation a;
#if defined(CUDA) || defined(HIP)
        if (type == pool_type::STANDARD) {
            gpuMallocDirect(&(a.ptr), req_size);
        } else if (type == pool_type::SHARED) {
#if defined(GPU_SHMEM)
            gpuMallocSharedDirect(&(a.ptr), req_size);
#else
            // Fallback or Error if SHARED requested but SHMEM not compiled
            hila::out << "Error: SHARED pool requested but GPU_SHMEM not defined\n";
            hila::terminate(1);
#endif
        }
#else
        a.ptr = memalloc(req_size);
#endif // GPU_SHMEM
        a.size = req_size;
        a.in_use = true;
        blocklist.push_back(a);

        p.n_true_allocs++;
        p.total_size += req_size;

#ifdef POOL_DEBUG
        hila::out << "GPU MEMORY: request " << req_size << " NEW allocation, current total "
                  << total_size << '\n';
#endif // CUDA OR HIP
        return a.ptr;
    }
}

void hila::memory_pool::free(void *ptr) {

    // search the list backwards for the memory block
    // it's more likely to find the ptr at the end
    for (int i = blocklist.size() - 1; i >= 0; --i) {
        if (blocklist[i].ptr == ptr) {
            if (blocklist[i].in_use == false) {
                hila::out << "Memory free error - ptr " << ptr << " already freed\n";
                hila::terminate(1);
            }
            blocklist[i].in_use = false;

            if (i < p.first_free_block)
                p.first_free_block = i;

#ifdef POOL_DEBUG
            hila::out << "GPU MEMORY: FREE block of size " << it->size << ", current total "
                      << total_size << '\n';
#endif
            return;
        }
    }
    // did not find!  serious error, quit
    hila::out << "Memory free error - unknown pointer  " << ptr << '\n';
    hila::terminate(1);
}


/// Release free memory to the system - avoids extra allocations
void hila::memory_pool::purge() {

    auto ptr = blocklist.begin();
    for (auto it = blocklist.begin(); it != blocklist.end(); it++) {
        if (it->in_use == false) {
#if defined(CUDA) || defined(HIP)
            if (type == pool_type::STANDARD) {
                gpuFreeDirect(it->ptr);
            } else if (type == pool_type::SHARED) {
#if defined(GPU_SHMEM)
                gpuFreeSharedDirect(it->ptr);
#else
                // This path should technically be unreachable if alloc() failed correctly,
                // but we keep it for safety.
                hila::out << "Error: Attempting to purge SHARED block without GPU_SHMEM\n";
                hila::terminate(1);
#endif // GPU_SHMEM
            }
#else
            free(it->ptr);
#endif // CUDA OR HIP
            p.total_size -= it->size;

#ifdef POOL_DEBUG
            hila::out << "GPU MEMORY: Purging " << it->size << ", bytes, total size " << total_size
                      << '\n';
#endif
        } else {
            // copy blocks in use to the beginning of blocklist
            if (ptr != it)
                *ptr = *it;
            ptr++;
        }
    }
    blocklist.resize(ptr - blocklist.begin());
}


#ifdef GPU_MEMORY_POOL

static hila::memory_pool gpu_pool;

void gpu_memory_pool_alloc(void **p, size_t req_size) {
    *p = gpu_pool.alloc(req_size);
}

void gpu_memory_pool_free(void *ptr) {
    gpu_pool.free(ptr);
}

void gpu_memory_pool_purge() {
    gpu_pool.purge();
}

void gpu_memory_pool_report() {
    auto p = gpu_pool.status();
    if_rank0 () {
        hila::out << "\nGPU Memory pool statistics from node 0:\n";
        hila::out << "   Total pool size " << ((double)p.total_size) / (1024 * 1024) << " MB in "
                  << gpu_pool.size() << " blocks\n";
        hila::out << "   # of allocations " << p.n_allocs << "  real allocs "
                  << std::setprecision(2) << ((double)p.n_true_allocs) / p.n_allocs * 100 << "%\n";
        hila::out << "   Average block list search " << (double)p.blocklist_avg_search / p.n_allocs
                  << " steps\n\n";
    }
}

#ifdef GPU_SHMEM
static hila::memory_pool gpu_shared_pool(hila::pool_type::SHARED);

void gpu_shared_memory_pool_alloc(void **p, size_t req_size) {
    *p = gpu_shared_pool.alloc(req_size);
}

void gpu_shared_memory_pool_free(void *ptr) {
    gpu_shared_pool.free(ptr);
}

void gpu_shared_memory_pool_purge() {
    gpu_shared_pool.purge();
}

void gpu_shared_memory_pool_report() {
    auto p = gpu_shared_pool.status();
    if_rank0 () {
        hila::out << "\nGPU SHARED Memory pool statistics from node 0:\n";
        hila::out << "   Total pool size " << ((double)p.total_size) / (1024 * 1024) << " MB in "
                  << gpu_shared_pool.size() << " blocks\n";
        hila::out << "   # of allocations " << p.n_allocs << "  real allocs "
                  << std::setprecision(2) << ((double)p.n_true_allocs) / p.n_allocs * 100 << "%\n";
        hila::out << "   Average block list search " << (double)p.blocklist_avg_search / p.n_allocs
                  << " steps\n\n";
    }
}
#endif // GPU_SHMEM

#endif // GPU_MEMORY_POOL
#endif // !HILAPP
