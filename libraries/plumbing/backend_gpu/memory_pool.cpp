///////////////////////////////////////////
/// gpu_malloc.cpp - simple list-based alloc program for cuda/hip

#include "plumbing/defs.h"
#include "plumbing/lattice.h"
#include "plumbing/field.h"
#include "plumbing/backend_gpu/defs.h"
#include <list>
#include <iomanip>

// no real need for HILAPP to go through here 
#if defined(GPU_MEMORY_POOL) && !defined(HILAPP)

#if defined(HIP)
#define gpuMallocDirect(a, b) GPU_CHECK(hipMalloc(a, b))
#define gpuFreeDirect(a) GPU_CHECK(hipFree(a))
#elif defined(CUDA)
#define gpuMallocDirect(a, b) GPU_CHECK(cudaMalloc(a, b))
#define gpuFreeDirect(a) GPU_CHECK(cudaFree(a))
#else
static_assert(0 && "HIP or CUDA must be defined");
#endif


// keep relatively large min allocation
#define MIN_ALLOC_SIZE 128

struct allocation {
    void *ptr;
    size_t size;
};

static size_t total_size = 0;
static size_t n_allocs = 0;
static size_t n_true_allocs = 0;
static double free_list_avg_size = 0;
static double free_list_avg_search = 0;

static std::list<allocation> free_list = {};
static std::list<allocation> in_use_list = {};

void gpu_memory_pool_alloc(void **p, size_t req_size) {

    if (req_size < MIN_ALLOC_SIZE) {
        req_size = MIN_ALLOC_SIZE;
    }

    n_allocs++;
    free_list_avg_size += free_list.size();

    // do we have free stuff?  Simple linear search - list should not be too large
    bool found_match = false;
    auto ptr = free_list.begin();
    int steps = 0;
    for (auto it = free_list.begin(); it != free_list.end(); it++) {
        steps++;
        if (it->size == req_size) {
            found_match = true;
            ptr = it;
            break; // perfect match, that's it
        }

        // allow allocated blocks at most twice larger
        if (it->size > req_size && it->size < 2 * req_size) {
            if (!found_match || ptr->size > it->size) {
                ptr = it;
            }
            found_match = true;
        }
    }

    free_list_avg_search += steps;

    // got it, move to in_use_list to the beginning (faster to find)
    if (found_match) {
        *p = ptr->ptr;
        in_use_list.splice(in_use_list.begin(), free_list, ptr);

    } else {

        // did not find free memory - allocate
        // alloc failure caught by gpuMalloc
        allocation a;
        gpuMallocDirect(&(a.ptr), req_size);
        a.size = req_size;
        in_use_list.push_front(a);
        *p = a.ptr;

        n_true_allocs++;
        total_size += req_size;
    }
}

void gpu_memory_pool_free(void *ptr) {

    // search the list for the memory block
    for (auto it = in_use_list.begin(); it != in_use_list.end(); it++) {
        if (it->ptr == ptr) {
            // found the allocation, move to free list to the beginning

            free_list.splice(free_list.begin(), in_use_list, it);

            return;
        }
    }

    // did not find!  serious error, quit
    hila::out << "Memory free error - unknown pointer  " << ptr << '\n';
    hila::terminate(1);
}

/// Release free memory to the system - avoids extra allocations
void gpu_memory_pool_purge() {

    for (auto it = free_list.begin(); it != free_list.end(); it++) {
        gpuFreeDirect(it->ptr);

        total_size -= it->size;
    }

    free_list.clear();
}

void gpu_memory_pool_report() {
    if (hila::myrank() == 0) {
        hila::out << "\nGPU Memory pool statistics from node 0:\n";
        hila::out << "   Total pool size " << ((double)total_size) / (1024 * 1024)
                     << " MB\n";
        hila::out << "   # of allocations " << n_allocs << "  real allocs "
                     << std::setprecision(2) << ((double)n_true_allocs) / n_allocs * 100
                     << "%\n";
        hila::out << "   Average free list search "
                     << free_list_avg_search / n_allocs << " steps\n";
        hila::out << "   Average free list size " << free_list_avg_size / n_allocs
                     << " items\n\n";
    }
}

#endif // GPU_MEMORY_POOL
