///////////////////////////////////////////
/// gpu_malloc.cpp - simple list-based alloc program for cuda/hip

#include "plumbing/defs.h"
#include "plumbing/lattice.h"
#include "plumbing/field.h"
#include "plumbing/backend_gpu/defs.h"
#include <list>
#include <iomanip>

///////////////////////////////////////////////////////////////////////
/// GPU memory manager
/// Allocates a slab of memory, which it then splits out in blocks as requested.
/// On free, merges adjoining free blocks to a larger block
/// Optimized for relatively large allocations and roughly "fifo"-type
/// alloc/free cycle.
///
/// Allocates first GPU_MEMORY_POOL_FRACTION part of the GPU memory
/// If this becomes full, tries to allocate half of the remaining memory
/// NOTE: do not allocate all of memory on GPU, for example FFTs need
/// own memory
///
/// Implements 3 doubly linked lists:
///    memory blocks allocated (ordered, all blocks are here)
///    free memory blocks
///    in use memory blocks
///////////////////////////////////////////////////////////////////////

#if defined(GPU_MEMORY_POOL)

#ifndef GPU_MEMORY_POOL_FRACTION
#define GPU_MEMORY_POOL_FRACTION 0.5
#endif

#if !defined(CUDA) && !defined(HIP)
static_assert(0 && "HIP or CUDA must be defined");
#endif

// keep relatively large min allocation
#define ALLOC_ALIGNMENT 256

// Use doubly 2-directional linked list to keep track
struct block {
    void *ptr;
    block *next, *prev;
    block *up, *down;
    size_t size;
    bool is_free;
};

static size_t max_used_size = 0;
static size_t current_used_size = 0;
static size_t total_pool_size = 0;
static size_t n_allocs = 0;
static size_t free_list_size;
static size_t gpu_total_memory = 0;

static double free_list_avg_size = 0;
static double free_list_avg_search = 0;

static block *memory_blocks; // keeps track of all of memory
static block *free_blocks;
static block *in_use_blocks;

static block *unused_blocks;

//////////////////////////////////////////////////////////////////////////
// Memory block descriptor manager
// make some memory alloc tracking blocks - last.up == nullptr
block *alloc_more_block_descriptors() {
    block *p;
    int n_blocks = 1000;

    p = (block *)memalloc(n_blocks * sizeof(block));

    for (int i = 0; i < n_blocks - 1; i++) {
        p[i].up = &p[i + 1];
    }

    p[n_blocks - 1].up = nullptr;

    return p;
}

block *get_block_descriptor() {
    if (unused_blocks == nullptr) {
        unused_blocks = alloc_more_block_descriptors();
    }
    block *ret = unused_blocks;
    unused_blocks = unused_blocks->up;
    ret->ptr = nullptr;
    ret->up = ret->down = ret->next = ret->prev = nullptr;
    return ret;
}

void release_block_descriptor(block *p) {
    p->up = unused_blocks;
    unused_blocks = p;
}

/////////////////////////////////////////////////////////////////////////
// Alloc memory slabs
/////////////////////////////////////////////////////////////////////////
void *gpu_memory_allocate(size_t m_alloc) {
    const int kb = 1024;
    const int mb = kb * kb;

    // ensure size is multiple of alignment
    if (m_alloc % ALLOC_ALIGNMENT != 0)
        m_alloc = m_alloc - m_alloc % ALLOC_ALIGNMENT + ALLOC_ALIGNMENT;

    double fraction = (double)m_alloc / gpu_total_memory;

    output0 << "GPU memory: allocating " << m_alloc / mb << " MB out of total "
            << gpu_total_memory / mb << "(" << (int)(fraction * 100) << "%)\n";
    total_pool_size += m_alloc;

    void *b = nullptr;

#ifndef HILAPP
#if defined(CUDA)
    GPU_CHECK(cudaMalloc(&b, m_alloc));
#elif defined(HIP)
    GPU_CHECK(hipMalloc(&b, m_alloc));
#endif
#endif

    return b;
}

/////////////////////////////////////////////////////////////////////////
// Init memory; allocate the slab
/////////////////////////////////////////////////////////////////////////
void gpu_memory_pool_init() {

#ifndef HILAPP
#if defined(CUDA)
    cudaDeviceProp props;
    int my_device;
    GPU_CHECK(cudaGetDevice(&my_device));
    GPU_CHECK(cudaGetDeviceProperties(&props, my_device));
#elif defined(HIP)
    hipDeviceProp_t props;
    int my_device;
    GPU_CHECK(hipGetDevice(&my_device));
    GPU_CHECK(hipGetDeviceProperties(&props, my_device));
#endif

    gpu_total_memory = props.totalGlobalMem;
#endif // HILAPP

    size_t m_alloc = gpu_total_memory * GPU_MEMORY_POOL_FRACTION;
    // ensure size is multiple of alignment
    m_alloc = m_alloc - m_alloc % ALLOC_ALIGNMENT + ALLOC_ALIGNMENT;

    block *b = get_block_descriptor();
    memory_blocks = b;

    b->ptr = gpu_memory_allocate(m_alloc);

    b->up = b->down = b->next = b->prev = nullptr;
    b->size = m_alloc;
    b->is_free = true;

    // one huge block free
    free_blocks = b;
    in_use_blocks = nullptr;
    free_list_size = 1;
    max_used_size = 0;
    current_used_size = 0;
}

//////////////////////////////////////////////////////////////////////
// Manage free and in-use lists: remove_from_list and
// insert_to_list_head.  New list member is at the head, anticipating
// reuse soon

void remove_from_list(block *p, block **head) {
    if (p->next != nullptr)
        p->next->prev = p->prev;
    if (p->prev != nullptr)
        p->prev->next = p->next;
    else
        *head = p->next;
}

void insert_to_list_head(block *p, block **head) {
    if (*head != nullptr)
        (*head)->prev = p;
    p->next = *head;
    p->prev = nullptr;
    *head = p;
}

///////////////////////////////////////////////////////////////////////
//  Manage active block lists:

// Merge block with free block below.  If orig block was free, remove
// from free list
void merge_block_down_free(block *p) {
    block *pdown = p->down;
    pdown->size += p->size;
    pdown->up = p->up;
    if (p->up != nullptr)
        p->up->down = pdown;

    if (p->is_free) {
        // remove from free list
        remove_from_list(p, &free_blocks);
        free_list_size--;
    }

    release_block_descriptor(p);
}

//  Merge block with free block above.
void merge_block_up_free(block *p) {
    block *pup = p->up;
    pup->size += p->size;
    pup->ptr = p->ptr; // set the ptr to base
    pup->down = p->down;
    if (p->down != nullptr)
        p->down->up = pup;
    else
        memory_blocks = pup;

    if (p->is_free) {
        remove_from_list(p, &free_blocks);
        free_list_size--;
    }

    release_block_descriptor(p);
}

// split free block, insert new block above p and return it
// return the upper slice of the free block
block *split_free_block(block *p, size_t req_size) {
    block *b = get_block_descriptor();
    b->up = p->up;
    b->down = p;
    p->up = b;
    if (b->up != nullptr)
        b->up->down = b;

    // set data pointer as appropriate
    b->ptr = static_cast<char *>(p->ptr) + p->size - req_size;
    b->size = req_size;
    p->size -= req_size;

    b->is_free = false;
    // Free list pointers remain OK (p is there)

    return b;
}

// Make memory block free, insert to free list
void mark_block_free(block *p) {
    p->is_free = true;
    insert_to_list_head(p, &free_blocks);
    free_list_size++;
}

// Add block to the tail of block list
void add_block_to_top(block *p) {
    p->up = nullptr;
    if (memory_blocks == nullptr) {
        memory_blocks = p;
        p->down = nullptr;
    } else {
        block *b;
        for (b = memory_blocks; b->up != nullptr; b = b->up)
            ;
        b->up = p;
        p->down = b;
    }
}

/////////////////////////////////////////////////////////////////////
// Memory allocator; similar interface to cudaMalloc
/////////////////////////////////////////////////////////////////////

void gpu_memory_pool_alloc(void **ret, size_t req_size) {

    static bool initialized = false;

    if (!initialized) {
        initialized = true;
        gpu_memory_pool_init();
    }

    // make req_size to be multiple of alignment
    size_t align_mod = req_size % ALLOC_ALIGNMENT;
    if (align_mod > 0)
        req_size = req_size - align_mod + ALLOC_ALIGNMENT;

    // output0 << "REQ SIZE " << req_size << '\n';

    n_allocs++;
    free_list_avg_size += free_list_size;

    // do we have free stuff?  Simple linear search - list should not be too large
    bool found_match = false;
    int steps = 0;
    block *ptr = free_blocks;
    for (block *p = free_blocks; p != nullptr; p = p->next) {
        steps++;
        if (p->size == req_size) {
            found_match = true; // perfect match, use it
            ptr = p;
            break;
        }

        if (p->size > req_size) {
            // find smallest free block which is OK
            if (!found_match || ptr->size > p->size) {
                ptr = p;
            }
            found_match = true;
        }
    }

    free_list_avg_search += steps;

    // got it, split a piece out of it
    if (found_match) {
        if (ptr->size > req_size) {
            ptr = split_free_block(ptr, req_size);
        } else {
            // now ptr->size == req_size
            // rm from free list
            remove_from_list(ptr, &free_blocks);
            free_list_size--;
            ptr->is_free = false;
        }

        // move to in-use list
        insert_to_list_head(ptr, &in_use_blocks);

        current_used_size += req_size;
        if (current_used_size > max_used_size)
            max_used_size = current_used_size;

        *ret = ptr->ptr;
        return;

    } else {
        // try to allocate more?
        if (total_pool_size < 0.9 * gpu_total_memory) {
            size_t m_alloc = 0.5 * (gpu_total_memory - total_pool_size);
            m_alloc = (m_alloc > req_size) ? m_alloc : req_size;
            // leave 5% of total memory
            if (m_alloc + total_pool_size < 0.95 * gpu_total_memory) {
                // put an "empty" block as a separator (non-mergeable)
                block *p = get_block_descriptor();
                p->size = 0;
                p->is_free = false;
                add_block_to_top(p);

                // and new memory block
                p = get_block_descriptor();
                p->ptr = gpu_memory_allocate(m_alloc);
                p->size = m_alloc;
                add_block_to_top(p);
                mark_block_free(p);

                gpu_memory_pool_alloc(ret, req_size);
                return;
            }
        }
    }

    hila::output << "Out of memory in GPU pool, request size " << req_size << '\n';
    hila::terminate(0);
}

//////////////////////////////////////////////////////////////////////
// And release memory.  Pointer must be exactly the same!
//////////////////////////////////////////////////////////////////////

void gpu_memory_pool_free(void *ptr) {

    // search the list for the memory block
    for (block *f = in_use_blocks; f != nullptr; f = f->next) {
        if (f->ptr == ptr) {
            // found the allocation, remove from in_use
            remove_from_list(f, &in_use_blocks);

            current_used_size -= f->size;

            // Are neighbour blocks also free?
            block *down = f->down;
            block *up = f->up;
            if (down != nullptr && down->is_free) {
                merge_block_down_free(f);
                if (up != nullptr && up->is_free) {
                    merge_block_down_free(up);
                }
            } else if (up != nullptr && up->is_free) {
                merge_block_up_free(f);
            } else {
                // no merging now
                mark_block_free(f);
            }

            return;
        }
    }

    // did not find!  serious error, quit
    hila::output << "Memory free error - unknown pointer  " << ptr << '\n';
    hila::terminate(1);
}

/// Release free memory to the system - avoids extra allocations
void gpu_memory_pool_purge() {}

void gpu_memory_pool_report() {
    if (hila::myrank() == 0) {
        hila::output << "\nGPU Memory pool statistics from node 0:\n";
        hila::output << "   Total pool size "
                     << ((double)total_pool_size) / (1024 * 1024) << " MB\n";
        hila::output << "   # of allocations " << n_allocs << '\n';
        hila::output << "   Average free list search "
                     << free_list_avg_search / n_allocs << " steps\n";
        hila::output << "   Average free list size " << free_list_avg_size / n_allocs
                     << " items\n";
        hila::output << "   Maximum memory use " << max_used_size / (1024 * 1024)
                     << " MB\n\n";
    }
}

#endif // GPU_MEMORY_POOL
