#ifndef HILA_MEMORY_POOL_H
#define HILA_MEMORY_POOL_H

///////////////////////////////////////////
/// simple list-based alloc program for cuda/hip

#include "plumbing/defs.h"
#include <iomanip>


namespace hila {

struct memory_pool_status {
    size_t total_size = 0;
    size_t n_allocs = 0;
    size_t free_allocs = 0;
    size_t in_use_allocs = 0;
    size_t n_true_allocs = 0;
    size_t first_free_block = 0;
    size_t blocklist_avg_search = 0;
};

class memory_pool {
  private:
    struct allocation {
        void *ptr;
        size_t size;
        bool in_use;
    };

    memory_pool_status p;

    std::vector<allocation> blocklist = {};

  public:
    void *alloc(size_t req_size);
    void free(void *ptr);
    void purge();
    const memory_pool_status &status() const {
        return p;
    }
    size_t size() const {
        return blocklist.size();
    }

}; // class memory_pool
} // namespace hila

#endif