#ifndef FIELD_STORAGEH
#define FIELD_STORAGEH

// Pointer to field data and accessors. Only this is passed to the
// CUDA kernels and other accelerators and it only contains a minimal
// amount of data.

template <typename T>
class field_storage {
  public:

    void * fieldbuf = NULL;

    void allocate_field( lattice_struct * lattice );
    void free_field();
    #pragma transformer loop_function
    auto get(const int i, const int field_alloc_size) const;
    #pragma transformer loop_function
    template<typename A>
    inline void set(const A &value, const int i, const int field_alloc_size);

    void gather_comm_elements(char * buffer, lattice_struct::comm_node_struct to_node, parity par, lattice_struct * lattice) const;
    /// Place boundary elements from neighbour
    void place_comm_elements(char * buffer, lattice_struct::comm_node_struct from_node, parity par, lattice_struct * lattice);
    /// Place boundary elements from local lattice (used in vectorized version)
    void set_local_boundary_elements(parity par, lattice_struct * lattice);

};


/*
Import backend
*/

#ifdef VECTORIZED

#include "../plumbing/backend_vector/field_storage_backend.h"

#elif defined(CUDA)

#include "../plumbing/backend_cuda/field_storage_backend.h"

#elif defined(VANILLA)

#include "../plumbing/backend_cpu/field_storage_backend.h"

#endif

#endif //field storage