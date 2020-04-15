#ifndef FIELD_STORAGEH
#define FIELD_STORAGEH

// Pointer to field data and accessors. Only this is passed to the
// CUDA kernels and other accelerators and it only contains a minimal
// amount of data.

template <typename T>
class field_storage {
  public:

    void * RESTRICT fieldbuf = nullptr;

    void allocate_field( lattice_struct * lattice );
    void free_field();

    #pragma transformer loop_function
    inline auto get(const int i, const int field_alloc_size) const;

    template<typename A>
    #pragma transformer loop_function
    inline void set(const A &value, const int i, const int field_alloc_size);

    void gather_comm_elements( char * RESTRICT buffer, const lattice_struct::comm_node_struct & to_node, 
                               parity par, const lattice_struct * RESTRICT lattice) const;
    void gather_elements( char * RESTRICT buffer, const unsigned * RESTRICT index_list, int n, 
                          const lattice_struct * RESTRICT lattice) const;
    /// Place boundary elements from neighbour
    void place_comm_elements( char * RESTRICT buffer, const lattice_struct::comm_node_struct & from_node,
                              parity par, const lattice_struct * RESTRICT lattice);
    void place_elements( char * RESTRICT buffer, const unsigned * RESTRICT index_list, int n,
                         const lattice_struct * RESTRICT lattice);
    /// Place boundary elements from local lattice (used in vectorized version)
    void set_local_boundary_elements(direction dir, parity par, lattice_struct * RESTRICT lattice);

    // HACK:
    #if defined(VANILLA)
    T * RESTRICT get_buffer() {
      return static_cast<T*>(fieldbuf);
    }
    #endif
};


template<typename T>
void field_storage<T>::gather_comm_elements(char * RESTRICT buffer, 
                                            const lattice_struct::comm_node_struct & to_node, 
                                            parity par, const lattice_struct * RESTRICT lattice) const {
  int n;
  const unsigned * index_list = to_node.get_sitelist(par,n);
  gather_elements(buffer, index_list, n, lattice);
}

template<typename T>
void field_storage<T>::place_comm_elements(char * RESTRICT buffer, 
                                           const lattice_struct::comm_node_struct & from_node, 
                                           parity par, const lattice_struct * RESTRICT lattice) {
  int n;
  const unsigned * index_list = from_node.get_sitelist(par,n);
  place_elements(buffer, index_list, n, lattice);
}



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