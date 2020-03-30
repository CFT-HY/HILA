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
    auto get(const int i, const int field_alloc_size) const;

    template<typename A>
    #pragma transformer loop_function
    inline void set(const A &value, const int i, const int field_alloc_size);

    void gather_comm_elements(char * buffer, lattice_struct::comm_node_struct to_node, parity par, lattice_struct * lattice) const;
    void gather_elements(char * buffer, std::vector<unsigned> index_list, lattice_struct * lattice) const;
    /// Place boundary elements from neighbour
    void place_comm_elements(char * buffer, lattice_struct::comm_node_struct from_node, parity par, lattice_struct * lattice);
    void place_elements(char * buffer, std::vector<unsigned> index_list, lattice_struct * lattice);
    /// Place boundary elements from local lattice (used in vectorized version)
    void set_local_boundary_elements(direction dir, parity par, lattice_struct * lattice);

    // HACK:
    #if defined(VANILLA)
    T * RESTRICT get_buffer() {
      return static_cast<T*>(fieldbuf);
    }
    #endif
};


template<typename T>
void field_storage<T>::gather_comm_elements(char * buffer, lattice_struct::comm_node_struct to_node, parity par, lattice_struct * lattice) const {
  std::vector<unsigned> index_list = to_node.get_site_list(par);
  gather_elements(buffer, index_list, lattice);
}

template<typename T>
void field_storage<T>::place_comm_elements(char * buffer, lattice_struct::comm_node_struct from_node, parity par, lattice_struct * lattice){
  std::vector<unsigned> index_list(from_node.n_sites(par));
  for (int j=0; j<from_node.n_sites(par); j++) {
    index_list[j] = from_node.offset(par)+j;
  }
  place_elements(buffer, index_list, lattice);
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