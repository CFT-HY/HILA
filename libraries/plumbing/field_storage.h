#ifndef FIELD_STORAGEH
#define FIELD_STORAGEH


#include "../plumbing/globals.h"

#include "../plumbing/defs.h"
#include "../plumbing/field.h"
#include "../plumbing/backend_vector/vector_types.h"

#include "../plumbing/has_unary_minus.h"


// Pointer to field data and accessors. Only this is passed to the
// CUDA kernels and other accelerators and it only contains a minimal
// amount of data.

template <typename T>
class field_struct;

template <typename T>
class field_storage {
  public:

    T * RESTRICT fieldbuf = nullptr;

    void allocate_field( lattice_struct * lattice );
    void free_field();

#ifndef VECTORIZED
    #pragma hila loop_function
    inline auto get(const int i, const int field_alloc_size) const;

    template<typename A>
    #pragma hila loop_function
    inline void set(const A &value, const int i, const int field_alloc_size);
#else
    #pragma hila loop_function
    inline T get_element(const int i) const;

    #pragma hila loop_function
    inline void set_element(const T &value, const int i);

    // in vector code, write only 1 element to field at site index idx
    template <typename vecT>
    #pragma hila loop_function
    inline void set_vector( const vecT & val, const int idx );

    template <typename vecT>
    #pragma hila loop_function
    inline vecT get_vector( const int idx ) const;

    void gather_comm_vectors( T * RESTRICT buffer, const lattice_struct::comm_node_struct & to_node, 
      parity par, const vectorized_lattice_struct<vector_info<T>::vector_size> * RESTRICT vlat, bool antiperiodic) const;

    void place_recv_elements(const T * RESTRICT buffer, direction d, parity par,
                             const vectorized_lattice_struct<vector_info<T>::vector_size> * RESTRICT vlat) const;

#endif


    void gather_comm_elements( T * RESTRICT buffer, const lattice_struct::comm_node_struct & to_node, 
                               parity par, const lattice_struct * RESTRICT lattice) const;
    void gather_elements( T * RESTRICT buffer, const unsigned * RESTRICT index_list, int n, 
                          const lattice_struct * RESTRICT lattice) const;



    void gather_elements_negated( T * RESTRICT buffer, const unsigned * RESTRICT index_list, int n, 
                                  const lattice_struct * RESTRICT lattice) const;


    /// Place boundary elements from neighbour
    void place_comm_elements( direction d, parity par, T * RESTRICT buffer, 
                              const lattice_struct::comm_node_struct & from_node,
                              const lattice_struct * RESTRICT lattice);
    void place_elements( T * RESTRICT buffer, const unsigned * RESTRICT index_list, int n,
                         const lattice_struct * RESTRICT lattice);
    /// Place boundary elements from local lattice (used in vectorized version)
#ifndef VECTORIZED
    void set_local_boundary_elements(direction dir, parity par, lattice_struct * RESTRICT lattice);
#else
    void set_local_boundary_elements(direction dir, parity par, const lattice_struct * RESTRICT lattice, bool antiperiodic);
#endif

    T * RESTRICT get_buffer() {
      return static_cast<T*>(fieldbuf);
    }
};


template<typename T>
void field_storage<T>::gather_comm_elements(T * RESTRICT buffer, 
                                            const lattice_struct::comm_node_struct & to_node, 
                                            parity par, const lattice_struct * RESTRICT lattice) const {
  int n;
  const unsigned * index_list = to_node.get_sitelist(par,n);

  gather_elements(buffer, index_list, n, lattice);
}

#ifndef VANILLA

template<typename T>
void field_storage<T>::place_comm_elements(direction d, parity par, T * RESTRICT buffer, 
                                           const lattice_struct::comm_node_struct & from_node, 
                                           const lattice_struct * RESTRICT lattice) {
  int n;
  const unsigned * index_list = from_node.get_sitelist(par,n);
  place_elements(buffer, index_list, n, lattice);
}

#endif

/*
Import backend
*/

#ifdef VECTORIZED

#include "backend_vector/field_storage_backend.h"

#elif defined(CUDA)

#include "backend_cuda/field_storage_backend.h"

#elif defined(VANILLA)

#include "backend_cpu/field_storage_backend.h"

#elif 
Something must be defined!

#endif

#endif //field storage