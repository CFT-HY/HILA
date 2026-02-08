/** @file field_storage.h */

#ifndef FIELD_STORAGEH
#define FIELD_STORAGEH

#include "plumbing/defs.h"
#include "plumbing/coordinates.h"
#include "plumbing/field.h"
#include "plumbing/backend_vector/vector_types.h"

#include "plumbing/has_unary_minus.h"


template <typename T>
class field_struct;

#if defined(CUDA) || defined(HIP)
#define DEVICE __device__
#else
#define DEVICE
#endif

/**
 * @brief The field_storage struct contains minimal information for using
 * the field in a loop. It is communicated to CUDA kernels and other
 * accelerator functions.
 */
template <typename T>
class field_storage {
  public:
    // The actual data content of the class: a pointer to the field and a
    // list of neighbour pointers.
    T *RESTRICT fieldbuf = nullptr;
    
    const unsigned *RESTRICT neighbours[NDIRS];

    void allocate_field(const Lattice lattice);
    void free_field();
    void copy_field(const Lattice lattice, const field_storage<T> &source);

#ifndef VECTORIZED
    // Get an element in a loop
    DEVICE inline auto get(const unsigned i, const unsigned field_alloc_size) const;

    // template <typename A>
    DEVICE inline void set(const T &value, const unsigned i, const unsigned field_alloc_size);

    // Get a single element outside loops
    auto get_element(const unsigned i, const Lattice lattice) const;
    template <typename A>
    void set_element(A &value, const unsigned i, const Lattice lattice);

#else
    inline T get_element(const unsigned i) const;

    inline void set_element(const T &value, const unsigned i);

    // in vector code, write only 1 element to field at site index idx
    template <typename vecT>
    inline void set_vector(const vecT &val, const unsigned idx);

    template <typename vecT>
    inline vecT get_vector(const unsigned idx) const;

    void gather_comm_vectors(
        T *RESTRICT buffer, const lattice_struct::comm_node_struct &to_node, Parity par,
        const vectorized_lattice_struct<hila::vector_info<T>::vector_size> *RESTRICT vlat,
        bool antiperiodic) const;

    void place_recv_elements(
        const T *RESTRICT buffer, Direction d, Parity par,
        const vectorized_lattice_struct<hila::vector_info<T>::vector_size> *RESTRICT vlat) const;

#endif

    void gather_comm_elements(T *RESTRICT buffer, const lattice_struct::comm_node_struct &to_node,
                              Parity par, const Lattice lattice, bool antiperiodic) const;

    void gather_elements(T *RESTRICT buffer, const unsigned *RESTRICT index_list, int n,
                         const Lattice lattice) const;

    void gather_elements_negated(T *RESTRICT buffer, const unsigned *RESTRICT index_list, int n,
                                 const Lattice lattice) const;

    /// Place boundary elements from neighbour
    void place_comm_elements(Direction d, Parity par, T *RESTRICT buffer,
                             const lattice_struct::comm_node_struct &from_node,
                             const Lattice lattice);
    void place_elements(T *RESTRICT buffer, const unsigned *RESTRICT index_list, int n,
                        const Lattice lattice);
    /// Place boundary elements from local lattice (used in vectorized version)
    void set_local_boundary_elements(Direction dir, Parity par, const Lattice lattice,
                                     bool antiperiodic);

    // Allocate buffers for mpi communication
    T *allocate_mpi_buffer(unsigned n);
    void free_mpi_buffer(T *buffer);

    T *RESTRICT get_buffer() {
        return static_cast<T *>(fieldbuf);
    }
};

/*
Import backend
*/

#ifdef VECTORIZED

#include "plumbing/backend_vector/field_storage_backend.h"

#elif defined(CUDA) || defined(HIP)

#include "plumbing/backend_gpu/field_storage_backend.h"

#elif defined(VANILLA)

#include "plumbing/backend_cpu/field_storage_backend.h"

#else
Something must be defined !

#endif

#endif // field storage