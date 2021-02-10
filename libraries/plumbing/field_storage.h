#ifndef FIELD_STORAGEH
#define FIELD_STORAGEH

#include "plumbing/globals.h"

#include "plumbing/defs.h"
#include "plumbing/field.h"
#include "plumbing/backend_vector/vector_types.h"

#include "plumbing/has_unary_minus.h"

// Forward declare the field struct
template <typename T> class field_struct;

#ifdef __CUDACC__
#define DEVICE __device__
#else
#define DEVICE
#endif

////////////////////////////////////////////////////////////////////////
/// The field_storage struct contains minimal information for using
/// the field in a loop. It is communicated to CUDA kernels and other
/// accelerator functions.
////////////////////////////////////////////////////////////////////////
template <typename T> class field_storage {
  public:
    // The actual data content of the class: a pointer to the field and a
    // list of neighbour pointers.
    T *RESTRICT fieldbuf = nullptr;
    const unsigned *RESTRICT neighbours[NDIRS];

    void allocate_field(lattice_struct *lattice);
    void free_field();

#ifndef VECTORIZED
    // Get an element in a loop
    DEVICE inline auto get(const unsigned i, const unsigned field_alloc_size) const;

    template <typename A>
    DEVICE inline void set(const A &value, const unsigned i,
                           const unsigned field_alloc_size);

    // Get a single element outside loops
    inline auto get_element(const unsigned i,
                            const lattice_struct *RESTRICT lattice) const;
    template <typename A>
    inline void set_element(A &value, const unsigned i,
                            const lattice_struct *RESTRICT lattice);

#else
    inline T get_element(const unsigned i) const;

    inline void set_element(const T &value, const unsigned i);

    // in vector code, write only 1 element to field at site index idx
    template <typename vecT> inline void set_vector(const vecT &val, const unsigned idx);

    template <typename vecT> inline vecT get_vector(const unsigned idx) const;

    void gather_comm_vectors(
        T *RESTRICT buffer, const lattice_struct::comm_node_struct &to_node, parity par,
        const vectorized_lattice_struct<vector_info<T>::vector_size> *RESTRICT vlat,
        bool antiperiodic) const;

    void place_recv_elements(const T *RESTRICT buffer, direction d, parity par,
                             const vectorized_lattice_struct<vector_info<T>::vector_size>
                                 *RESTRICT vlat) const;

#endif

    void gather_comm_elements(T *RESTRICT buffer,
                              const lattice_struct::comm_node_struct &to_node, parity par,
                              const lattice_struct *RESTRICT lattice,
                              bool antiperiodic) const;
    void gather_elements(T *RESTRICT buffer, const unsigned *RESTRICT index_list, int n,
                         const lattice_struct *RESTRICT lattice) const;

    void gather_elements_negated(T *RESTRICT buffer, const unsigned *RESTRICT index_list,
                                 int n, const lattice_struct *RESTRICT lattice) const;

    /// Place boundary elements from neighbour
    void place_comm_elements(direction d, parity par, T *RESTRICT buffer,
                             const lattice_struct::comm_node_struct &from_node,
                             const lattice_struct *RESTRICT lattice);
    void place_elements(T *RESTRICT buffer, const unsigned *RESTRICT index_list, int n,
                        const lattice_struct *RESTRICT lattice);
    /// Place boundary elements from local lattice (used in vectorized version)
    void set_local_boundary_elements(direction dir, parity par,
                                     const lattice_struct *RESTRICT lattice,
                                     bool antiperiodic);

    // Allocate buffers for mpi communication
    T *allocate_mpi_buffer(unsigned n);
    void free_mpi_buffer(T *buffer);

    T *RESTRICT get_buffer() { return static_cast<T *>(fieldbuf); }
};

/*
Import backend
*/

#ifdef VECTORIZED

#include "plumbing/backend_vector/field_storage_backend.h"

#elif defined(CUDA)

#include "plumbing/backend_cuda/field_storage_backend.h"

#elif defined(VANILLA)

#include "plumbing/backend_cpu/field_storage_backend.h"

#elif
Something must be defined !

#endif

#endif // field storage