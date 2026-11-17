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
                              Parity par, const Lattice lattice, bool antiperiodic
#if defined(CUDA) || defined(HIP)
                              ,
                              gpuStream_t stream
#endif
    ) const;

    void gather_elements(T *RESTRICT buffer, const unsigned *RESTRICT index_list, int n,
                         const Lattice lattice) const;

    void gather_elements_negated(T *RESTRICT buffer, const unsigned *RESTRICT index_list, int n,
                                 const Lattice lattice) const;

    /// Place boundary elements from neighbour
    void place_comm_elements(Direction d, Parity par, T *RESTRICT buffer,
                             const lattice_struct::comm_node_struct &from_node,
                             const Lattice lattice
#if defined(CUDA) || defined(HIP)
                             ,
                             gpuStream_t stream
#endif
    );
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


// these are here for now, used to map the input args of kernelized code to the indexes in user
// written code.
namespace hila {

/**
 * @brief Select the storage of a field by an index value. Used by hilapp in GPU kernels for
 * a field selected with an index that depends on a variable defined inside the site loop, e.g.
 * U[mu] inside foralldir(mu).
 *
 * @details hilapp passes the storages of all fields the reference can be to the kernel:
 * s0, ss... are the storages for the possible index values V0, Vs... found by its analysis,
 * e.g. select_storage<0, 1, 2, 3>(mu, U0, U1, U2, U3). The index is the same for all
 * threads. The last storage is returned without a check, the analysis guarantees that the index
 * has one of the values.
 */
template <int64_t V0, int64_t... Vs, typename S, typename... Ss>
DEVICE inline const S &select_storage(const int64_t idx, const S &s0, const Ss &...ss) {
    static_assert(sizeof...(Vs) == sizeof...(Ss), "a storage for each index value needed");
    if constexpr (sizeof...(Vs) == 0) {
        (void)idx;
        return s0;
    } else {
        return (idx == V0) ? s0 : select_storage<Vs...>(idx, ss...);
    }
}

/**
 * @brief The position of the value `v` in the list of values V0, Vs... Used by hilapp
 * (FETCH_FIELD_REFS_LOOP_LOCAL) to find the element of a field with a loop local index read into an
 * array at the start of the loop, e.g. value_position<0, 1, 2, 3>(mu) for U[mu][X].
 *
 * @details The last position is returned without a check, the analysis of hilapp guarantees that
 * v is one of the values. A constant if the loop over the index is unrolled.
 */
template <int64_t V0, int64_t... Vs>
DEVICE constexpr inline int value_position(const int64_t v) {
    if constexpr (sizeof...(Vs) == 0) {
        (void)v;
        return 0;
    } else {
        return (v == V0) ? 0 : 1 + value_position<Vs...>(v);
    }
}

/**
 * @brief The position of the pair (i, d) in the list of pairs (I0, D0), Rest... given as
 * I0, D0, I1, D1, ... Like value_position(), used for (index value, direction) pairs of the
 * references to a field with a loop local index and direction, e.g.
 * pair_position<0, 7, 1, 6>(mu, -mu) for U[mu][X - mu].
 */
template <int64_t I0, int64_t D0, int64_t... Rest>
DEVICE constexpr inline int pair_position(const int64_t i, const int64_t d) {
    static_assert(sizeof...(Rest) % 2 == 0, "pairs of (index value, direction) needed");
    if constexpr (sizeof...(Rest) == 0) {
        (void)i;
        (void)d;
        return 0;
    } else {
        return (i == I0 && d == D0) ? 0 : 1 + pair_position<Rest...>(i, d);
    }
}
} // namespace hila

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