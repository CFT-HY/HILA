#ifndef VECTOR_BACKEND
#define VECTOR_BACKEND

#include "../defs.h"
#include "../lattice.h"
#include "../field_storage.h"
#include "vector_types.h"
#include "../coordinates.h"
#include "defs.h"

/// Replaces basetypes with vectors in a given templated class

/// First base definition for replace_type, which recursively looks for the
/// base type and replaces it in the end
/// General template, never matched
template <typename A, int vector_size, class Enable = void>
struct vectorize_struct {};

/// A is a basic type, so just return the matching vector type
template <typename A, int vector_size>
struct vectorize_struct<A, vector_size, typename std::enable_if_t<hila::is_arithmetic<A>::value>> {
    using type = typename hila::vector_base_type<A, vector_size>::type;
};

// B is a templated class, so construct a vectorized type
template <template <typename B> class C, typename B, int vector_size>
struct vectorize_struct<C<B>, vector_size> {
    using vectorized_B = typename vectorize_struct<B, vector_size>::type;
    using type = C<vectorized_B>;
};

template <template <int a, typename B> class C, int a, typename B, int vector_size>
struct vectorize_struct<C<a, B>, vector_size> {
    using vectorized_B = typename vectorize_struct<B, vector_size>::type;
    using type = C<a, vectorized_B>;
};

template <template <int a, int b, typename B> class C, int a, int b, typename B, int vector_size>
struct vectorize_struct<C<a, b, B>, vector_size> {
    using vectorized_B = typename vectorize_struct<B, vector_size>::type;
    using type = C<a, b, vectorized_B>;
};

/// Match coordinate vectors explicitly
// template<>
// struct vectorize_struct<CoordinateVector, 4> {
//   using type = std::array<Vec4i, NDIM>;
// };

// template<>
// struct vectorize_struct<CoordinateVector, 8> {
//   using type = std::array<Vec8i, NDIM>;
// };

// template<>
// struct vectorize_struct<CoordinateVector, 16> {
//   using type = std::array<Vec16i, NDIM>;
// };

/// Short version of mapping type to longest possible vector
template <typename T>
using vector_type = typename vectorize_struct<T, hila::vector_info<T>::vector_size>::type;

template <typename T>
void field_storage<T>::allocate_field(const lattice_struct &lattice) {
    if constexpr (hila::is_vectorizable_type<T>::value) {
        fieldbuf = (T *)memalloc(
            lattice.backend_lattice->get_vectorized_lattice<hila::vector_info<T>::vector_size>()
                ->field_alloc_size() *
            sizeof(T));
    } else {
        fieldbuf = (T *)memalloc(sizeof(T) * lattice.field_alloc_size());
    }
}

template <typename T>
void field_storage<T>::free_field() {
#pragma acc exit data delete (fieldbuf)
    if (fieldbuf != nullptr)
        free(fieldbuf);
    fieldbuf = nullptr;
}

// get and set a full vector T

template <typename T>
template <typename vecT>
inline vecT field_storage<T>::get_vector(const unsigned i) const {
    using vectortype = typename hila::vector_info<T>::type;
    using basetype = typename hila::vector_info<T>::base_type;
    constexpr size_t elements = hila::vector_info<T>::elements;
    constexpr size_t vector_size = hila::vector_info<T>::vector_size;
    // using vectorized_type = vector_type<T>;

    static_assert(sizeof(vecT) == sizeof(T) * vector_size);
    // assert (((int64_t)fieldbuf) % ((vector_size)*sizeof(basetype)) == 0);

    vecT value;
    basetype *vp = (basetype *)(fieldbuf) + i * elements * vector_size;
    vectortype *valuep = (vectortype *)(&value);
    for (unsigned e = 0; e < elements; e++) {
        valuep[e].load_a(vp + e * vector_size);
    }
    return value;
}

// note: here i is the vector index

template <typename T>
template <typename vecT>
inline void field_storage<T>::set_vector(const vecT &value, const unsigned i) {
    using vectortype = typename hila::vector_info<T>::type;
    using basetype = typename hila::vector_info<T>::base_type;
    constexpr size_t elements = hila::vector_info<T>::elements;
    constexpr size_t vector_size = hila::vector_info<T>::vector_size;

    static_assert(sizeof(vecT) == sizeof(T) * vector_size);

    basetype *vp = (basetype *)(fieldbuf) + i * elements * vector_size;
    vectortype *valuep = (vectortype *)(&value);
    for (unsigned e = 0; e < elements; e++) {
        valuep[e].store_a(vp + e * vector_size);
    }
}

/// set_element scatters one individual T-element to vectorized store,
/// using the "site" index idx.

template <typename T>
inline void field_storage<T>::set_element(const T &value, const unsigned idx) {
    static_assert(hila::vector_info<T>::is_vectorizable);
    using basetype = typename hila::vector_info<T>::base_type;
    constexpr size_t elements = hila::vector_info<T>::elements;
    constexpr size_t vector_size = hila::vector_info<T>::vector_size;

    // "base" of the vector is (idx/vector_size)*elements; index in vector is idx %
    // vector_size
    basetype *RESTRICT b =
        ((basetype *)(fieldbuf)) + (idx / vector_size) * vector_size * elements + idx % vector_size;
    const basetype *RESTRICT vp = (basetype *)(&value);
    for (unsigned e = 0; e < elements; e++) {
        b[e * vector_size] = vp[e];
    }
}

/// get_element gathers one T-element from vectorized store
/// again, idx is the "site" index

template <typename T>
inline T field_storage<T>::get_element(const unsigned idx) const {
    static_assert(hila::vector_info<T>::is_vectorizable);
    using basetype = typename hila::vector_info<T>::base_type;
    constexpr size_t elements = hila::vector_info<T>::elements;
    constexpr size_t vector_size = hila::vector_info<T>::vector_size;

    static_assert(sizeof(T) == sizeof(basetype) * elements);

    T value;
    // "base" of the vector is (idx/vector_size)*elements; index in vector is idx %
    // vector_size
    const basetype *RESTRICT b =
        (basetype *)(fieldbuf) + (idx / vector_size) * vector_size * elements + idx % vector_size;
    basetype *RESTRICT vp = (basetype *)(&value); // does going through address slow down?
    for (unsigned e = 0; e < elements; e++) {
        vp[e] = b[e * vector_size];
    }
    return value;
}

/// Fetch elements from the field to buffer using sites in index_list
template <typename T>
void field_storage<T>::gather_elements(T *RESTRICT buffer, const unsigned *RESTRICT index_list,
                                       int n, const lattice_struct &lattice) const {

    for (unsigned j = 0; j < n; j++) {
        buffer[j] = get_element(index_list[j]);
    }
}

#ifdef SPECIAL_BOUNDARY_CONDITIONS

template <typename T>
void field_storage<T>::gather_elements_negated(T *RESTRICT buffer,
                                               const unsigned *RESTRICT index_list, int n,
                                               const lattice_struct &lattice) const {
    if constexpr (hila::has_unary_minus<T>::value) {
        for (unsigned j = 0; j < n; j++) {
            buffer[j] = -get_element(index_list[j]); /// requires unary - !!
        }
    } else {
        // sizeof(T) here to prevent compile time evaluation of assert
        assert(sizeof(T) < 1 && "Antiperiodic boundary conditions require that unary - "
                                "-operator is defined!");
    }
}

#endif

/// Vectorized implementation of setting elements
template <typename T>
void field_storage<T>::place_elements(T *RESTRICT buffer, const unsigned *RESTRICT index_list,
                                      int n, const lattice_struct &lattice) {
    for (unsigned j = 0; j < n; j++) {
        set_element(buffer[j], index_list[j]);
    }
}

template <typename T>
void field_storage<T>::set_local_boundary_elements(Direction dir, Parity par,
                                                   const lattice_struct &lattice,
                                                   bool antiperiodic) {

#ifndef SPECIAL_BOUNDARY_CONDITIONS
    assert(!antiperiodic && "antiperiodic only with SPECIAL_BOUNDARY_CONDITIONS");
#endif

    if constexpr (hila::is_vectorizable_type<T>::value) {

        // do the boundary vectorized copy

        constexpr size_t vector_size = hila::vector_info<T>::vector_size;
        constexpr size_t elements = hila::vector_info<T>::elements;
        using vectortype = typename hila::vector_info<T>::type;
        using basetype = typename hila::vector_info<T>::base_type;

        // hila::out0 << "Vectorized boundary dir " << dir << " parity " << (int)par << "
        // bc " << (int)antiperiodic << '\n';

        const auto vector_lattice =
            lattice.backend_lattice
                ->template get_vectorized_lattice<hila::vector_info<T>::vector_size>();
        // The halo copy and permutation is only necessary if vectorization
        // splits the lattice in this Direction or local boundary is copied
        if (vector_lattice->is_boundary_permutation[abs(dir)] ||
            vector_lattice->only_local_boundary_copy[dir]) {

            unsigned start = 0;
            unsigned end = vector_lattice->n_halo_vectors[dir];
            if (par == ODD)
                start = vector_lattice->n_halo_vectors[dir] / 2;
            if (par == EVEN)
                end = vector_lattice->n_halo_vectors[dir] / 2;
            unsigned offset = vector_lattice->halo_offset[dir];

            /// Loop over the boundary sites - i is the vector index
            /// location where the vectors are copied from are in halo_index

            if (vector_lattice->is_boundary_permutation[abs(dir)]) {

                // hila::out0 << "its permutation\n";
                const int *RESTRICT perm = vector_lattice->boundary_permutation[dir];

                basetype *fp = static_cast<basetype *>(static_cast<void *>(fieldbuf));
                for (unsigned idx = start; idx < end; idx++) {
                    /// get ptrs to target and source vec elements
                    basetype *RESTRICT t = fp + (idx + offset) * (elements * vector_size);
                    basetype *RESTRICT s =
                        fp + vector_lattice->halo_index[dir][idx] * (elements * vector_size);

                    if (!antiperiodic) {
                        for (unsigned e = 0; e < elements * vector_size; e += vector_size)
                            for (unsigned i = 0; i < vector_size; i++)
                                t[e + i] = s[e + perm[i]];
                    } else {
#ifdef SPECIAL_BOUNDARY_CONDITIONS
                        for (unsigned e = 0; e < elements * vector_size; e += vector_size)
                            for (unsigned i = 0; i < vector_size; i++)
                                t[e + i] = -s[e + perm[i]];
#endif
                    }
                }
            } else {
                //  hila::out0 << "its not permutation, go for copy: bc " <<
                //  (int)antiperiodic << '\n';
                if (!antiperiodic) {
                    // no boundary permutation, straight copy for all vectors
                    for (unsigned idx = start; idx < end; idx++) {
                        std::memcpy(fieldbuf + (idx + offset) * vector_size,
                                    fieldbuf + vector_lattice->halo_index[dir][idx] * vector_size,
                                    sizeof(T) * vector_size);
                    }
                } else {
#ifdef SPECIAL_BOUNDARY_CONDITIONS
                    basetype *fp = static_cast<basetype *>(static_cast<void *>(fieldbuf));
                    for (unsigned idx = start; idx < end; idx++) {
                        /// get ptrs to target and source vec elements
                        basetype *RESTRICT t = fp + (idx + offset) * (elements * vector_size);
                        basetype *RESTRICT s =
                            fp + vector_lattice->halo_index[dir][idx] * (elements * vector_size);
                        for (unsigned e = 0; e < elements * vector_size; e++)
                            t[e] = -s[e];
                    }
#endif
                }
            }
        }

    } else {
        // now the field is not vectorized.  Std. access copy
        // needed only if b.c. is not periodic

#ifdef SPECIAL_BOUNDARY_CONDITIONS
        if (antiperiodic) {
            // need to copy or do something w. local boundary
            unsigned n, start = 0;
            if (par == ODD) {
                n = lattice.special_boundaries[dir].n_odd;
                start = lattice.special_boundaries[dir].n_even;
            } else {
                if (par == EVEN)
                    n = lattice.special_boundaries[dir].n_even;
                else
                    n = lattice.special_boundaries[dir].n_total;
            }
            unsigned offset = lattice.special_boundaries[dir].offset + start;

            gather_elements_negated(fieldbuf + offset,
                                    lattice.special_boundaries[dir].move_index + start, n, lattice);
        }
#endif
    }
}

// gather full vectors from fieldbuf to buffer, for communications
template <typename T>
void field_storage<T>::gather_comm_vectors(
    T *RESTRICT buffer, const lattice_struct::comm_node_struct &to_node, Parity par,
    const vectorized_lattice_struct<hila::vector_info<T>::vector_size> *RESTRICT vlat,
    bool antiperiodic) const {

    // Use sitelist in to_node, but use only every vector_size -index.  These point to
    // the beginning of the vector
    constexpr size_t vector_size = hila::vector_info<T>::vector_size;
    constexpr size_t elements = hila::vector_info<T>::elements;
    using basetype = typename hila::vector_info<T>::base_type;

    int n;
    const unsigned *index_list = to_node.get_sitelist(par, n);

    assert(n % vector_size == 0);

    if (!antiperiodic) {
        for (unsigned i = 0; i < n; i += vector_size) {
            std::memcpy(buffer + i, fieldbuf + index_list[i], sizeof(T) * vector_size);

            // check that indices are really what they should -- REMOVE
            for (unsigned j = 0; j < vector_size; j++)
                assert(index_list[i] + j == index_list[i + j]);
        }
    } else {
        // copy this as elements
        for (unsigned i = 0; i < n; i += vector_size) {
            basetype *RESTRICT t = static_cast<basetype *>(static_cast<void *>(buffer + i));
            basetype *RESTRICT s =
                static_cast<basetype *>(static_cast<void *>(fieldbuf + index_list[i]));
            for (unsigned e = 0; e < elements * vector_size; e++)
                t[e] = -s[e];
        }
    }
}

// Place the received MPI elements to halo (neighbour) buffer
template <typename T>
void field_storage<T>::place_recv_elements(
    const T *RESTRICT buffer, Direction d, Parity par,
    const vectorized_lattice_struct<hila::vector_info<T>::vector_size> *RESTRICT vlat) const {

    constexpr size_t vector_size = hila::vector_info<T>::vector_size;
    constexpr size_t elements = hila::vector_info<T>::elements;
    using basetype = typename hila::vector_info<T>::base_type;

    unsigned start = 0;
    if (par == ODD)
        start = vlat->recv_list_size[d] / 2;
    unsigned n = vlat->recv_list_size[d];
    if (par != ALL)
        n /= 2;

    // remove const  --  the payload of the buffer remains const, but the halo  bits are
    // changed
    T *targetbuf = const_cast<T *>(fieldbuf);

    for (unsigned i = 0; i < n; i++) {
        unsigned idx = vlat->recv_list[d][i + start];

        basetype *RESTRICT t = ((basetype *)targetbuf) +
                               (idx / vector_size) * vector_size * elements + idx % vector_size;
        const basetype *RESTRICT vp = (basetype *)(&buffer[i]);

        for (unsigned e = 0; e < elements; e++) {
            t[e * vector_size] = vp[e];
        }
    }
}

template <typename T>
void field_storage<T>::free_mpi_buffer(T *buffer) {
    std::free(buffer);
}

template <typename T>
T *field_storage<T>::allocate_mpi_buffer(unsigned n) {
    return (T *)memalloc(n * sizeof(T));
}

#endif