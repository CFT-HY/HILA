#ifndef VANILLA_BACKEND
#define VANILLA_BACKEND

#include "../lattice.h"
#include "../field_storage.h"

#ifndef layout_SOA
// Array of Structures layout (standard)

template <typename T>
inline auto field_storage<T>::get(const unsigned i,
                                  const unsigned field_alloc_size) const {
    return fieldbuf[i];
}

template <typename T>
// template <typename A>
inline void field_storage<T>::set(const T &value, const unsigned i,
                                  const unsigned field_alloc_size) {
    fieldbuf[i] = value;
}

template <typename T> void field_storage<T>::allocate_field(lattice_struct *lattice) {
    fieldbuf = (T *)memalloc(sizeof(T) * lattice->field_alloc_size());
    if (fieldbuf == nullptr) {
        std::cout << "Failure in Field memory allocation\n";
        exit(1);
    }
#pragma acc enter data create(fieldbuf)
}

template <typename T> void field_storage<T>::free_field() {
#pragma acc exit data delete (fieldbuf)
    if (fieldbuf != nullptr)
        free(fieldbuf);
    fieldbuf = nullptr;
}

#else
// Now structure of arrays

template <typename T>
inline auto field_storage<T>::get(const unsigned i,
                                  const unsigned field_alloc_size) const {
    assert(i < field_alloc_size);
    T value;
    real_t *value_f = (real_t *)(&value);
    for (unsigned e = 0; e < (sizeof(T) / sizeof(real_t)); e++) {
        value_f[e] = ((real_t *)fieldbuf)[e * field_alloc_size + i];
    }
    return value;
}

template <typename T>
// template <typename A>
inline void field_storage<T>::set(const T &value, const unsigned i,
                                  const unsigned field_alloc_size) {
    assert(i < field_alloc_size);
    real_t *value_f = (real_t *)(&value);
    for (unsigned e = 0; e < (sizeof(T) / sizeof(real_t)); e++) {
        ((real_t *)fieldbuf)[e * field_alloc_size + i] = value_f[e];
    }
}

template <typename T> void field_storage<T>::allocate_field(lattice_struct *lattice) {
    constexpr static unsigned t_elements = sizeof(T) / sizeof(real_t);
    fieldbuf = memalloc(t_elements * sizeof(real_t) * lattice->field_alloc_size());
    if (fieldbuf == nullptr) {
        std::cout << "Failure in Field memory allocation\n";
        exit(1);
    }
#pragma acc enter data create(fieldbuf)
}

template <typename T> void field_storage<T>::free_field() {
#pragma acc exit data delete (fieldbuf)
    if (fieldbuf != nullptr)
        free(fieldbuf);
    fieldbuf = nullptr;
}

#endif // layout

template <typename T>
void field_storage<T>::gather_elements(T *RESTRICT buffer,
                                       const unsigned *RESTRICT index_list, int n,
                                       const lattice_struct *RESTRICT lattice) const {
    for (unsigned j = 0; j < n; j++) {
        unsigned index = index_list[j];
        buffer[j] = get(index, lattice->field_alloc_size());
        // std::memcpy( buffer + j, (char *) (&element), sizeof(T) );
    }
}

template <typename T>
void field_storage<T>::gather_elements_negated(
    T *RESTRICT buffer, const unsigned *RESTRICT index_list, int n,
    const lattice_struct *RESTRICT lattice) const {
    if constexpr (has_unary_minus<T>::value) {
        for (unsigned j = 0; j < n; j++) {
            unsigned index = index_list[j];
            buffer[j] = -get(index, lattice->field_alloc_size()); /// requires unary - !!
            // std::memcpy( buffer + j, (char *) (&element), sizeof(T) );
        }
    } else {
        // sizeof(T) here to prevent compile time evaluation of assert
        assert(sizeof(T) < 1 && "Antiperiodic boundary conditions require that unary - "
                                "-operator is defined!");
    }
}

template <typename T>
void field_storage<T>::place_elements(T *RESTRICT buffer,
                                      const unsigned *RESTRICT index_list, int n,
                                      const lattice_struct *RESTRICT lattice) {
    for (unsigned j = 0; j < n; j++) {
        set(buffer[j], index_list[j], lattice->field_alloc_size());
    }
}

template <typename T>
void field_storage<T>::set_local_boundary_elements(Direction dir, Parity par,
                                                   const lattice_struct *RESTRICT lattice,
                                                   bool antiperiodic) {
    // Only need to do something for antiperiodic boundaries

 #ifdef SPECIAL_BOUNDARY_CONDITIONS

    if (antiperiodic) {
        unsigned n, start = 0;
        if (par == ODD) {
            n = lattice->special_boundaries[dir].n_odd;
            start = lattice->special_boundaries[dir].n_even;
        } else {
            if (par == EVEN)
                n = lattice->special_boundaries[dir].n_even;
            else
                n = lattice->special_boundaries[dir].n_total;
        }
        unsigned offset = lattice->special_boundaries[dir].offset + start;
        gather_elements_negated(fieldbuf + offset,
                                lattice->special_boundaries[dir].move_index + start, n,
                                lattice);


    }

#else
    assert(!antiperiodic && "antiperiodic BC possible only if SPECIAL_BOUNDARY_CONDITIONS defined");
#endif

}

template <typename T>
void field_storage<T>::gather_comm_elements(
    T *RESTRICT buffer, const lattice_struct::comm_node_struct &to_node, Parity par,
    const lattice_struct *RESTRICT lattice, bool antiperiodic) const {
    int n;
    const unsigned *index_list = to_node.get_sitelist(par, n);

    if (antiperiodic) {
        gather_elements_negated(buffer, index_list, n, lattice);
    } else {
        gather_elements(buffer, index_list, n, lattice);
    }
}

template <typename T>
auto field_storage<T>::get_element(const unsigned i,
                                   const lattice_struct *RESTRICT lattice) const {
    return this->get(i, lattice->field_alloc_size());
}

template <typename T>
template <typename A>
void field_storage<T>::set_element(A &value, const unsigned i,
                                   const lattice_struct *RESTRICT lattice) {
    this->set(value, i, lattice->field_alloc_size());
}

template <typename T> void field_storage<T>::free_mpi_buffer(T *buffer) {
    std::free(buffer);
}

template <typename T> T *field_storage<T>::allocate_mpi_buffer(unsigned n) {
    return (T *)memalloc(n * sizeof(T));
}

#endif
