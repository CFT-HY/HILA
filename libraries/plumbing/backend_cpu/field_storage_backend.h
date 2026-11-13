#ifndef VANILLA_BACKEND
#define VANILLA_BACKEND

#include "../lattice.h"
#include "../field_storage.h"

// Array of Structures layout (standard)

template <typename T>
inline auto field_storage<T>::get(const unsigned i, const unsigned field_alloc_size) const {
    return fieldbuf[i];
}

template <typename T>
// template <typename A>
inline void field_storage<T>::set(const T &value, const unsigned i,
                                  const unsigned field_alloc_size) {
    fieldbuf[i] = value;
}

template <typename T>
void field_storage<T>::allocate_field(const Lattice lattice) {
    fieldbuf = (T *)memalloc(sizeof(T) * lattice->mynode.field_alloc_size);
    if (fieldbuf == nullptr) {
        std::cout << "Failure in Field memory allocation\n";
        exit(1);
    }
#pragma acc enter data create(fieldbuf)
}

template <typename T>
void field_storage<T>::free_field() {
#pragma acc exit data delete (fieldbuf)
    if (fieldbuf != nullptr)
        free(fieldbuf);
    fieldbuf = nullptr;
}

template <typename T>
void field_storage<T>::copy_field(const Lattice lattice, const field_storage<T> &source) {
    memcpy(fieldbuf, source.fieldbuf, sizeof(T) * lattice->mynode.volume);
}


template <typename T>
void field_storage<T>::gather_elements(T *RESTRICT buffer, const unsigned *RESTRICT index_list,
                                       int n, const Lattice lattice) const {

// decorate with pragma omp, should not matter if omp not used
#pragma omp parallel for
    for (unsigned j = 0; j < n; j++) {
        unsigned index = index_list[j];
        buffer[j] = get(index, lattice->mynode.field_alloc_size);
        // std::memcpy( buffer + j, (char *) (&element), sizeof(T) );
    }
}

template <typename T>
void field_storage<T>::place_elements(T *RESTRICT buffer, const unsigned *RESTRICT index_list,
                                      int n, const Lattice lattice) {
#pragma omp parallel for
    for (unsigned j = 0; j < n; j++) {
        set(buffer[j], index_list[j], lattice->mynode.field_alloc_size);
    }
}

template <typename T>
void field_storage<T>::set_local_boundary_elements(Direction dir, Parity par,
                                                   const Lattice lattice) {
    // nothing to do in the standard layout: the neighbours wrap on the node, and the
    // b.c. are applied when the neighbour values are read (hila::bc_nb_value)
}

template <typename T>
void field_storage<T>::gather_comm_elements(T *RESTRICT buffer,
                                            const lattice_struct::comm_node_struct &to_node,
                                            Parity par, const Lattice lattice) const {
    int n;
    const unsigned *index_list = to_node.get_sitelist(par, n);
    gather_elements(buffer, index_list, n, lattice);
}

template <typename T>
inline auto field_storage<T>::get_element(const unsigned i, const Lattice lattice) const {
    return this->get(i, lattice->mynode.field_alloc_size);
}

template <typename T>
template <typename A>
inline void field_storage<T>::set_element(A &value, const unsigned i,
                                          const Lattice lattice) {
    this->set(value, i, lattice->mynode.field_alloc_size);
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
