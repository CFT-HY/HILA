#ifndef VANILLA_BACKEND
#define VANILLA_BACKEND

template<typename T>
void field_storage<T>::gather_comm_elements(char * buffer, lattice_struct::comm_node_struct to_node, parity par) const {
    for (int j=0; j<to_node.n_sites(par); j++) {
        T element = get(to_node.site_index(j, par), lattice->field_alloc_size());
        std::memcpy( buffer + j*sizeof(T), (char *) (&element), sizeof(T) );
    }
}
template<typename T>
void field_storage<T>::place_comm_elements(char * buffer, lattice_struct::comm_node_struct from_node, parity par){
    for (int j=0; j<from_node.n_sites(par); j++) {
        T element = *((T*) ( buffer + j*sizeof(T) ));
        set(element, from_node.offset(par)+j);
    }
}
template<typename T>
void field_storage<T>::set_local_boundary_elements(parity par){}

#endif