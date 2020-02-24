#ifndef VANILLA_BACKEND
#define VANILLA_BACKEND


#ifdef layout_SOA

template<typename T> 
inline auto field_storage<T>::get(const int i, const int field_alloc_size) const {
  assert( idx < field_alloc_size);
  T value;
  real_t *value_f = static_cast<real_t *>(static_cast<void *>(&value));
  for (int i=0; i<(sizeof(T)/sizeof(real_t)); i++) {
     value_f[i] = fieldbuf[i*field_alloc_size + idx];
   }
  return value; 
}

template<typename T>
template<typename A>
inline void field_storage<T>::set(const A &value, const int i, const int field_alloc_size){
  assert( idx < field_alloc_size);
  real_t *value_f = static_cast<real_t *>(static_cast<void *>(&value));
  for (int i=0; i<(sizeof(T)/sizeof(real_t)); i++) {
    fieldbuf[i*field_alloc_size + idx] = value_f[i];
  }
}

template<typename T>
void field_storage<T>::allocate_field(lattice_struct * lattice) {
  constexpr static int t_elements = sizeof(T) / sizeof(real_t);
  fieldbuf = (real_t*) allocate_field_mem( t_elements*sizeof(real_t) * lattice->field_alloc_size );
  #pragma acc enter data create(fieldbuf)
}

template<typename T>
void field_storage<T>::free_field() {
  #pragma acc exit data delete(fieldbuf)
  free_field_mem((void *)fieldbuf);
  fieldbuf = nullptr;
}

#else 

template<typename T> 
inline auto field_storage<T>::get(const int i, const int field_alloc_size) const 
{
      T value = fieldbuf[i];
      return value;
}

template<typename T>
template<typename A>
inline void field_storage<T>::set(const A &value, const int i, const int field_alloc_size)
{
  fieldbuf[i] = value;
}

template<typename T>
inline void field_storage<T>::set(const T &value, unsigned int i)
{
  fieldbuf[i] = value;
}

template<typename T>
void field_storage<T>::allocate_field(lattice_struct * lattice) {
    fieldbuf = (T*) allocate_field_mem( sizeof(T) *lattice->field_alloc_size() );
    #pragma acc enter data create(fieldbuf)
}

template<typename T>
void field_storage<T>::free_field() {
  #pragma acc exit data delete(fieldbuf)
  free_field_mem((void *)fieldbuf);
  fieldbuf = nullptr;
}

#endif //layout soa 



template<typename T>
void field_storage<T>::gather_comm_elements(char * buffer, lattice_struct::comm_node_struct to_node, parity par, lattice_struct * lattice) const {
    for (int j=0; j<to_node.n_sites(par); j++) {
        T element = get(to_node.site_index(j, par), lattice->field_alloc_size());
        std::memcpy( buffer + j*sizeof(T), (char *) (&element), sizeof(T) );
    }
}
template<typename T>
void field_storage<T>::place_comm_elements(char * buffer, lattice_struct::comm_node_struct from_node, parity par, lattice_struct * lattice){
    for (int j=0; j<from_node.n_sites(par); j++) {
        T element = *((T*) ( buffer + j*sizeof(T) ));
        set(element, from_node.offset(par)+j);
    }
}
template<typename T>
void field_storage<T>::set_local_boundary_elements(parity par, lattice_struct * lattice){}


#endif