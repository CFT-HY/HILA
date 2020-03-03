#ifndef VECTOR_BACKEND
#define VECTOR_BACKEND

#include "../defs.h"
#include "../field_storage.h"

#ifdef AVX512
#define FIELD_ALIGNMENT 64
#elif AVX
#define FIELD_ALIGNMENT 32
#else
#define FIELD_ALIGNMENT 16
#endif


template<typename T>
void field_storage<T>::allocate_field( lattice_struct * lattice ) {
  constexpr int vector_size = field_info<T>::vector_size;
  vectorized_lattice_struct * vlat = lattice->backend_lattice->get_vectorized_lattice(vector_size);

  int size = sizeof(T) * vector_size * vlat->field_alloc_size();
  if (size % FIELD_ALIGNMENT) 
    size = size - (size % FIELD_ALIGNMENT) + FIELD_ALIGNMENT;
  fieldbuf = (T*) aligned_alloc( FIELD_ALIGNMENT, size);
}

template<typename T>
void field_storage<T>::free_field() {
  #pragma acc exit data delete(fieldbuf)
  if(fieldbuf != nullptr)
    free(fieldbuf);
  fieldbuf = nullptr;
}


template<typename T>
auto field_storage<T>::get(const int i, const int field_alloc_size) const
{
  using vectortype = typename field_info<T>::base_vector_type;
  using basetype = typename field_info<T>::base_type;
  constexpr int elements = field_info<T>::elements;
  constexpr int vector_size = field_info<T>::vector_size;

  typename field_info<T>::vector_type value;
  basetype *vp = (basetype *) (fieldbuf) + i*elements*vector_size;
  vectortype *valuep = (vectortype *)(&value);
  for( int e=0; e<elements; e++ ){
    valuep[e].load(vp+e*vector_size);
  }
  return value;
}



template<typename T>
template<typename A>
inline void field_storage<T>::set(const A &value, const int i, const int field_alloc_size) 
{
  using vectortype = typename field_info<T>::base_vector_type;
  using basetype = typename field_info<T>::base_type;
  constexpr int elements = field_info<T>::elements;
  constexpr int vector_size = field_info<T>::vector_size;

  basetype *vp = (basetype *) (fieldbuf) + i*elements*vector_size;
  vectortype *valuep = (vectortype *)(&value);
  for( int e=0; e<elements; e++ ){
    valuep[e].store((vp + e*vector_size));
  }
}









/// Vectorized implementation of fetching boundary elements
/* Gathers sites at the boundary that need to be communicated to neighbours */
template<typename T>
void field_storage<T>::gather_comm_elements(char * buffer, lattice_struct::comm_node_struct to_node, parity par, lattice_struct * lattice) const {
  vectorized_lattice_struct * vlat = lattice->backend_lattice->get_vectorized_lattice(field_info<T>::vector_size);
  for (int j=0; j<to_node.n_sites(par); j++) {
    int index = to_node.site_index(j, par);
    int v_index = vlat->vector_index[index];
    auto element = get(vlat->lattice_index[index], vlat->field_alloc_size());
    auto pvector = (typename field_info<T>::base_vector_type*) (&element);

    for( int e=0; e<field_info<T>::elements; e++ ){
      auto basenumber = pvector[e].extract(v_index);
      auto * number_buffer = (typename field_info<T>::base_type *) buffer;

      number_buffer[j*field_info<T>::elements + e] = basenumber;
    }
  }
}

/// Vectorized implementation of setting boundary elements
/* Sets the values the neighbour elements from the communication buffer */
template<typename T>
void field_storage<T>::place_comm_elements(char * buffer, lattice_struct::comm_node_struct from_node, parity par, lattice_struct * lattice){
  vectorized_lattice_struct * vlat = lattice->backend_lattice->get_vectorized_lattice(field_info<T>::vector_size);
  for (int j=0; j<from_node.n_sites(par); j++) {
    int index = from_node.offset(par)+j;
    int v_index = vlat->vector_index[index];
    auto element = get(vlat->lattice_index[index], vlat->field_alloc_size());
    auto pvector = (typename field_info<T>::base_vector_type*) (&element);

    for( int e=0; e<field_info<T>::elements; e++ ){
      auto number_buffer = (typename field_info<T>::base_type *) buffer;
      auto basenumber = number_buffer[j*field_info<T>::elements + e];

      pvector[e].insert(v_index, basenumber);
    }
    set(element, vlat->lattice_index[index], vlat->field_alloc_size());
  }
}

template<typename T>
void field_storage<T>::set_local_boundary_elements(parity par, lattice_struct * lattice){
  constexpr int vector_size = field_info<T>::vector_size;
  constexpr int elements = field_info<T>::elements;
  using vectortype = typename field_info<T>::base_vector_type;
  using basetype = typename field_info<T>::base_type;
  vectorized_lattice_struct * vlat = lattice->backend_lattice->get_vectorized_lattice(vector_size);
  // Loop over the boundary sites
  for( vectorized_lattice_struct::halo_site hs: vlat->halo_sites )
  if(par == ALL || par == hs.par ) {
    int *perm = vlat->boundary_permutation[hs.dir];
    auto temp = get(hs.nb_index, vlat->field_alloc_size());
    auto dest = (basetype *) (fieldbuf) + elements*vector_size*(vlat->sites + hs.halo_index);
    vectortype * e = (vectortype*) &temp;
    basetype * d = (basetype*) dest;
    for( int v=0; v<elements; v++ ){
      basetype t1[vector_size], t2[vector_size];
      e[v].store(&(t1[0]));
      for( int i=0; i<vector_size; i++ )
       d[v*vector_size+i] =  t1[perm[i]];
    }
  }
}

#endif