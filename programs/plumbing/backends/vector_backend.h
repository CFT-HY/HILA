#ifndef VECTOR_BACKEND
#define VECTOR_BACKEND

/// Vectorized implementation of fetching boundary elements
/* Gathers sites at the boundary that need to be communicated to neighbours */
template<typename T>
void field_storage<T>::gather_comm_elements(char * buffer, lattice_struct::comm_node_struct to_node, parity par) const {
  for (int j=0; j<to_node.n_sites(par); j++) {
    int index = to_node.site_index(j, par);
    int v_index = lattice->vector_index[index];
    auto element = get(lattice->lattice_index[index], lattice->field_alloc_size());
    auto pvector = (typename field_info<T>::vector_type*) (&element);

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
void field_storage<T>::place_comm_elements(char * buffer, lattice_struct::comm_node_struct from_node, parity par){
  for (int j=0; j<from_node.n_sites(par); j++) {
    int index = from_node.offset(par)+j;
    int v_index = lattice->vector_index[index];
    auto element = get(lattice->lattice_index[index], lattice->field_alloc_size());
    auto pvector = (typename field_info<T>::vector_type*) (&element);

    for( int e=0; e<field_info<T>::elements; e++ ){
      auto number_buffer = (typename field_info<T>::base_type *) buffer;
      auto basenumber = number_buffer[j*field_info<T>::elements + e];

      pvector[e].insert(v_index, basenumber);
    }
    set(element, lattice->lattice_index[index]);
  }
}

template<typename T>
void field_storage<T>::set_local_boundary_elements(parity par){
  constexpr int vector_size = field_info<T>::vector_size;
  constexpr int elements = field_info<T>::elements;
  using vectortype = typename field_info<T>::vector_type;
  using basetype = typename field_info<T>::base_type;
  // Loop over the boundary sites
  for( vectorized_lattice_struct::halo_site hs: lattice->halo_sites )
  if(par == ALL || par == hs.par ) {
    int *perm = lattice->boundary_permutation[hs.dir];
    //vectorized::permute<vector_size>(lattice->boundary_permutation[hs.dir], &temp, elements);
    auto temp = get(hs.nb_index);
    auto dest = payload.fieldbuf + lattice->sites + hs.halo_index;
    vectortype * e = (vectortype*) &temp;
    basetype * d = (basetype*) dest;
    for( int v=0; v<elements; v++ ){
      basetype t1[vector_size], t2[vector_size];
      e[v].store(&(t1[0]));
      for( int i=0; i<vector_size; i++ )
         d[v*vector_size+i] =  t1[perm[i]];
    }
    //set(temp, lattice->sites + hs.halo_index);
  }
}

#endif