#ifndef VECTOR_BACKEND
#define VECTOR_BACKEND

#include "../defs.h"
#include "../field_storage.h"

#ifndef VECTOR_SIZE
#define VECTOR_SIZE 32
#endif


/// Utility for selecting a vector type
template<typename T>
struct vector_base_type_struct {};

/// Specializations of the vector type selector
#if VECTOR_SIZE == 32
template<>
struct vector_base_type_struct<double> {
  using type = Vec4d;
  static constexpr int vector_size = 4;
};

template<>
struct vector_base_type_struct<float> {
  using type = Vec8f;
  static constexpr int vector_size = 8;
};

template<>
struct vector_base_type_struct<int> {
  using type = Vec8i;
  static constexpr int vector_size = 8;
};

template<>
struct vector_base_type_struct<coordinate_vector> {
  using type = Vec8i;
  static constexpr int vector_size = 8;
};


#elif VECTOR_SIZE == 64
template<typename T>
struct vector_base_type_struct {
};

template<>
struct vector_base_type_struct<double> {
  using type = Vec8d;
  static constexpr int vector_size = 8;
};

template<>
struct vector_base_type_struct<float> {
  using type = Vec16f;
  static constexpr int vector_size = 16;
};

template<>
struct vector_base_type_struct<int> {
  using type = Vec16i;
  static constexpr int vector_size = 16;
};

template<>
struct vector_base_type_struct<coordinate_vector> {
  using type = Vec16i;
  static constexpr int vector_size = 16;
};

#endif


/// Get the vector type for any larger type (such as cmplx, matrix...)
template<typename T>
struct vector_info{
  using base_type = typename basetypestruct<T>::type; 
  using type = typename vector_base_type_struct<base_type>::type;

  static constexpr int vector_size = vector_base_type_struct<base_type>::vector_size;
  static constexpr int elements = sizeof(T)/sizeof(base_type);
  static constexpr int base_type_size = sizeof(base_type);
};



/// Replaces basetypes with vectors in a given templated class

/// First base definition for replace_type, which recursively looks for the
/// base type and replaces it in the end
/// Here is not a basic type, so we need to replace it
template<typename A, class Enable = void>
struct vectorize_struct{};

/// A is a basic type, so just return the matching vector type
template<typename A>
struct vectorize_struct<A, typename std::enable_if_t<is_arithmetic<A>::value>> {
  using type = typename vector_info<A>::type;
};


// B is a templated class, so construct a vectorized type
template<template<typename B> class C, typename B>
struct vectorize_struct<C<B>>{
  using vectorized_B = typename vectorize_struct<B>::type;
  using type = C<vectorized_B>;
};

template<template<int a, typename B> class C, int a, typename B>
struct vectorize_struct<C<a,B>>{
  using vectorized_B = typename  vectorize_struct<B>::type;
  using type = C<a, vectorized_B>;
};

template<template<int a, int b, typename B> class C, int a, int b,  typename B>
struct vectorize_struct<C<a,b,B>>{
  using vectorized_B = typename  vectorize_struct<B>::type;
  using type = C<a, b, vectorized_B>;
};





template<typename T>
void field_storage<T>::allocate_field( lattice_struct * lattice ) {
  constexpr int vector_size = vector_info<T>::vector_size;
  vectorized_lattice_struct * vlat = lattice->backend_lattice->get_vectorized_lattice(vector_size);

  int size = sizeof(T) * vector_size * vlat->field_alloc_size();
  if (size % VECTOR_SIZE) 
    size = size - (size % VECTOR_SIZE) + VECTOR_SIZE;
  fieldbuf = (T*) aligned_alloc( VECTOR_SIZE, size);
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
  using vectortype = typename vector_info<T>::type;
  using basetype = typename vector_info<T>::base_type;
  using vectorized_type = typename vectorize_struct<T>::type;
  constexpr int elements = vector_info<T>::elements;
  constexpr int vector_size = vector_info<T>::vector_size;

  vectorized_type value;
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
  using vectortype = typename vector_info<T>::type;
  using basetype = typename vector_info<T>::base_type;
  constexpr int elements = vector_info<T>::elements;
  constexpr int vector_size = vector_info<T>::vector_size;

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
  vectorized_lattice_struct * vlat = lattice->backend_lattice->get_vectorized_lattice(vector_info<T>::vector_size);
  for (int j=0; j<to_node.n_sites(par); j++) {
    int index = to_node.site_index(j, par);
    int v_index = vlat->vector_index[index];
    auto element = get(vlat->lattice_index[index], vlat->field_alloc_size());
    auto pvector = (typename vector_info<T>::type*) (&element);

    for( int e=0; e<vector_info<T>::elements; e++ ){
      auto basenumber = pvector[e].extract(v_index);
      auto * number_buffer = (typename vector_info<T>::base_type *) buffer;

      number_buffer[j*vector_info<T>::elements + e] = basenumber;
    }
  }
}

/// Vectorized implementation of setting boundary elements
/* Sets the values the neighbour elements from the communication buffer */
template<typename T>
void field_storage<T>::place_comm_elements(char * buffer, lattice_struct::comm_node_struct from_node, parity par, lattice_struct * lattice){
  vectorized_lattice_struct * vlat = lattice->backend_lattice->get_vectorized_lattice(vector_info<T>::vector_size);
  for (int j=0; j<from_node.n_sites(par); j++) {
    int index = from_node.offset(par)+j;
    int v_index = vlat->vector_index[index];
    auto element = get(vlat->lattice_index[index], vlat->field_alloc_size());
    auto pvector = (typename vector_info<T>::type*) (&element);

    for( int e=0; e<vector_info<T>::elements; e++ ){
      auto number_buffer = (typename vector_info<T>::base_type *) buffer;
      auto basenumber = number_buffer[j*vector_info<T>::elements + e];

      pvector[e].insert(v_index, basenumber);
    }
    set(element, vlat->lattice_index[index], vlat->field_alloc_size());
  }
}

template<typename T>
void field_storage<T>::set_local_boundary_elements(parity par, lattice_struct * lattice){
  constexpr int vector_size = vector_info<T>::vector_size;
  constexpr int elements = vector_info<T>::elements;
  using vectortype = typename vector_info<T>::type;
  using basetype = typename vector_info<T>::base_type;
  vectorized_lattice_struct * vlat = lattice->backend_lattice->get_vectorized_lattice(vector_size);
  // Loop over the boundary sites
  for( vectorized_lattice_struct::halo_site hs: vlat->halo_sites )
  if(par == ALL || par == hs.par ) {
    int *perm = vlat->boundary_permutation[hs.dir];
    basetype * s = (basetype *) (fieldbuf) + elements*vector_size*hs.nb_index;
    basetype * d = (basetype *) (fieldbuf) + elements*vector_size*(vlat->sites + hs.halo_index);
    for( int e=0; e<elements; e++ ){
      for( int i=0; i<vector_size; i++ )
       d[e*vector_size+i] =  s[e*vector_size + perm[i]];
    }
  }
}

#endif