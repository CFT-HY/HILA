#ifndef VECTOR_BACKEND
#define VECTOR_BACKEND

#include "../defs.h"
#include "../lattice.h"
#include "../field_storage.h"
#include "vector_types.h"




/// Replaces basetypes with vectors in a given templated class

/// First base definition for replace_type, which recursively looks for the
/// base type and replaces it in the end
/// General template, never matched
template<typename A, int vector_size, class Enable = void>
struct vectorize_struct{};

/// A is a basic type, so just return the matching vector type
template<typename A, int vector_size>
struct vectorize_struct<A, vector_size, typename std::enable_if_t<is_arithmetic<A>::value>> {
  using type = typename vector_base_type<A,vector_size>::type;
};

/// B is a templated class, so construct a vectorized type
template<template<typename B> class C, typename B, int vector_size>
struct vectorize_struct<C<B>, vector_size>{
  using vectorized_B = typename vectorize_struct<B, vector_size>::type;
  using type = C<vectorized_B>;
};

/// C is a templated class with an integer parameter
template<template<int a, typename B> class C, int a, typename B, int vector_size>
struct vectorize_struct<C<a,B>, vector_size>{
  using vectorized_B = typename  vectorize_struct<B, vector_size>::type;
  using type = C<a, vectorized_B>;
};

/// C is a templated class with two integer parameters
template<template<int a, int b, typename B> class C, int a, int b,  typename B, int vector_size>
struct vectorize_struct<C<a,b,B>, vector_size>{
  using vectorized_B = typename  vectorize_struct<B, vector_size>::type;
  using type = C<a, b, vectorized_B>;
};


/// Match coordinate vectors explicitly
template<>
struct vectorize_struct<CoordinateVector, 4> {
  using type = std::array<Vec4i, NDIM>;
};

/// Match coordinate vectors explicitly
template<>
struct vectorize_struct<CoordinateVector, 8> {
  using type = std::array<Vec8i, NDIM>;
};

/// Match coordinate vectors explicitly
template<>
struct vectorize_struct<CoordinateVector, 16> {
  using type = std::array<Vec16i, NDIM>;
};



/// Short version of mapping type to longest possible vector
template <typename T>
using vector_type = typename vectorize_struct<T,vector_info<T>::vector_size>::type;






template<typename T>
void field_storage<T>::allocate_field( lattice_struct * lattice ) {
  constexpr int vector_size = vector_info<T>::vector_size;
  vectorized_lattice_struct<vector_size> * vlat = lattice->backend_lattice->get_vectorized_lattice<vector_size>();

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
  constexpr int elements = vector_info<T>::elements;
  constexpr int vector_size = vector_info<T>::vector_size;
  using vectorized_type = vector_type<T>;

  vectorized_type value;
  basetype *vp = (basetype *) (fieldbuf) + i*elements*vector_size;
  vectortype *valuep = (vectortype *)(&value);
  for( int e=0; e<elements; e++ ){
    valuep[e].load_a(vp+e*vector_size);
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
    valuep[e].store_a(vp + e*vector_size);
  }
}







/// Vectorized implementation of fetching elements
template<typename T>
void field_storage<T>::gather_elements(char * RESTRICT buffer, const unsigned * RESTRICT index_list, 
                                       int n, const lattice_struct * RESTRICT lattice) const {
  constexpr int vector_size = vector_info<T>::vector_size;
  vectorized_lattice_struct<vector_size> * vlat = lattice->backend_lattice->get_vectorized_lattice<vector_info<T>::vector_size>();
  for (int j=0; j<n; j++) {
    int index = index_list[j];
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



/// Vectorized implementation of setting elements
template<typename T>
void field_storage<T>::place_elements(char * RESTRICT buffer, const unsigned * RESTRICT index_list, int n,
                                      const lattice_struct * RESTRICT lattice) {
  constexpr int vector_size = vector_info<T>::vector_size;
  const vectorized_lattice_struct<vector_size> * RESTRICT vlat = 
      lattice->backend_lattice->get_vectorized_lattice<vector_info<T>::vector_size>();
  for (int j=0; j<n; j++) {
    int index = index_list[j];
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
void field_storage<T>::set_local_boundary_elements(direction dir, parity par, lattice_struct * lattice){
  constexpr int vector_size = vector_info<T>::vector_size;
  constexpr int elements = vector_info<T>::elements;
  using vectortype = typename vector_info<T>::type;
  using basetype = typename vector_info<T>::base_type;
  const vectorized_lattice_struct<vector_size> * RESTRICT vlat = 
    lattice->backend_lattice->get_vectorized_lattice<vector_size>();

  // The halo copy and permutation is only necessary if vectorization
  // splits the lattice in this direction
  if( vlat->split[dir] > 1 ) {

    // Loop over the boundary sites
    for( parity p : loop_parities(par)){
      int par_int = (int) p - 1; 
      const int * RESTRICT perm = vlat->boundary_permutation[dir];
      auto hs = vlat->halo_sites[par_int][dir];

      for( int idx = 0; idx < hs.nb_index.size(); idx++ ) {
        basetype * fp = static_cast<basetype *>(fieldbuf);
        basetype * RESTRICT s = fp + elements*vector_size*hs.nb_index[idx];
        basetype * RESTRICT d = fp + elements*vector_size*(vlat->sites + hs.first_index+idx);
        for( int e=0; e<elements; e++ ){
          for( int i=0; i<vector_size; i++ )
           d[e*vector_size+i] =  s[e*vector_size + perm[i]];
        }
      }
    }
  }
}

#endif