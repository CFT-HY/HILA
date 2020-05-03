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
template<typename A, int vector_size, class Enable = void>
struct vectorize_struct{};

/// A is a basic type, so just return the matching vector type
template<typename A, int vector_size>
struct vectorize_struct<A, vector_size, typename std::enable_if_t<is_arithmetic<A>::value>> {
  using type = typename vector_base_type<A,vector_size>::type;
};

// B is a templated class, so construct a vectorized type
template<template<typename B> class C, typename B, int vector_size>
struct vectorize_struct<C<B>, vector_size>{
  using vectorized_B = typename vectorize_struct<B, vector_size>::type;
  using type = C<vectorized_B>;
};

template<template<int a, typename B> class C, int a, typename B, int vector_size>
struct vectorize_struct<C<a,B>, vector_size>{
  using vectorized_B = typename  vectorize_struct<B, vector_size>::type;
  using type = C<a, vectorized_B>;
};

template<template<int a, int b, typename B> class C, int a, int b,  typename B, int vector_size>
struct vectorize_struct<C<a,b,B>, vector_size>{
  using vectorized_B = typename  vectorize_struct<B, vector_size>::type;
  using type = C<a, b, vectorized_B>;
};


/// Match coordinate vectors explicitly
template<>
struct vectorize_struct<coordinate_vector, 4> {
  using type = std::array<Vec4i, NDIM>;
};

template<>
struct vectorize_struct<coordinate_vector, 8> {
  using type = std::array<Vec8i, NDIM>;
};

template<>
struct vectorize_struct<coordinate_vector, 16> {
  using type = std::array<Vec16i, NDIM>;
};



/// Short version of mapping type to longest possible vector
template <typename T>
using vector_type = typename vectorize_struct<T,vector_info<T>::vector_size>::type;


template<typename T>
void field_storage<T>::allocate_field( lattice_struct * lattice ) {
  fieldbuf = (T *)memalloc( lattice->backend_lattice->
                            get_vectorized_lattice<vector_info<T>::vector_size>()->field_alloc_size()*sizeof(T) );
}

template<typename T>
void field_storage<T>::free_field() {
  #pragma acc exit data delete(fieldbuf)
  if(fieldbuf != nullptr)
    free(fieldbuf);
  fieldbuf = nullptr;
}

// get and set a full vector T

template<typename T>
template<typename vecT>
#pragma transformer loop_function
inline vecT field_storage<T>::get_vector(const int i) const {
  using vectortype = typename vector_info<T>::type;
  using basetype = typename vector_info<T>::base_type;
  constexpr int elements = vector_info<T>::elements;
  constexpr int vector_size = vector_info<T>::vector_size;
  // using vectorized_type = vector_type<T>;

  static_assert( sizeof(vecT) == sizeof(T) * vector_size );

  vecT value;
  basetype *vp = (basetype *) (fieldbuf) + i*elements*vector_size;
  vectortype *valuep = (vectortype *)(&value);
  for( int e=0; e<elements; e++ ){
    valuep[e].load_a(vp+e*vector_size);
  }
  return value;
}


// note: here i is the vector index

template<typename T>
template<typename vecT>
#pragma transformer loop_function
inline void field_storage<T>::set_vector(const vecT &value, const int i)  {
  using vectortype = typename vector_info<T>::type;
  using basetype = typename vector_info<T>::base_type;
  constexpr int elements = vector_info<T>::elements;
  constexpr int vector_size = vector_info<T>::vector_size;

  static_assert( sizeof(vecT) == sizeof(T) * vector_size );

  basetype *vp = (basetype *) (fieldbuf) + i*elements*vector_size;
  vectortype *valuep = (vectortype *)(&value);
  for( int e=0; e<elements; e++ ){
    valuep[e].store_a(vp + e*vector_size);
  }
}


/// set_element scatters one individual T-element to vectorized store,
/// using the "site" index idx.

template<typename T>
#pragma transformer loop_function
inline void field_storage<T>::set_element(const T &value, const int idx) {
  static_assert( vector_info<T>::is_vectorizable );
  using basetype = typename vector_info<T>::base_type;
  constexpr int elements = vector_info<T>::elements;
  constexpr int vector_size = vector_info<T>::vector_size;

  // "base" of the vector is (idx/vector_size)*elements; index in vector is idx % vector_size
  basetype * RESTRICT b = ((basetype *) (fieldbuf)) + (idx/vector_size)*vector_size*elements + idx % vector_size;
  const basetype * RESTRICT vp = (basetype *)(&value);
  for( int e=0; e<elements; e++ ){
    b[e*vector_size] = vp[e];
  }
}

/// get_element fetches one T-element from vectorized store
/// again, idx is the "site" index

template<typename T>
#pragma transformer loop_function
inline T field_storage<T>::get_element(const int idx) const {
  static_assert( vector_info<T>::is_vectorizable );
  using basetype = typename vector_info<T>::base_type;
  constexpr int elements = vector_info<T>::elements;
  constexpr int vector_size = vector_info<T>::vector_size;

  T value;
  // "base" of the vector is (idx/vector_size)*elements; index in vector is idx % vector_size
  const basetype * RESTRICT b = (basetype *) (fieldbuf) + (idx/vector_size)*vector_size*elements + idx % vector_size;
  basetype * RESTRICT vp = (basetype *)(&value);   // does going through address slow down?
  for( int e=0; e<elements; e++ ){
    vp[e] = b[e*vector_size];
  }
  return value;
}


/// Fetch elements from the field to buffer using sites in index_list
template<typename T>
void field_storage<T>::gather_elements(T * RESTRICT buffer, const unsigned * RESTRICT index_list, 
                                       int n, const lattice_struct * RESTRICT lattice) const {

  for (int j=0; j<n; j++) {
    buffer[j] = get_element(index_list[j]);
  }
}



/// Vectorized implementation of setting elements
template<typename T>
void field_storage<T>::place_elements(T * RESTRICT buffer, const unsigned * RESTRICT index_list, int n,
                                      const lattice_struct * RESTRICT lattice) {
  for (int j=0; j<n; j++) {
    set_element(buffer[j],index_list[j]);
  }
}


template<typename T>
void field_storage<T>::set_local_boundary_elements(direction dir, parity par, lattice_struct * lattice){
  constexpr int vector_size = vector_info<T>::vector_size;
  constexpr int elements = vector_info<T>::elements;
  using vectortype = typename vector_info<T>::type;
  using basetype = typename vector_info<T>::base_type;

  const auto vector_lattice = 
      lattice->backend_lattice->get_vectorized_lattice<vector_info<T>::vector_size>();
  // The halo copy and permutation is only necessary if vectorization
  // splits the lattice in this direction
  if ( vector_lattice->is_boundary_permutation[abs(dir)]) {

    hila::output << "BOundary elems to dir " << dir << " parity " << (int)par << '\n';

    int start = 0;
    int end   = vector_lattice->n_halo_vectors[dir];
    if (par == ODD)  start = vector_lattice->n_halo_vectors[dir]/2;
    if (par == EVEN) end   = vector_lattice->n_halo_vectors[dir]/2;
    int offset = vector_lattice->halo_offset[dir];

    /// Loop over the boundary sites - i is the vector index
    /// location where the vectors are copied from are in halo_index

    const int * RESTRICT perm = vector_lattice->boundary_permutation[dir];

    basetype * fp = static_cast<basetype *>(static_cast<void *>(fieldbuf));
    for (int idx=start; idx<end; idx++) {
      /// get ptrs to target and source vec elements
      basetype * RESTRICT t = fp + (idx+offset)*(elements*vector_size);
      basetype * RESTRICT s = fp + vector_lattice->halo_index[dir][idx]*(elements*vector_size);

      for (int e=0; e<elements*vector_size; e+=vector_size)
        for (int i=0; i<vector_size; i++) 
          t[e + i] = s[e + perm[i]];

    }
 
  }
}

#endif