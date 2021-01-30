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
using vector_type = typename vectorize_struct<T,vector_info<T>::vector_size>::type;


template<typename T>
void field_storage<T>::allocate_field( lattice_struct * lattice ) {
  if constexpr (is_vectorizable_type<T>::value) {
    fieldbuf = (T *)memalloc( lattice->backend_lattice->
                              get_vectorized_lattice<vector_info<T>::vector_size>()->field_alloc_size()*sizeof(T) );
  } else {
    fieldbuf = (T*)memalloc( sizeof(T) * lattice->field_alloc_size() );
  }
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
inline vecT field_storage<T>::get_vector(const int i) const {
  using vectortype = typename vector_info<T>::type;
  using basetype = typename vector_info<T>::base_type;
  constexpr int elements = vector_info<T>::elements;
  constexpr int vector_size = vector_info<T>::vector_size;
  // using vectorized_type = vector_type<T>;

  static_assert( sizeof(vecT) == sizeof(T) * vector_size );
  // assert (((int64_t)fieldbuf) % ((vector_size)*sizeof(basetype)) == 0);

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
inline T field_storage<T>::get_element(const int idx) const {
  static_assert( vector_info<T>::is_vectorizable );
  using basetype = typename vector_info<T>::base_type;
  constexpr int elements = vector_info<T>::elements;
  constexpr int vector_size = vector_info<T>::vector_size;

  static_assert( sizeof(T) == sizeof(basetype)*elements );

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

template<typename T>
void field_storage<T>::gather_elements_negated( T * RESTRICT buffer, 
                        const unsigned * RESTRICT index_list, int n,
                        const lattice_struct * RESTRICT lattice) const {
  if constexpr (has_unary_minus<T>::value) {
    for (int j=0; j<n; j++) {
      buffer[j] = - get_element(index_list[j]);    /// requires unary - !!
    }
  } else {
    // sizeof(T) here to prevent compile time evaluation of assert
    assert(sizeof(T) < 1 && 
    "Antiperiodic boundary conditions require that unary - -operator is defined!");
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
void field_storage<T>::set_local_boundary_elements(direction dir, parity par,const lattice_struct * lattice, 
                                                   bool antiperiodic) {

  if constexpr (is_vectorizable_type<T>::value) {

    // do the boundary vectorized copy

    constexpr int vector_size = vector_info<T>::vector_size;
    constexpr int elements = vector_info<T>::elements;
    using vectortype = typename vector_info<T>::type;
    using basetype = typename vector_info<T>::base_type;

    // output0 << "Vectorized boundary dir " << dir << " parity " << (int)par << " bc " << (int)antiperiodic << '\n';

    const auto vector_lattice = 
        lattice->backend_lattice->template get_vectorized_lattice<vector_info<T>::vector_size>();
    // The halo copy and permutation is only necessary if vectorization
    // splits the lattice in this direction or local boundary is copied
    if ( vector_lattice->is_boundary_permutation[abs(dir)] || vector_lattice->only_local_boundary_copy[dir] ) {

      int start = 0;
      int end   = vector_lattice->n_halo_vectors[dir];
      if (par == ODD)  start = vector_lattice->n_halo_vectors[dir]/2;
      if (par == EVEN) end   = vector_lattice->n_halo_vectors[dir]/2;
      int offset = vector_lattice->halo_offset[dir];
    
      /// Loop over the boundary sites - i is the vector index
      /// location where the vectors are copied from are in halo_index

      if (vector_lattice->is_boundary_permutation[abs(dir)]) {

        // output0 << "its permutation\n";
        const int * RESTRICT perm = vector_lattice->boundary_permutation[dir];

        basetype * fp = static_cast<basetype *>(static_cast<void *>(fieldbuf));
        for (int idx=start; idx<end; idx++) {
          /// get ptrs to target and source vec elements
          basetype * RESTRICT t = fp + (idx+offset)*(elements*vector_size);
          basetype * RESTRICT s = fp + vector_lattice->halo_index[dir][idx]*(elements*vector_size);

          if (!antiperiodic) {
            for (int e=0; e<elements*vector_size; e+=vector_size)
              for (int i=0; i<vector_size; i++) 
                t[e + i] = s[e + perm[i]];
          } else {
            for (int e=0; e<elements*vector_size; e+=vector_size)
              for (int i=0; i<vector_size; i++) 
                t[e + i] = -s[e + perm[i]];
          }
        }
      } else {
        //  output0 << "its not permutation, go for copy: bc " << (int)antiperiodic << '\n';
        if (!antiperiodic) {
          // no boundary permutation, straight copy for all vectors
          for (int idx=start; idx<end; idx++) {
            std::memcpy( fieldbuf + (idx+offset)*vector_size, 
                         fieldbuf + vector_lattice->halo_index[dir][idx]*vector_size, 
                         sizeof(T) * vector_size );
          }
        } else {
          basetype * fp = static_cast<basetype *>(static_cast<void *>(fieldbuf));
          for (int idx=start; idx<end; idx++) {
            /// get ptrs to target and source vec elements
            basetype * RESTRICT t = fp + (idx+offset)*(elements*vector_size);
            basetype * RESTRICT s = fp + vector_lattice->halo_index[dir][idx]*(elements*vector_size);
            for (int e=0; e<elements*vector_size; e++)
              t[e] = -s[e];
          }
        }
      }
    }

  } else {
    // now the field is not vectorized.  Std. access copy
    // needed only if b.c. is not periodic 

    if (antiperiodic) {
      // need to copy or do something w. local boundary
      int n, start = 0;
      if (par == ODD) {
        n = lattice->special_boundaries[dir].n_odd;
        start = lattice->special_boundaries[dir].n_even;
      } else {
        if (par == EVEN) n = lattice->special_boundaries[dir].n_even;
        else n = lattice->special_boundaries[dir].n_total;
      }
      int offset = lattice->special_boundaries[dir].offset + start;

      gather_elements_negated( fieldbuf + offset,
                lattice->special_boundaries[dir].move_index + start, n, lattice);
    }
  }
}


// gather full vectors from fieldbuf to buffer, for communications
template <typename T>
void field_storage<T>::gather_comm_vectors( T * RESTRICT buffer, const lattice_struct::comm_node_struct & to_node, 
          parity par, const vectorized_lattice_struct<vector_info<T>::vector_size> * RESTRICT vlat,
          bool antiperiodic) const {
  
  // Use sitelist in to_node, but use only every vector_size -index.  These point to the beginning of the vector
  constexpr int vector_size = vector_info<T>::vector_size;
  constexpr int elements = vector_info<T>::elements;
  using basetype = typename vector_info<T>::base_type;

  int n;
  const unsigned * index_list = to_node.get_sitelist(par,n);
  
  assert(n % vector_size == 0);

  if (!antiperiodic) {
    for (int i=0; i<n; i+=vector_size) {
      std::memcpy( buffer + i, fieldbuf + index_list[i], sizeof(T)*vector_size);

      // check that indices are really what they should -- REMOVE
      for (int j=0; j<vector_size; j++) assert( index_list[i]+j == index_list[i+j] );
    }
  } else {
    // copy this as elements
    for (int i=0; i<n; i+=vector_size) {
      basetype * RESTRICT t = static_cast<basetype *>(static_cast<void *>(buffer + i));
      basetype * RESTRICT s = static_cast<basetype *>(static_cast<void *>(fieldbuf + index_list[i]));
      for (int e=0; e<elements*vector_size; e++)
        t[e] = -s[e];
    }
  }
}



// Place the received MPI elements to halo (neighbour) buffer
template <typename T>
void field_storage<T>::place_recv_elements(const T * RESTRICT buffer, direction d, parity par,
                    const vectorized_lattice_struct<vector_info<T>::vector_size> * RESTRICT vlat) const {

  constexpr int vector_size = vector_info<T>::vector_size;
  constexpr int elements = vector_info<T>::elements;
  using basetype = typename vector_info<T>::base_type;


  int start = 0;
  if (par == ODD) start = vlat->recv_list_size[d]/2;
  int n = vlat->recv_list_size[d];
  if (par != ALL) n /= 2;

  // remove const  --  the payload of the buffer remains const, but the halo  bits are changed
  T * targetbuf = const_cast<T*>(fieldbuf);

  for (int i=0; i<n; i++) {
    int idx = vlat->recv_list[d][i+start];

    basetype * RESTRICT t = ((basetype *)targetbuf) + (idx/vector_size)*vector_size*elements + idx % vector_size;
    const basetype * RESTRICT vp = (basetype *)(&buffer[i]);

    for( int e=0; e<elements; e++ ){
      t[e*vector_size] = vp[e];
    }
  }
}



template <typename T>
void field_storage<T>::free_mpi_buffer( T * buffer){
  std::free(buffer);
}

template <typename T>
T * field_storage<T>::allocate_mpi_buffer( int n ){
  return (T *)memalloc( n * sizeof(T) );
}


#endif