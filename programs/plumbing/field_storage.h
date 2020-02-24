#ifndef FIELD_STORAGEH
#define FIELD_STORAGEH

// Pointer to field data and accessors. Only this is passed to the
// CUDA kernels and other accelerators and it only contains a minimal
// amount of data.

#ifndef VECTORIZED // field storage for vectorized types is not defined yet

template <typename T>
class field_storage {
  public:

    #ifdef layout_SOA 
    constexpr static int t_elements = sizeof(T) / sizeof(real_t);
    real_t * fieldbuf = NULL;
    #else
    T * fieldbuf = NULL;
    #endif

    void allocate_field( lattice_struct * lattice );
    void free_field();
    #pragma transformer loop_function
    inline T get(const int i, const int field_alloc_size) const;
    #pragma transformer loop_function
    inline void set(const T &value, const int i, const int field_alloc_size);
    #pragma transformer loop_function
    inline void set(const T &value, unsigned int i);

    void gather_comm_elements(char * buffer, lattice_struct::comm_node_struct to_node, parity par, lattice_struct * lattice) const;
    /// Place boundary elements from neighbour
    void place_comm_elements(char * buffer, lattice_struct::comm_node_struct from_node, parity par, lattice_struct * lattice);
    /// Place boundary elements from local lattice (used in vectorized version)
    void set_local_boundary_elements(parity par, lattice_struct * lattice);

};



#ifdef layout_SOA

template<typename T> 
inline T field_storage<T>::get(const int i, const int field_alloc_size) const {
  assert( idx < field_alloc_size);
  T value;
  real_t *value_f = static_cast<real_t *>(static_cast<void *>(&value));
  for (int i=0; i<(sizeof(T)/sizeof(real_t)); i++) {
     value_f[i] = fieldbuf[i*field_alloc_size + idx];
   }
  return value; 
}

template<typename T>
inline void field_storage<T>::set(const T &value, const int i, const int field_alloc_size){
  assert( idx < field_alloc_size);
  real_t *value_f = static_cast<real_t *>(static_cast<void *>(&value));
  for (int i=0; i<(sizeof(T)/sizeof(real_t)); i++) {
    fieldbuf[i*field_alloc_size + idx] = value_f[i];
  }
}

template<typename T>
inline void field_storage<T>::set(const T &value, unsigned int i) {}

template<typename T>
void field_storage<T>::allocate_field(lattice_struct * lattice) {
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
inline T field_storage<T>::get(const int i, const int field_alloc_size) const 
{
      T value = fieldbuf[i];
      return value;
}

template<typename T>
inline void field_storage<T>::set(const T &value, const int i, const int field_alloc_size)
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
    fieldbuf = (T*) allocate_field_mem( sizeof(T) *lattice->field_alloc_size);
    #pragma acc enter data create(fieldbuf)
}

template<typename T>
void field_storage<T>::free_field() {
  #pragma acc exit data delete(fieldbuf)
  free_field_mem((void *)fieldbuf);
  fieldbuf = nullptr;
}

#endif //layout soa 

#endif //vectorized 


/*
Import backend
*/

#ifdef VECTORIZED

#include "../plumbing/backends/vector_backend.h"

#elif CUDA

#include "../plumbing/backends/cuda_backend.h"

#else

#include "../plumbing/backends/vanilla_backend.h"

#endif

#endif //field storage