#ifndef CUDA_BACKEND
#define CUDA_BACKEND

#include "../defs.h"
#include "../field_storage.h"
#include <tgmath.h>

/* CUDA implementations */
template<typename T>
void field_storage<T>::allocate_field(lattice_struct * lattice) {
  auto status = cudaMalloc( (void **)&fieldbuf, sizeof(T) * lattice->field_alloc_size() );
  check_cuda_error(status, "Allocate field memory");
  if (fieldbuf == nullptr) {
    std::cout << "Failure in field memory allocation\n";
  }
  assert(fieldbuf != nullptr);
}

template<typename T>
void field_storage<T>::free_field() {
  if (fieldbuf != nullptr){
    auto status = cudaFree(fieldbuf);
    check_cuda_error(status, "Free field memory");
  }
  fieldbuf = nullptr;
}



// Only attempt to compile with CUDA compiler.
// Hilapp will skip these.
#ifdef __CUDACC__


// These are used in device code. Can be called directly in a kernel.
template<typename T> 
__device__ auto field_storage<T>::get(const int i, const int field_alloc_size) const {
  assert( i < field_alloc_size);
  using base_type = typename base_type_struct<T>::type;
  constexpr int n_elements = sizeof(T) / sizeof(base_type);
  T value;
  base_type * value_f = (base_type *)&value;
  base_type * fp = (base_type *)(fieldbuf);
  for (int e=0; e<n_elements; e++) {
    value_f[e] = fp[e*field_alloc_size + i];
  }
  return value;
}


template<typename T>
template<typename A>
__device__ inline void field_storage<T>::set(const A &value, const int i, const int field_alloc_size){
  assert( i < field_alloc_size);
  using base_type = typename base_type_struct<T>::type;
  constexpr int n_elements =sizeof(T) / sizeof(base_type);
  const base_type * value_f = (base_type *)&value;
  base_type * fp = (base_type *)(fieldbuf);
  for (int e=0; e<n_elements; e++) {
    fp[e*field_alloc_size + i] = value_f[e];
  }
}



/// Get a single element from the field outside a loop. Slow, should only be used for setup
template <typename T>
__global__ void get_element_kernel( field_storage<T> field, char *buffer, unsigned i, const int field_alloc_size )
{
  *((T*) buffer) = field.get(i, field_alloc_size);
}


template<typename T>
auto field_storage<T>::get_element( const int i, const lattice_struct * RESTRICT lattice) const {
  char * d_buffer;
  T value;
  
  // Call the kernel to collect the element
  cudaMalloc( (void **)&(d_buffer), sizeof(T));
  get_element_kernel<<< 1, 1 >>>(*this, d_buffer, i, lattice->field_alloc_size() );
  
  // Copy the result to the host
  cudaMemcpy((char *) (&value), d_buffer, sizeof(T), cudaMemcpyDeviceToHost );
  cudaFree(d_buffer);
  return value;
}


/// Set a single element from outside a loop. Slow, should only be used for setup
template <typename T>
__global__ void set_element_kernel( field_storage<T> field, char *buffer, unsigned i, const int field_alloc_size )
{
  field.set( (T*) buffer, i, field_alloc_size);
}

template<typename T>
template<typename A>
void field_storage<T>::set_element(A &value, const int i, const lattice_struct * RESTRICT lattice) {
  char * d_buffer;
  T t_value = value;

  // Allocate space and copy the buffer to the device
  cudaMalloc( (void **)&(d_buffer), sizeof(T));
  cudaMemcpy( d_buffer, (char*) &t_value, sizeof(T), cudaMemcpyHostToDevice );

  // call the kernel to set correct indexes
  set_element_kernel<<< 1, 1 >>>( *this, d_buffer, i, lattice->field_alloc_size() );
  cudaFree(d_buffer);
}





/// A kernel that gathers elements
template <typename T>
__global__ void gather_elements_kernel( field_storage<T> field, char *buffer, unsigned * site_index, const int n, const int field_alloc_size )
{
  int Index = threadIdx.x + blockIdx.x * blockDim.x;
  if( Index < n ) {
    ((T*) buffer)[Index] = field.get(site_index[Index], field_alloc_size);
  }
}

/// CUDA implementation of gather_elements without CUDA aware MPI
template<typename T>
void field_storage<T>::gather_elements( T * RESTRICT buffer, 
                                        const unsigned * RESTRICT index_list, int n,
                                        const lattice_struct * RESTRICT lattice) const {
  unsigned *d_site_index;
  char * d_buffer;
  
  // Copy the list of boundary site indexes to the device
  cudaMalloc( (void **)&(d_site_index), n*sizeof(unsigned));
  cudaMemcpy( d_site_index, index_list, n*sizeof(unsigned), cudaMemcpyHostToDevice );

  // Call the kernel to build the list of elements
  cudaMalloc( (void **)&(d_buffer), n*sizeof(T));
  int N_blocks = n/N_threads + 1;
  gather_elements_kernel<<< N_blocks, N_threads >>>(*this, d_buffer, d_site_index, n, lattice->field_alloc_size() );
  
  // Copy the result to the host
  cudaMemcpy( (char *) buffer, d_buffer, n*sizeof(T), cudaMemcpyDeviceToHost );

  cudaFree(d_site_index);
  cudaFree(d_buffer);
}


/// A kernel that gathers elements negated
// requires unary - 
template <typename T>
__global__ void gather_elements_negated_kernel( field_storage<T> field, char *buffer, unsigned * site_index, const int n, const int field_alloc_size )
{
  int Index = threadIdx.x + blockIdx.x * blockDim.x;
  if( Index < n ) {
    ((T*) buffer)[Index] = - field.get(site_index[Index], field_alloc_size);
  }
}

/// CUDA implementation of gather_elements_negated without CUDA aware MPI
// The only difference to the above is the kernel that is called
template<typename T>
void field_storage<T>::gather_elements_negated(T * RESTRICT buffer, 
                                        const unsigned * RESTRICT index_list, int n,
                                        const lattice_struct * RESTRICT lattice) const {
  unsigned *d_site_index;
  char * d_buffer;
  
  // Copy the list of boundary site indexes to the device
  cudaMalloc( (void **)&(d_site_index), n*sizeof(unsigned));
  cudaMemcpy( d_site_index, index_list, n*sizeof(unsigned), cudaMemcpyHostToDevice );

  // Call the kernel to build the list of elements
  cudaMalloc( (void **)&(d_buffer), n*sizeof(T));
  int N_blocks = n/N_threads + 1;
  gather_elements_negated_kernel<<< N_blocks, N_threads >>>(*this, d_buffer, d_site_index, n, lattice->field_alloc_size() );
  
  // Copy the result to the host
  cudaMemcpy( (char *) buffer, d_buffer, n*sizeof(T), cudaMemcpyDeviceToHost );

  cudaFree(d_site_index);
  cudaFree(d_buffer);
}




/// A kernel that scatters the elements
template <typename T>
__global__ void place_elements_kernel( field_storage<T> field, char *buffer, unsigned * site_index, const int n, const int field_alloc_size )
{
  int Index = threadIdx.x + blockIdx.x * blockDim.x;
  if( Index < n ) {
    field.set( ((T*) buffer)[Index], site_index[Index], field_alloc_size);
  }
}

/// CUDA implementation of place_elements without CUDA aware MPI
template<typename T>
void field_storage<T>::place_elements(T * RESTRICT buffer, 
                                      const unsigned * RESTRICT index_list, int n,
                                      const lattice_struct * RESTRICT lattice) {
  unsigned *d_site_index;
  char * d_buffer;

  // Allocate space and copy the buffer to the device
  cudaMalloc( (void **)&(d_buffer), n*sizeof(T));
  cudaMemcpy( d_buffer, (char*) buffer, n*sizeof(T), cudaMemcpyHostToDevice );

  // Copy the list of boundary site indexes to the device
  cudaMalloc( (void **)&(d_site_index), n*sizeof(unsigned));
  cudaMemcpy( d_site_index, index_list, n*sizeof(unsigned), cudaMemcpyHostToDevice );

  // Call the kernel to place the elements 
  int N_blocks = n/N_threads + 1;
  place_elements_kernel<<< N_blocks, N_threads >>>( *this, d_buffer, d_site_index, n, lattice->field_alloc_size() );

  cudaFree(d_buffer);
  cudaFree(d_site_index);
}





template<typename T>
void field_storage<T>::set_local_boundary_elements(direction dir, parity par, lattice_struct * RESTRICT lattice){}



#elif defined(HILAPP)

// Hilapp requires a dummy implementation for functions with auto return type
template<typename T> 
auto field_storage<T>::get(const int i, const int field_alloc_size) const {
  T value; return value;
}
template<typename T>
auto field_storage<T>::get_element( const int i, const lattice_struct * RESTRICT lattice) const {
  T value; return value;
}

#endif //__CUDACC__

#endif