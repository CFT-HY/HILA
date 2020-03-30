#ifndef CUDA_BACKEND
#define CUDA_BACKEND

#include "../defs.h"
#include "../field_storage.h"


/* CUDA implementations */
template<typename T>
void field_storage<T>::allocate_field(lattice_struct * lattice) {
  constexpr static int t_elements = sizeof(T) / sizeof(real_t);
  cudaMalloc( &fieldbuf, t_elements*sizeof(real_t) * lattice->field_alloc_size() );
  check_cuda_error("Allocate field memory");
  if (fieldbuf == nullptr) {
    std::cout << "Failure in field memory allocation\n";
    exit(1);
  }
}

template<typename T>
void field_storage<T>::free_field() {
  if (fieldbuf != nullptr){
    cudaFree(fieldbuf);
    check_cuda_error("Free field memory");
  }
  fieldbuf = nullptr;
}




#ifdef __CUDACC__

template<typename T> 
__device__ __host__ auto field_storage<T>::get(const int i, const int field_alloc_size) const {
  assert( i < field_alloc_size);
  T value;
  real_t *value_f = static_cast<real_t *>(static_cast<void *>(&value));
  real_t *fp = static_cast<real_t *>(fieldbuf);
  for (int e=0; e<(sizeof(T)/sizeof(real_t)); e++) {
     value_f[e] = fp[e*field_alloc_size + i];
   }
  return value;
}


template<typename T>
template<typename A>
__device__ __host__ inline void field_storage<T>::set(const A &value, const int i, const int field_alloc_size){
  assert( i < field_alloc_size);
  real_t *value_f = (real_t *)(&value);
  real_t *fp = static_cast<real_t *>(fieldbuf);
  for (int e=0; e<(sizeof(T)/sizeof(real_t)); e++) {
    fp[e*field_alloc_size + i] = value_f[e];
  }
}

/// A kernel that gathers elements
template <typename T>
__global__ void gather_elements_kernel( field_storage<T> field, char *buffer, unsigned * site_index, const int sites, const int field_alloc_size )
{
  int Index = threadIdx.x + blockIdx.x * blockDim.x;
  if( Index < sites ) {
    ((T*) buffer)[Index] = field.get(site_index[Index], field_alloc_size);
  }
}

/// CUDA implementation of gather_elements without CUDA aware MPI
template<typename T>
void field_storage<T>::gather_elements(char * buffer, std::vector<unsigned> index_list, lattice_struct * lattice) const {
  unsigned *d_site_index;
  char * d_buffer;
  int sites = index_list.size();
  
  // Copy the list of boundary site indexes to the device
  cudaMalloc( (void **)&(d_site_index), sites*sizeof(unsigned));
  cudaMemcpy( d_site_index, index_list.data(), sites*sizeof(unsigned), cudaMemcpyHostToDevice );

  // Call the kernel to build the list of elements
  cudaMalloc( (void **)&(d_buffer), sites*sizeof(T));
  int N_blocks = sites/N_threads + 1; 
  gather_elements_kernel<<< N_blocks, N_threads >>>(*this, d_buffer, d_site_index, sites, lattice->field_alloc_size() );
  
  // Copy the result to the host
  cudaMemcpy( buffer, d_buffer, sites*sizeof(T), cudaMemcpyDeviceToHost );

  cudaFree(d_site_index);
  cudaFree(d_buffer);
}




/// A kernel that scatters the elements
template <typename T>
__global__ void scatter_elements_kernel( field_storage<T> field, char *buffer, unsigned * site_index, const int sites, const int field_alloc_size )
{
  int Index = threadIdx.x + blockIdx.x * blockDim.x;
  if( Index < sites ) {
    field.set( ((T*) buffer)[Index], site_index[Index], field_alloc_size);
  }
}

/// CUDA implementation of gather_elements without CUDA aware MPI
template<typename T>
void field_storage<T>::place_elements(char * buffer, std::vector<unsigned> index_list, lattice_struct * lattice) {
  unsigned *d_site_index;
  char * d_buffer;
  int sites = index_list.size();

  // Allocate space and copy the buffer to the device
  cudaMalloc( (void **)&(d_buffer), sites*sizeof(T));
  cudaMemcpy( d_buffer, buffer, sites*sizeof(T), cudaMemcpyHostToDevice );

  // Copy the list of boundary site indexes to the device
  cudaMalloc( (void **)&(d_site_index), sites*sizeof(unsigned));
  cudaMemcpy( d_site_index, index_list.data(), sites*sizeof(unsigned), cudaMemcpyHostToDevice );

  // Call the kernel to place the elements 
  int N_blocks = sites/N_threads + 1;
  scatter_elements_kernel<<< N_blocks, N_threads >>>( *this, d_buffer, d_site_index, sites, lattice->field_alloc_size() );

  cudaFree(d_buffer);
  cudaFree(d_site_index);
}



template<typename T>
void field_storage<T>::set_local_boundary_elements(direction dir, parity par, lattice_struct * lattice){}

#endif //__CUDACC__

#endif