#ifndef CUDA_BACKEND
#define CUDA_BACKEND

/* CUDA implementations */


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
  fieldbuf = (real_t*) allocate_field_mem( t_elements*sizeof(real_t) * lattice->field_alloc_size() );
  #pragma acc enter data create(fieldbuf)

}

template<typename T>
void field_storage<T>::free_field() {
  #pragma acc exit data delete(fieldbuf)
  free_field_mem((void *)fieldbuf);
  fieldbuf = nullptr;
}



/// A kernel that gathers neighbour elements for communication (using the getter)
template <typename T>
__global__ void gather_comm_elements_kernel( field_storage<T> field, char *buffer, int * site_index, const int sites, const int field_alloc_size )
{
  int Index = threadIdx.x + blockIdx.x * blockDim.x;
  if( Index < sites ) {
    ((T*) buffer)[Index] = field.get(site_index[Index], field_alloc_size);
  }
}

/// CUDA implementation of gather_comm_elements without CUDA aware MPI
/// Gathers sites at the boundary that need to be communicated to neighbours
template<typename T>
void field_storage<T>::gather_comm_elements(char * buffer, lattice_struct::comm_node_struct to_node, parity par, lattice_struct * lattice) const {
  int *site_index, *d_site_index;
  char * d_buffer;
  int sites = to_node.n_sites(par);
  
  // Copy the list of boundary site indexes to the device
  site_index = (int *)std::malloc( sites*sizeof(int) );
  cudaMalloc( (void **)&(d_site_index), sites*sizeof(int));
  for (int j=0; j<sites; j++) {
    site_index[j] = to_node.site_index(j, par);
  }
  cudaMemcpy( d_site_index, site_index, sites*sizeof(int), cudaMemcpyHostToDevice );
  std::free(site_index);

  // Call the kernel to build the list of elements
  cudaMalloc( (void **)&(d_buffer), sites*sizeof(T));
  int N_blocks = sites/N_threads + 1; 
  gather_comm_elements_kernel<<< N_blocks, N_threads >>>(*this, d_buffer, d_site_index, sites, lattice->field_alloc_size() );
  
  // Copy the result to the host
  cudaMemcpy( buffer, d_buffer, sites*sizeof(T), cudaMemcpyDeviceToHost );

  cudaFree(d_site_index);
  cudaFree(d_buffer);
}


/// A kernel that scatters the neighbour sites received from a neihbour into 
/// it's proper place (using the setter)
template <typename T>
__global__ void scatter_comm_elements_kernel( field_storage<T> field, char *buffer, const int offset, const int sites, const int field_alloc_size )
{
  int Index = threadIdx.x + blockIdx.x * blockDim.x;
  if( Index < sites ) {
    field.set( ((T*) buffer)[Index], offset + Index, field_alloc_size);
  }
}

/// CUDA implementation of gather_comm_elements without CUDA aware MPI
/// Sets the values the neighbour elements from the communication buffer 
template<typename T>
void field_storage<T>::place_comm_elements(char * buffer, lattice_struct::comm_node_struct from_node, parity par){
  char * d_buffer;
  int sites = from_node.n_sites(par);

  // Allocate space and copy the buffer to the device
  cudaMalloc( (void **)&(d_buffer), sites*sizeof(T));
  cudaMemcpy( d_buffer, buffer, sites*sizeof(T), cudaMemcpyHostToDevice );

  // Call the kernel to place the elements 
  int N_blocks = sites/N_threads + 1; 
  scatter_comm_elements_kernel<<< N_blocks, N_threads >>>( *this, d_buffer, from_node.offset(par), sites, lattice->field_alloc_size() );

  cudaFree(d_buffer);
}

template<typename T>
void field_storage<T>::set_local_boundary_elements(parity par){}

#endif