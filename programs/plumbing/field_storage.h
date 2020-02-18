#ifndef FIELD_STORAGEH
#define FIELD_STORAGEH

// Pointer to field data and accessors. Only this is passed to the
// CUDA kernels and other accelerators and it only contains a minimal
// amount of data.
#ifndef VECTORIZED

#ifdef layout_SOA
template <typename T>
class field_storage {
  public:
    // Structure of arrays implementation
    constexpr static int t_elements = sizeof(T) / sizeof(real_t);
    real_t * fieldbuf = NULL;

    void allocate_field( const int field_alloc_size ) {
      fieldbuf = (real_t*) allocate_field_mem( t_elements*sizeof(real_t) * field_alloc_size );
      #pragma acc enter data create(fieldbuf)
    }

    void free_field() {
      #pragma acc exit data delete(fieldbuf)
      free_field_mem((void *)fieldbuf);
      fieldbuf = nullptr;
    }

    /// Get a single element in a field
    /// With CUDA this only works in a loop
    #pragma transformer loop_function
    inline T get(const int idx, const int field_alloc_size) const {
      assert( idx < field_alloc_size);
      T value;
      real_t *value_f = static_cast<real_t *>(static_cast<void *>(&value));
      for (int i=0; i<(sizeof(T)/sizeof(real_t)); i++) {
         value_f[i] = fieldbuf[i*field_alloc_size + idx];
      }
      return value; 
    }

    /// Set a single element in a field
    /// With CUDA this only works in a loop  
    #pragma transformer loop_function
    inline void set(T value, const int idx, const int field_alloc_size) {
      assert( idx < field_alloc_size);
      real_t *value_f = static_cast<real_t *>(static_cast<void *>(&value));
      for (int i=0; i<(sizeof(T)/sizeof(real_t)); i++) {
        fieldbuf[i*field_alloc_size + idx] = value_f[i];
      }
    }
};

#else

template <typename T>
class field_storage {
  public:

    // Array of structures implementation
    T * fieldbuf = NULL;

    void allocate_field( const int field_alloc_size ) {
      fieldbuf = (T*) allocate_field_mem( sizeof(T) * field_alloc_size);
      #pragma acc enter data create(fieldbuf)
    }

    void free_field() {
      #pragma acc exit data delete(fieldbuf)
      free_field_mem((void *)fieldbuf);
      fieldbuf = nullptr;
    }

    #pragma transformer loop_function
     inline T get(const int i, const int field_alloc_size) const
    {
      // There is some problem with directly assigning intrinsic vectors, at least.
      // This is a universal workaround, but could be fixed by assigning element
      // by element
      T value = fieldbuf[i];
      return value;
    }

    #pragma transformer loop_function
    inline void set(const T &value, const int i, const int field_alloc_size) 
    {
      fieldbuf[i] = value;
    }

    #pragma transformer loop_function
    inline void set(const T &value, unsigned int i) 
    {
      fieldbuf[i] = value;
    }

    void gather_comm_elements(char * buffer, lattice_struct::comm_node_struct to_node, parity par) const;
    /// Place boundary elements from neighbour
    void place_comm_elements(char * buffer, lattice_struct::comm_node_struct from_node, parity par);
    /// Place boundary elements from local lattice (used in vectorized version)
    void set_local_boundary_elements(parity par);
};

#endif
#endif



#ifdef VECTORIZED

#ifdef USE_MPI
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

#endif

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



#elif CUDA


/* CUDA implementations */

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
void field_storage<T>::gather_comm_elements(char * buffer, lattice_struct::comm_node_struct to_node, parity par) const {
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


#else


/// The standard (Non-CUDA) implementation of gather_comm_elements
/* Gathers sites at the boundary that need to be communicated to neighbours */
template<typename T>
void field_storage<T>::gather_comm_elements(char * buffer, lattice_struct::comm_node_struct to_node, parity par) const {
  for (int j=0; j<to_node.n_sites(par); j++) {
    T element = get(to_node.site_index(j, par), lattice->field_alloc_size());
    std::memcpy( buffer + j*sizeof(T), (char *) (&element), sizeof(T) );
  }
}
/// The standard (Non-CUDA) implementation of place_comm_elements
/* Sets the values the neighbour elements from the communication buffer */
template<typename T>
void field_storage<T>::place_comm_elements(char * buffer, lattice_struct::comm_node_struct from_node, parity par){
  for (int j=0; j<from_node.n_sites(par); j++) {
    T element = *((T*) ( buffer + j*sizeof(T) ));
    set(element, from_node.offset(par)+j);
  }
}

template<typename T>
void field_storage<T>::set_local_boundary_elements(parity par){}


#endif


#endif //field storage