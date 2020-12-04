#ifndef FFT_H
#define FFT_H

#include "plumbing/field.h"
#include "datatypes/cmplx.h"

#ifdef __CUDACC__
#ifdef USE_MPI
#include <cuda.h>
#include <cuda_runtime.h>
#include "cufft.h"


/// Gather one element column from the mpi buffer
template <typename complex_type>
__global__ void gather_column( cufftDoubleComplex *data, complex_type* mpi_recv_buffer, int elements, int nnodes, int node_column_size, int block_size, int e )
{
  int t = threadIdx.x + blockIdx.x * blockDim.x;
  if( t < node_column_size ) {
    for(int s=0; s<nnodes; s++){ // Cycle over sender nodes to collect the data
      complex_type * field_elem = (complex_type*)(mpi_recv_buffer + block_size*s);
      data[t+node_column_size*s].x = field_elem[e+elements*t].re;
      data[t+node_column_size*s].y = field_elem[e+elements*t].im;
    }
  }
}

/// Place results in the MPI buffer
template <typename complex_type>
__global__ void scatter_column( cufftDoubleComplex *data, complex_type* mpi_recv_buffer, int elements, int nnodes, int node_column_size, int block_size, int e )
{
  // Put the transformed data back in place
  int t = threadIdx.x + blockIdx.x * blockDim.x;
  if( t < node_column_size ) {
    for(int s=0; s<nnodes; s++){
      complex_type * field_elem = (complex_type*)(mpi_recv_buffer + block_size*s);
      field_elem[e+elements*t].re = data[t+node_column_size*s].x;
      field_elem[e+elements*t].im = data[t+node_column_size*s].y;
    }
  }
}


/// Run Fast Fourier Transform on the field to each direction
// This is done by collecting a column of elements to each node,
// running the Fourier transform on the column and redistributing
// the result
// Input and result are passed by reference. They may be the same.
// The field must be complex and the underlying complex type is supplied
// by the complex_type template argument
template<typename T, typename complex_type>
inline void FFT_field_complex(field<T> & input, field<T> & result){

  lattice_struct * lattice = input.fs->lattice;
  field<T> * read_pointer = &input; // Read from input on first time, then work in result
  size_t local_volume = lattice->local_volume();
  int elements = sizeof(T)/sizeof(complex_type);

  // Make store the result is allocated and mark it changed 
  result.check_alloc();
  result.mark_changed(ALL);

  // Allocate buffers for the MPI
  char * mpi_send_buffer, * mpi_recv_buffer;
  cudaMalloc( (void **)&(mpi_send_buffer), local_volume*sizeof(T));
  cudaMalloc( (void **)&(mpi_recv_buffer), local_volume*sizeof(T));

  // Run transform in all directions
  foralldir(dir){
    // Get the number of sites per column on this node and on all nodes
    size_t node_column_size = lattice->local_size(dir);
    size_t column_size = lattice->size(dir);
    // Count columns on this rank
    int cols = 1;
    foralldir(d2) if(d2!=dir) cols *= lattice->local_size(d2);

    // Get the MPI column in this direction
    lattice_struct::mpi_column_struct mpi_column = lattice->get_mpi_column(dir);
    MPI_Comm column_communicator = mpi_column.column_communicator;
    std::vector<int> nodelist = mpi_column.nodelist;
    int my_column_rank = mpi_column.my_column_rank;
    int nnodes = nodelist.size(); // Nodes in this column
    int cpn = cols/nnodes; // Columns per node
    // Amount of data that is communicated from a single node to another
    int block_size = cpn*node_column_size*sizeof(T);
    // Size of a full column
    int col_size = node_column_size*sizeof(T);

    // Variables needed for constructing the columns of sites
    std::vector<node_info> allnodes = lattice->nodelist();
    coordinate_vector min = allnodes[lattice->node_rank()].min;
    coordinate_vector size = allnodes[lattice->node_rank()].size;

    // CUFFT buffers
    cufftHandle plan;
    cufftDoubleComplex *data;
    int BATCH=1;
    cudaMalloc((void**)&data, sizeof(cufftDoubleComplex)*column_size*BATCH);
    cufftPlan1d(&plan, column_size, CUFFT_Z2Z, BATCH); //Z2Z for double, C2C for float


    // Construct lists of sites for each column
    std::vector<std::vector<unsigned>> sitelist(cols);
    for( int c=0; c<cols; c++ ) {
      // Get the index in the other directions
      int cc = c;
      coordinate_vector thiscol=min;
      foralldir(d2) if(d2!=dir) {
        thiscol[d2] += cc%size[d2];
        cc/=size[d2];
      }

      // And build the list of sites
      sitelist[c].resize(node_column_size);
      for(int i=0; i<node_column_size; i++ ){
        thiscol[dir] = min[dir] + i;
        sitelist[c][i] = lattice->site_index(thiscol);
      }
    }


    // Gather a number of columns to each node
    for( int r=0; r<nnodes; r++ ){
      char * sendbuf = mpi_send_buffer + block_size*r;
      for( int l=0; l<cpn; l++ ) {
        int c = r*cpn+l;
        unsigned *d_site_index;
        int n = sitelist[c].size();
        cudaMalloc( (void **)&(d_site_index), n*sizeof(unsigned));
        cudaMemcpy( d_site_index, sitelist[c].data(), n*sizeof(unsigned), cudaMemcpyHostToDevice );

        int N_blocks = n/N_threads + 1;
        gather_elements_kernel<<< N_blocks, N_threads >>>(read_pointer->fs->payload, (T*)(sendbuf + col_size*l), d_site_index, n, lattice->field_alloc_size() );
        
        cudaFree(d_site_index);
      }
      MPI_Gather( sendbuf, cpn*node_column_size*sizeof(T), MPI_BYTE, 
                  mpi_recv_buffer, cpn*node_column_size*sizeof(T), MPI_BYTE,
                  r, column_communicator);
    }
    
    // now that we have columns, run FFT on each
    for( int l=0; l<cpn; l++ ) { // Columns
      for( int e=0; e<elements; e++ ){ // Complex elements / field element
        int N_blocks = node_column_size/N_threads + 1;
        gather_column<<< N_blocks, N_threads >>>( data, 
          (complex_type*) (mpi_recv_buffer + col_size*l),
          elements, nnodes, node_column_size, block_size, e );

        // Run the fft
        cufftExecZ2Z(plan, data, data, CUFFT_FORWARD);

        scatter_column<<< N_blocks, N_threads >>>( data, 
          (complex_type*) (mpi_recv_buffer + col_size*l),
          elements, nnodes, node_column_size, block_size, e );
      }
    }

    // Now reverse the gather operation. After this each node will have its original local sites
    for( int s=0; s<nnodes; s++ ){
      char * sendbuf = mpi_send_buffer + block_size*s;
      MPI_Scatter( mpi_recv_buffer, cpn*node_column_size*sizeof(T), MPI_BYTE, 
                   sendbuf, cpn*node_column_size*sizeof(T), MPI_BYTE,
                   s, column_communicator);

      // Place the new data into field memory
      for( int l=0; l<cpn; l++ ) {
        int c = s*cpn+l;
        unsigned *d_site_index;
        int n = sitelist[c].size();
        cudaMalloc( (void **)&(d_site_index), n*sizeof(unsigned));
        cudaMemcpy( d_site_index, sitelist[c].data(), n*sizeof(unsigned), cudaMemcpyHostToDevice );

        int N_blocks = n/N_threads + 1;
        place_elements_kernel<<< N_blocks, N_threads >>>(read_pointer->fs->payload, (T*)(sendbuf + col_size*l), d_site_index, n, lattice->field_alloc_size() );
        
        cudaFree(d_site_index);
      }
    }

    read_pointer = &result; // From now on we work in result
    cufftDestroy(plan);
    cudaFree(data);
  }

  cudaFree(mpi_send_buffer);
  cudaFree(mpi_recv_buffer);
}

#else

// No MPI
template<typename T, typename C>
inline void FFT_field_complex(field<T> & input, field<T> & result){}

#endif

#else
// Not CUDA compiler
template<typename T, typename C>
inline void FFT_field_complex(field<T> & input, field<T> & result){}


#endif


template<>
inline void field<cmplx<double>>::FFT(){
  FFT_field_complex<cmplx<double>,cmplx<double>>(*this, *this);
}




/// Match a given type T to it's underlying complex type
template<typename T, class Enable = void>
struct complex_base{};

// Match to a complex type
template<>
struct complex_base<cmplx<float>>{
  using type = cmplx<float>;
};

template<>
struct complex_base<cmplx<double>>{
  using type = cmplx<double>;
};

// Match templated class B
template<template<typename B> class C, typename B>
struct complex_base<C<B>>{
  using type = typename complex_base<B>::type;
};

template<template<int a, typename B> class C, int a, typename B>
struct complex_base<C<a, B>>{
  using type = typename complex_base<B>::type;
};

template<template<int a,int b,typename B> class C, int a, int b, typename B>
struct complex_base<C<a,b,B>>{
  using type = typename complex_base<B>::type;
};



/// Run fourier transform on a complex field
// Called with any type T with a cmplx type nested in the lowest level
template<typename T>
void FFT_field(field<T> & input, field<T> & result){
  FFT_field_complex<T,typename complex_base<T>::type>(input, result);
}


#endif