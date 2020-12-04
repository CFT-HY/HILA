#ifndef FFT_H
#define FFT_H

#include "plumbing/field.h"
#include "datatypes/cmplx.h"
#include "fftw3.h"



#ifdef USE_MPI
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
  char * mpi_send_buffer = (char*) malloc(local_volume*sizeof(T));
  char * mpi_recv_buffer = (char*) malloc(local_volume*sizeof(T));

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

    // FFTW buffers
    fftw_complex *in, *out;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * column_size);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * column_size);
    fftw_plan plan = fftw_plan_dft_1d( column_size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);


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
        read_pointer->fs->payload.gather_elements((T*)(sendbuf + col_size*l), sitelist[c].data(), sitelist[c].size(),  lattice);
      }
      MPI_Gather( sendbuf, cpn*node_column_size*sizeof(T), MPI_BYTE, 
                  mpi_recv_buffer, cpn*node_column_size*sizeof(T), MPI_BYTE,
                  r, column_communicator);
    }
    
    // now that we have columns, run FFT on each
    for( int l=0; l<cpn; l++ ) { // Columns
      for( int e=0; e<elements; e++ ){ // Complex elements / field element

        for(int s=0; s<nnodes; s++){ // Cycle over sender nodes to collect the data
          complex_type * field_elem = (complex_type*)(mpi_recv_buffer + block_size*s + col_size*l);
          for(int t=0;t<node_column_size; t++){
            in[t+node_column_size*s][0] = field_elem[e+elements*t].re;
            in[t+node_column_size*s][1] = field_elem[e+elements*t].im;
          }
        }
        
        // Run the fft
        fftw_execute(plan);

        // Put the transformed data back in place
        for(int s=0; s<nnodes; s++){
          complex_type * field_elem = (complex_type*)(mpi_recv_buffer + block_size*s + col_size*l);
          for(int t=0;t<node_column_size; t++){
            field_elem[e+elements*t].re = out[t+node_column_size*s][0];
            field_elem[e+elements*t].im = out[t+node_column_size*s][1];
          }
        }

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
        result.fs->payload.place_elements((T*)(sendbuf + col_size*l), sitelist[c].data(), sitelist[c].size(), lattice);
      }
    }

    read_pointer = &result; // From now on we work in result
    fftw_destroy_plan(plan);
    fftw_free(in); fftw_free(out);
  }

  free(mpi_send_buffer);
  free(mpi_recv_buffer);
}


#else


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