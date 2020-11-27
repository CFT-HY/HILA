#ifndef FFT_H
#define FFT_H

#include "plumbing/field.h"
#include "datatypes/cmplx.h"
#include "fftw3.h"



#ifdef USE_MPI

struct mpi_column_struct {
  std::vector<int> nodelist;
  MPI_Comm column_communicator = nullptr;
  int my_column_rank;
  bool init = true;
};


/// Return a communicator and a list of nodes in this column
static mpi_column_struct get_mpi_column(direction dir){
  static mpi_column_struct mpi_column[NDIM];
  
  if( mpi_column[dir].init ) {
    // The communicator does not exist yet, so build it now

    std::vector<node_info> allnodes = lattice->nodelist();
    int myrank = lattice->node_rank();
    
    // Build a list of nodes in this column
    // All nodes should have these in the same order
    coordinate_vector min = allnodes[myrank].min;
    coordinate_vector size = allnodes[myrank].size;
    foralldir(d2){
      assert( min[d2] == lattice->min_coordinate()[d2] );
    }
    for( int rank=0; rank < allnodes.size(); rank++ ) {
      node_info node = allnodes[rank];
      bool in_column = true;
      foralldir(d2) if( d2 != dir && node.min[d2] != min[d2] ){
        in_column = false;
      }
      if( in_column ){
        if(rank==myrank)
          mpi_column[dir].my_column_rank = mpi_column[dir].nodelist.size();
        mpi_column[dir].nodelist.push_back(rank);
        //printf(" node %d: %d in my column\n", myrank, rank);
      }
    }

    MPI_Comm_split( MPI_COMM_WORLD, mpi_column[dir].nodelist[0], mpi_column[dir].my_column_rank, &mpi_column[dir].column_communicator );

    mpi_column[dir].init = false;

  }

  return mpi_column[dir];
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

  // Run transform in all directions
  foralldir(dir){
    // Get the number of sites per column on this node and on all nodes
    size_t sites_per_col = lattice->local_size(dir);
    size_t sites = lattice->size(dir);

    // Get the MPI column in this direction
    mpi_column_struct mpi_column = get_mpi_column(dir);
    MPI_Comm column_communicator = mpi_column.column_communicator;
    std::vector<int> nodelist = mpi_column.nodelist;
    int my_column_rank = mpi_column.my_column_rank;

    // Variables needed for constructing the columns of sites
    std::vector<node_info> allnodes = lattice->nodelist();
    int myrank = lattice->node_rank();
    int nnodes = nodelist.size();
    coordinate_vector min = allnodes[myrank].min;
    coordinate_vector size = allnodes[myrank].size;

    // FFTW buffers
    fftw_complex *in, *out;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sites);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sites);
    fftw_plan plan = fftw_plan_dft_1d( sites, in, out, FFTW_FORWARD, FFTW_ESTIMATE);


    // Count columns on this rank
    int cols = 1;
    foralldir(d2) if(d2!=dir) cols *= lattice->local_size(d2);
    int cpn = cols/nnodes;

    // Buffers for sending and receiving a column
    std::vector<std::vector<complex_type>> column(cols), send_buffer(cols);

    // Do transform in all columns
    int c=0;
    std::vector<std::vector<unsigned>> sitelist(cols);
    for( c=0; c<cols; c++ ) {
      column[c].resize(elements*sites);
      send_buffer[c].resize(elements*sites_per_col);

      // Build a column for each node and send the data
      int cc = c;
      int l = c/nnodes; // The column index on that node
      coordinate_vector thiscol=min;

      foralldir(d2) if(d2!=dir) {
        thiscol[d2] += cc%size[d2];
        cc/=size[d2];
      }

      // Build a list of sites matching this column
      sitelist[c].resize(sites_per_col);
      for(int i=0; i<sites_per_col; i++ ){
        thiscol[dir] = min[dir] + i;
        sitelist[c][i] = lattice->site_index(thiscol);
      }
    }

    char * mpi_send_buffer = (char*) malloc(local_volume*sizeof(T));
    char * mpi_recv_buffer = (char*) malloc(local_volume*sizeof(T));
    int block_size = cpn*sites_per_col*sizeof(T);
    int col_size = sites_per_col*sizeof(T);
    for( int r=0; r<nnodes; r++ ){
      char * sendbuf = mpi_send_buffer + block_size*r;
      for( int l=0; l<cpn; l++ ) {
        int c = r*cpn+l;
        read_pointer->fs->payload.gather_elements((T*)(sendbuf + col_size*l), sitelist[c].data(), sitelist[c].size(),  lattice);
      }
      MPI_Gather( sendbuf, cpn*sites_per_col*sizeof(T), MPI_BYTE, 
                  mpi_recv_buffer, cpn*sites_per_col*sizeof(T), MPI_BYTE,
                  r, column_communicator);
    }
    
    for( int l=0; l<cpn; l++ ) {
      int c = my_column_rank*cpn+l;
      // Run the FFT on my column

      for( int e=0; e<elements; e++ ){
        for(int s=0; s<nnodes; s++){
          complex_type * field_elem = (complex_type*)(mpi_recv_buffer + block_size*s + col_size*l);
          for(int t=0;t<sites_per_col; t++){
            in[t+sites_per_col*s][0] = field_elem[e+elements*t].re;
            in[t+sites_per_col*s][1] = field_elem[e+elements*t].im;
          }
        }

        fftw_execute(plan);

        for(int t=0;t<sites; t++){
          for(int s=0; s<nnodes; s++){
            complex_type * field_elem = (complex_type*)(mpi_recv_buffer + block_size*s + col_size*l);
            for(int t=0;t<sites_per_col; t++){
              field_elem[e+elements*t].re = out[t+sites_per_col*s][0];
              field_elem[e+elements*t].im = out[t+sites_per_col*s][1];
            }
          }
        }
      }
    }

    for( int s=0; s<nnodes; s++ ){
      char * sendbuf = mpi_send_buffer + block_size*s;
      MPI_Scatter( mpi_recv_buffer, cpn*sites_per_col*sizeof(T), MPI_BYTE, 
                   sendbuf, cpn*sites_per_col*sizeof(T), MPI_BYTE,
                   s, column_communicator);

      for( int l=0; l<cpn; l++ ) {
        int c = s*cpn+l;
        result.fs->payload.place_elements((T*)(sendbuf + col_size*l), sitelist[c].data(), sitelist[c].size(), lattice);
      }
    }

    read_pointer = &result; // From now on we work in result
    fftw_destroy_plan(plan);
    fftw_free(in); fftw_free(out);
  }

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