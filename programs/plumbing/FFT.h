#ifndef FFT_H
#define FFT_H

#include "../plumbing/field.h"
#include "../datatypes/cmplx.h"
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
// The must be complex and the underlying complex type is supplied
// by the complex_type template argument
template<typename T, typename complex_type>
inline void FFT_field_complex(field<T> & input, field<T> & result){

  lattice_struct * lattice = input.fs->lattice;
  field<T> * read_pointer = &input; // Read from input on first time, then work in result

  int elements = sizeof(T)/sizeof(complex_type);

  // Mark changed and make sure it's allocated
  result.mark_changed(ALL);

  // Run transform in all directions
  foralldir(dir){
    // Get the number of sites per column on this node and on all nodes
    size_t local_sites = lattice->local_size(dir);
    size_t sites = lattice->size(dir);

    // Get the MPI column in this direction
    mpi_column_struct mpi_column = get_mpi_column(dir);
    MPI_Comm column_communicator = mpi_column.column_communicator;
    std::vector<int> nodelist = mpi_column.nodelist;
    int my_column_rank = mpi_column.my_column_rank;

    // Buffers for sending and receiving a column
    std::vector<complex_type> column(elements*sites), send_buffer(elements*sites);

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

    // Do transform in all columns
    int c=0;
    int n;
    std::vector<std::vector<unsigned>> sitelist(nnodes);
    while( c < cols ) {

      // Build a column for each node and send the data
      for(n=0; n < nnodes && c+n < cols; n++ ){
        int root = (c+n)%nnodes; // The node that does the calculation
        int cc = c+n;
        coordinate_vector thiscol=min;

        foralldir(d2) if(d2!=dir) {
          thiscol[d2] += cc%size[d2];
          cc/=size[d2];
        }

        // Build a list of sites matching this column
        sitelist[n].resize(local_sites);
        for(int i=0; i<local_sites; i++ ){
          thiscol[dir] = min[dir] + i;
          sitelist[n][i] = lattice->site_index(thiscol);
        }

        // Collect the data to node n
        char * sendbuf = (char*) send_buffer.data()+root*elements*local_sites;
        read_pointer->fs->payload.gather_elements(sendbuf, sitelist[n], lattice);
        MPI_Gather( sendbuf, local_sites*sizeof(T), MPI_BYTE, 
                   column.data(), local_sites*sizeof(T), MPI_BYTE,
                   root, column_communicator);
      }

      if( my_column_rank < n ){
        // Run the FFT on my column

        for( int e=0; e<elements; e++ ){
          for(int t=0;t<sites; t++){
            in[t][0] = column[e+elements*t].re;
            in[t][1] = column[e+elements*t].im;
          }

          fftw_execute(plan);

          for(int t=0;t<sites; t++){
            column[e+elements*t].re = out[t][0];
            column[e+elements*t].im = out[t][1];
          }
        }
      }

      for(n=0; n < nnodes && c+n < cols; n++ ){
        int root = (c+n)%nnodes; // The node that does the calculation
        char * sendbuf = (char*) send_buffer.data()+root*elements*local_sites;

        MPI_Scatter( column.data(), local_sites*sizeof(T), MPI_BYTE, 
                   sendbuf, local_sites*sizeof(T), MPI_BYTE,
                   root, column_communicator);
        result.fs->payload.place_elements(sendbuf, sitelist[n], lattice);
      }

      c+=n;
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