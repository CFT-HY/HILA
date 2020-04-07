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
inline void FFT_field(field<cmplx<double>> & input, field<cmplx<double>> & result){
  lattice_struct * lattice = input.fs->lattice;
  field<cmplx<double>> * read_pointer = &input; // Read from input on first time, then work in result

  // Mark changed and make sure it's allocated
  result.mark_changed(ALL);

  // Run transform in all directions
  foralldir(dir){

    size_t local_sites = lattice->local_size(dir);
    size_t sites = lattice->size(dir);
    int nnodes = sites / local_sites;
    std::vector<node_info> allnodes = lattice->nodelist();
    int myrank = lattice->node_rank();
    coordinate_vector min = allnodes[myrank].min;
    coordinate_vector size = allnodes[myrank].size;
    
    mpi_column_struct mpi_column = get_mpi_column(dir);
    MPI_Comm column_communicator = mpi_column.column_communicator;
    std::vector<int> nodelist = mpi_column.nodelist;
    int my_column_rank = mpi_column.my_column_rank;


    // Count columns on this rank
    int cols = 1;
    foralldir(d2) if(d2!=dir) cols *= lattice->local_size(d2);
    //printf(" node %d: %d columns\n", myrank, cols);

    // Buffers for sending and receiving a column
    std::vector<cmplx<double>> column(sites), send_buffer(sites);

    // Do transform in all columns
    int c=0;
    while( c < cols ) {
      coordinate_vector thiscol=min;
      int cc = c;
      foralldir(d2) if(d2!=dir) {
        thiscol[d2] += cc%size[d2];
        cc/=size[d2];
      }

      // Build a list of sites matching this column
      coordinate_vector site = thiscol;
      std::vector<unsigned> sitelist(local_sites);
      for(int i=0; i<local_sites; i++ ){
        site[dir] = min[dir] + i;
        sitelist[i] = lattice->site_index(site);
      }

      // Print initial data for the column
      //printf("rank %d, col %d %d, col rank %d, send (",myrank,c,c%nodelist.size(),my_column_rank);
      //for(int t=0;t<local_sites; t++){
      //  cmplx<double> elem = fs->payload.get(sitelist[t],lattice->field_alloc_size());
      //  int t2 = t + (c%nodelist.size())*local_sites;
      //  printf(" (%g, %g) ", elem.re, elem.im );
      //}
      //printf(")\n");

      // Collect the data on this node
      char * sendbuf = (char*) send_buffer.data()+(c%nodelist.size())*local_sites;
      read_pointer->fs->payload.gather_elements(sendbuf, sitelist, lattice);

      // Send the data from each node to rank c in the column
      MPI_Gather( sendbuf, local_sites*sizeof(cmplx<double>), MPI_BYTE, 
                  column.data(), local_sites*sizeof(cmplx<double>), MPI_BYTE,
                  c%nodelist.size(), column_communicator);

      if(my_column_rank == c%nodelist.size()){
        fftw_complex *in, *out;
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sites);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sites);
        for(int t=0;t<sites; t++){
          in[t][0] = column[t].re;
          in[t][1] = column[t].im;
        }

        fftw_plan plan = fftw_plan_dft_1d( sites, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(plan);

        for(int t=0;t<sites; t++){
          column[t].re = out[t][0];
          column[t].im = out[t][1];
        }

        fftw_destroy_plan(plan);
        fftw_free(in); fftw_free(out);
      }


      MPI_Scatter( column.data(), local_sites*sizeof(cmplx<double>), MPI_BYTE, 
                  sendbuf, local_sites*sizeof(cmplx<double>), MPI_BYTE,
                  c%nodelist.size(), column_communicator);
      result.fs->payload.place_elements(sendbuf, sitelist, lattice);


      // Print result
      //printf("rank %d, col %d %d, col rank %d, recv (",myrank,c,c%nodelist.size(),my_column_rank);
      //for(int t=0;t<local_sites; t++){
      //  cmplx<double> elem = fs->payload.get(sitelist[t],lattice->field_alloc_size());
      //  int t2 = t + (c%nodelist.size())*local_sites;
      //  printf(" (%g, %g) ", elem.re, elem.im );
      //}
      //printf(")\n");
      c++;
    }

    read_pointer = &result; // From now on we work in result
  }

}


#else

inline void FFT_field(field<cmplx<double>> & input, field<cmplx<double>> & result){
}



#endif

template<>
inline void field<cmplx<double>>::FFT(){
  FFT_field(*this, *this);
}


#endif