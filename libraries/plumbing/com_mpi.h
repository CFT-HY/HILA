#ifndef COMM_MPI_H
#define COMM_MPI_H

#include "plumbing/globals.h"

#ifdef USE_MPI

#include "plumbing/lattice.h"

bool is_comm_initialized(void);

extern timer start_send_timer, 
             wait_send_timer, 
             post_receive_timer,
             wait_receive_timer,
             synchronize_timer,
             reduction_timer,
             broadcast_timer,
             send_timer,
             cancel_send_timer,
             cancel_receive_timer,
             sublattice_sync_timer;


///***********************************************************
/// Implementations of communication routines.
///

// broadcast templates

template <typename T>
void broadcast(T & var) {
  static_assert( std::is_trivial<T>::value, "broadcast(var) must use trivial type");
  broadcast_timer.start();
  MPI_Bcast(&var, sizeof(T), MPI_BYTE, 0, lattice->mpi_comm_lat);
  broadcast_timer.stop();
}

template <typename T>
void broadcast(std::vector<T> & list) {
  broadcast_timer.start();

  static_assert( std::is_trivial<T>::value, "broadcast(vector<T>) must have trivial T" ); 

  int size = list.size();
  MPI_Bcast(&size, sizeof(int), MPI_BYTE, 0, lattice->mpi_comm_lat);
  if (hila::myrank() != 0) {
    list.resize(size);
  }

  // move vectors directly to the storage
  MPI_Bcast((void *)list.data(), sizeof(T)*size, MPI_BYTE, 0, lattice->mpi_comm_lat);

  broadcast_timer.stop();
}

template <typename T>
void broadcast(T * var) {
  static_assert(sizeof(T) > 0 && "Do not use pointers to broadcast()-function");
}

template <typename T>
void broadcast_array(T * var, int n) {
  broadcast_timer.start();
  MPI_Bcast((void *)var, sizeof(T)*n, MPI_BYTE, 0, lattice->mpi_comm_lat);
  broadcast_timer.stop();
}

// DO string bcasts separately
void broadcast(std::string & r);
void broadcast(std::vector<std::string> &l);

// Reduction templates
// TODO: implement using custom MPI Ops!  Then no ambiguity

template <typename T>
void lattice_struct::reduce_node_sum(T * value, int N, bool distribute) {
  T work[N];
  MPI_Datatype dtype;

  if constexpr (std::is_same<number_type<T>,int>::value) {
    dtype = MPI_INT;
  } else if constexpr (std::is_same<number_type<T>,float>::value) {
    dtype = MPI_FLOAT;
  } else if constexpr (std::is_same<number_type<T>,double>::value) {
    dtype = MPI_DOUBLE;
  } else {
    static_assert(sizeof(T)<0 && "Unknown number_type in reduce_node_sum");
  }

  reduction_timer.start();
  if(distribute) {
    MPI_Allreduce( (void *)value, (void *) work, N*sizeof(T)/sizeof(number_type<T>), 
                   dtype, MPI_SUM, lattice->mpi_comm_lat );
    for(int i=0; i<N; i++)
      value[i] = work[i];
  } else {
    MPI_Reduce( (void *)value, (void *)work, N*sizeof(T)/sizeof(number_type<T>), dtype,
                MPI_SUM, 0 , lattice->mpi_comm_lat );
    if (hila::myrank() == 0) for(int i=0; i<N; i++)
      value[i] = work[i];

  }
  reduction_timer.stop();
}

// Product reduction template - so far only for int, float, dbl

template <typename T>
void lattice_struct::reduce_node_product(T * value, int N, bool distribute) {
  T work[N];
  MPI_Datatype dtype;

  if constexpr (std::is_same<T,int>::value) {
    dtype = MPI_INT;
  } else if constexpr (std::is_same<T,float>::value) {
    dtype = MPI_FLOAT;
  } else if constexpr (std::is_same<T,double>::value) {
    dtype = MPI_DOUBLE;
  } else {
    static_assert(sizeof(T)<0 && "Unknown number_type in reduce_node_product");
  }

  reduction_timer.start();
  if(distribute) {
    MPI_Allreduce( (void *)value, (void *) work, N, dtype, MPI_PROD, lattice->mpi_comm_lat );
    for(int i=0; i<N; i++)
      value[i] = work[i];
  } else {
    MPI_Reduce( (void *)value, (void *)work, N, dtype, MPI_PROD, 0 , lattice->mpi_comm_lat );
    if (hila::myrank() == 0) for(int i=0; i<N; i++)
      value[i] = work[i];

  }
  reduction_timer.stop();
}




#endif //USE_MPI

#endif //COMM_MPI_H