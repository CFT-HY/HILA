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

/* Integer reductions */
template <>
inline void lattice_struct::reduce_node_sum(int * value, int N, bool distribute){
  int work[N];
  reduction_timer.start();
  if(distribute) {
    MPI_Allreduce( value, work, N, MPI_INT, MPI_SUM, lattice->mpi_comm_lat );
    for(int i=0; i<N; i++)
      value[i] = work[i];
  } else {
    MPI_Reduce( value, work, N, MPI_INT, MPI_SUM, 0 , lattice->mpi_comm_lat );
    if (mynode() == 0) for(int i=0; i<N; i++)
      value[i] = work[i];

  }
  reduction_timer.end();
}

template <>
inline void lattice_struct::reduce_node_product(int * value, int N, bool distribute){
  int work[N];
  reduction_timer.start();
  if(distribute) {
    MPI_Allreduce( value, work, N, MPI_INT, MPI_PROD, lattice->mpi_comm_lat );
    for(int i=0; i<N; i++)
      value[i] = work[i];
  } else {
    MPI_Reduce( value, work, N, MPI_INT, MPI_PROD, 0 , lattice->mpi_comm_lat );
    if (mynode() == 0) for(int i=0; i<N; i++)
      value[i] = work[i];
  }
  reduction_timer.end();
}

/* Float reductions */
template <>
inline void lattice_struct::reduce_node_sum(float * value, int N, bool distribute){
  float work[N];
  reduction_timer.start();
  if(distribute) {
    MPI_Allreduce( value, work, N, MPI_FLOAT, MPI_SUM, lattice->mpi_comm_lat );
    for(int i=0; i<N; i++)
      value[i] = work[i];
  } else {
    MPI_Reduce( value, work, N, MPI_FLOAT, MPI_SUM, 0 , lattice->mpi_comm_lat );
    if (mynode() == 0) for(int i=0; i<N; i++)
      value[i] = work[i];
  }
  reduction_timer.end();
}

template <>
inline void lattice_struct::reduce_node_product(float * value, int N, bool distribute){
  float work[N];
  reduction_timer.start();
  if(distribute) {
    MPI_Allreduce( value, work, N, MPI_FLOAT, MPI_PROD, lattice->mpi_comm_lat );
    for(int i=0; i<N; i++)
      value[i] = work[i];
  } else {
    MPI_Reduce( value, work, N, MPI_FLOAT, MPI_PROD, 0 , lattice->mpi_comm_lat );
    if (mynode() == 0) for(int i=0; i<N; i++)
      value[i] = work[i];
  }
  reduction_timer.end();
}


/* Double precision reductions */
template <>
inline void lattice_struct::reduce_node_sum(double * value, int N, bool distribute){
  double work[N];
  reduction_timer.start();
  if(distribute) {
    MPI_Allreduce( value, work, N, MPI_DOUBLE, MPI_SUM, lattice->mpi_comm_lat );
    for(int i=0; i<N; i++)
      value[i] = work[i];
  } else {
    MPI_Reduce( value, work, N, MPI_DOUBLE, MPI_SUM, 0 , lattice->mpi_comm_lat );
    if (mynode() == 0) for(int i=0; i<N; i++)
      value[i] = work[i];
  }
  reduction_timer.end();
}

template <>
inline void lattice_struct::reduce_node_product(double * value, int N, bool distribute){
  double work[N];
  reduction_timer.start();
  if(distribute) {
    MPI_Allreduce( value, work, N, MPI_DOUBLE, MPI_PROD, lattice->mpi_comm_lat );
    for(int i=0; i<N; i++)
      value[i] = work[i];
  } else {
    MPI_Reduce( value, work, N, MPI_DOUBLE, MPI_PROD, 0 , lattice->mpi_comm_lat );
    if (mynode() == 0) for(int i=0; i<N; i++)
      value[i] = work[i];
  }
  reduction_timer.end();
}

template <typename T>
void broadcast(T & var) {
  broadcast_timer.start();
  MPI_Bcast(&var, sizeof(T), MPI_BYTE, 0, lattice->mpi_comm_lat);
  broadcast_timer.end();
} 

template <typename T>
void broadcast(T * var) {
  static_assert(sizeof(T) > 0 && "Do not use pointers to broadcast()-function");
}

template <typename T>
void broadcast_array(T * var, int n) {
  broadcast_timer.start();
  MPI_Bcast(var, sizeof(T)*n, MPI_BYTE, 0, lattice->mpi_comm_lat);
  broadcast_timer.end();
}

#endif //USE_MPI

#endif //COMM_MPI_H