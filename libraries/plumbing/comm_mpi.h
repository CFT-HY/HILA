#ifndef COMM_MPI_H
#define COMM_MPI_H

#include "plumbing/globals.h"

#ifdef USE_MPI


bool is_comm_initialized(void);

///***********************************************************
/// Implementations of communication routines.
///

/* Integer reductions */
template <>
inline void lattice_struct::reduce_node_sum(int * value, int N, bool distribute){
  int work[N];
  if(distribute) {
    MPI_Allreduce( value, work, N, MPI_INT, MPI_SUM, lattice->mpi_comm_lat );
    for(int i=0; i<N; i++)
      value[i] = work[i];
  } else {
    MPI_Reduce( value, work, N, MPI_INT, MPI_SUM, 0 , lattice->mpi_comm_lat );
    if (mynode() == 0) for(int i=0; i<N; i++)
      value[i] = work[i];

  }
}

template <>
inline void lattice_struct::reduce_node_product(int * value, int N, bool distribute){
  int work[N];
  if(distribute) {
    MPI_Allreduce( value, work, N, MPI_INT, MPI_PROD, lattice->mpi_comm_lat );
    for(int i=0; i<N; i++)
      value[i] = work[i];
  } else {
    MPI_Reduce( value, work, N, MPI_INT, MPI_PROD, 0 , lattice->mpi_comm_lat );
    if (mynode() == 0) for(int i=0; i<N; i++)
      value[i] = work[i];
  }
}

/* Float reductions */
template <>
inline void lattice_struct::reduce_node_sum(float * value, int N, bool distribute){
  float work[N];
  if(distribute) {
    MPI_Allreduce( value, work, N, MPI_FLOAT, MPI_SUM, lattice->mpi_comm_lat );
    for(int i=0; i<N; i++)
      value[i] = work[i];
  } else {
    MPI_Reduce( value, work, N, MPI_FLOAT, MPI_SUM, 0 , lattice->mpi_comm_lat );
    if (mynode() == 0) for(int i=0; i<N; i++)
      value[i] = work[i];
  }
}

template <>
inline void lattice_struct::reduce_node_product(float * value, int N, bool distribute){
  float work[N];
  if(distribute) {
    MPI_Allreduce( value, work, N, MPI_FLOAT, MPI_PROD, lattice->mpi_comm_lat );
    for(int i=0; i<N; i++)
      value[i] = work[i];
  } else {
    MPI_Reduce( value, work, N, MPI_FLOAT, MPI_PROD, 0 , lattice->mpi_comm_lat );
    if (mynode() == 0) for(int i=0; i<N; i++)
      value[i] = work[i];
  }
}


/* Double precision reductions */
template <>
inline void lattice_struct::reduce_node_sum(double * value, int N, bool distribute){
  double work[N];
  if(distribute) {
    MPI_Allreduce( value, work, N, MPI_DOUBLE, MPI_SUM, lattice->mpi_comm_lat );
    for(int i=0; i<N; i++)
      value[i] = work[i];
  } else {
    MPI_Reduce( value, work, N, MPI_DOUBLE, MPI_SUM, 0 , lattice->mpi_comm_lat );
    if (mynode() == 0) for(int i=0; i<N; i++)
      value[i] = work[i];
  }
}

template <>
inline void lattice_struct::reduce_node_product(double * value, int N, bool distribute){
  double work[N];
  if(distribute) {
    MPI_Allreduce( value, work, N, MPI_DOUBLE, MPI_PROD, lattice->mpi_comm_lat );
    for(int i=0; i<N; i++)
      value[i] = work[i];
  } else {
    MPI_Reduce( value, work, N, MPI_DOUBLE, MPI_PROD, 0 , lattice->mpi_comm_lat );
    if (mynode() == 0) for(int i=0; i<N; i++)
      value[i] = work[i];
  }
}

template <typename T>
void broadcast(T & var) {
  MPI_Bcast(&var, sizeof(T), MPI_BYTE, 0, lattice->mpi_comm_lat);
} 

template <typename T>
void broadcast(T * var) {
  static_assert(sizeof(T) > 0 && "Do not use pointers to broadcast()-function");
}

template <typename T>
void broadcast_array(T * var, int n) {
  MPI_Bcast(var, sizeof(T)*n, MPI_BYTE, 0, lattice->mpi_comm_lat);
}

#endif //USE_MPI

#endif //COMM_MPI_H