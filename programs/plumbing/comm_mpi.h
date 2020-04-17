#ifndef COMM_MPI_H
#define COMM_MPI_H

#include "../plumbing/globals.h"

#ifdef USE_MPI


///***********************************************************
/// Implementations of communication routines.
///

/* Integer reductions */
template <>
inline void lattice_struct::reduce_node_sum(int * value, int N, bool distribute){
  int work[N];
  if(distribute) {
    MPI_Allreduce( value, work, N, MPI_INT, MPI_SUM, mpi_comm_lat );
    for(int i=0; i<N; i++)
      value[i] = work[i];
  } else {
    MPI_Reduce( value, work, N, MPI_INT, MPI_SUM, 0 , mpi_comm_lat );
    if (mynode() == 0) for(int i=0; i<N; i++)
      value[i] = work[i];

  }
}

template <>
inline void lattice_struct::reduce_node_product(int * value, int N, bool distribute){
  int work[N];
  if(distribute) {
    MPI_Allreduce( value, work, N, MPI_INT, MPI_PROD, mpi_comm_lat );
    for(int i=0; i<N; i++)
      value[i] = work[i];
  } else {
    MPI_Reduce( value, work, N, MPI_INT, MPI_PROD, 0 , mpi_comm_lat );
    if (mynode() == 0) for(int i=0; i<N; i++)
      value[i] = work[i];
  }
}

/* Float reductions */
template <>
inline void lattice_struct::reduce_node_sum(float * value, int N, bool distribute){
  float work[N];
  if(distribute) {
    MPI_Allreduce( value, work, N, MPI_FLOAT, MPI_SUM, mpi_comm_lat );
    for(int i=0; i<N; i++)
      value[i] = work[i];
  } else {
    MPI_Reduce( value, work, N, MPI_FLOAT, MPI_SUM, 0 , mpi_comm_lat );
    if (mynode() == 0) for(int i=0; i<N; i++)
      value[i] = work[i];
  }
}

template <>
inline void lattice_struct::reduce_node_product(float * value, int N, bool distribute){
  float work[N];
  if(distribute) {
    MPI_Allreduce( value, work, N, MPI_FLOAT, MPI_PROD, mpi_comm_lat );
    for(int i=0; i<N; i++)
      value[i] = work[i];
  } else {
    MPI_Reduce( value, work, N, MPI_FLOAT, MPI_PROD, 0 , mpi_comm_lat );
    if (mynode() == 0) for(int i=0; i<N; i++)
      value[i] = work[i];
  }
}


/* Double precision reductions */
template <>
inline void lattice_struct::reduce_node_sum(double * value, int N, bool distribute){
  double work[N];
  if(distribute) {
    MPI_Allreduce( value, work, N, MPI_DOUBLE, MPI_SUM, mpi_comm_lat );
    for(int i=0; i<N; i++)
      value[i] = work[i];
  } else {
    MPI_Reduce( value, work, N, MPI_DOUBLE, MPI_SUM, 0 , mpi_comm_lat );
    if (mynode() == 0) for(int i=0; i<N; i++)
      value[i] = work[i];
  }
}

template <>
inline void lattice_struct::reduce_node_product(double * value, int N, bool distribute){
  double work[N];
  if(distribute) {
    MPI_Allreduce( value, work, N, MPI_DOUBLE, MPI_PROD, mpi_comm_lat );
    for(int i=0; i<N; i++)
      value[i] = work[i];
  } else {
    MPI_Reduce( value, work, N, MPI_DOUBLE, MPI_PROD, 0 , mpi_comm_lat );
    if (mynode() == 0) for(int i=0; i<N; i++)
      value[i] = work[i];
  }
}



#endif //USE_MPI

#endif //COMM_MPI_H