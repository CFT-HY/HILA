#ifndef COMM_MPI_H
#define COMM_MPI_H

#include "../plumbing/globals.h"
#include "../plumbing/field.h"

extern MPI_Comm mpi_comm_lat;


///***********************************************************
/// Implementations of communication routines.
///

/* Float reductions */
template <>
inline void lattice_struct::reduce_node_sum(float & value, bool distribute){
  float work;
  if(distribute) {
    MPI_Allreduce( &value, &work, 1, MPI_FLOAT, MPI_SUM, mpi_comm_lat );
    value = work;
  } else {
    MPI_Reduce( &value, &work, 1, MPI_FLOAT, MPI_SUM, 0 , mpi_comm_lat );
    if (mynode() == 0) value = work;
  }
}

template <>
inline void lattice_struct::reduce_node_product(float & value, bool distribute){
  float work;
  if(distribute) {
    MPI_Allreduce( &value, &work, 1, MPI_FLOAT, MPI_PROD, mpi_comm_lat );
    value = work;
  } else {
    MPI_Reduce( &value, &work, 1, MPI_FLOAT, MPI_PROD, 0 , mpi_comm_lat );
    if (mynode() == 0) value = work;
  }
}


/* Double precision reductions */
template <>
inline void lattice_struct::reduce_node_sum(double & value, bool distribute){
  double work;
  if(distribute) {
    MPI_Allreduce( &value, &work, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_lat );
    value = work;
  } else {
    MPI_Reduce( &value, &work, 1, MPI_DOUBLE, MPI_SUM, 0 , mpi_comm_lat );
    if (mynode() == 0) value = work;
  }
}

template <>
inline void lattice_struct::reduce_node_product(double & value, bool distribute){
  double work;
  if(distribute) {
    MPI_Allreduce( &value, &work, 1, MPI_DOUBLE, MPI_PROD, mpi_comm_lat );
    value = work;
  } else {
    MPI_Reduce( &value, &work, 1, MPI_DOUBLE, MPI_PROD, 0 , mpi_comm_lat );
    if (mynode() == 0) value = work;
  }
}



#endif //COMM_MPI_H