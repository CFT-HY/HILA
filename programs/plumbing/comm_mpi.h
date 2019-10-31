#ifndef COMM_MPI_H
#define COMM_MPI_H

#include "../plumbing/globals.h"
#include "../plumbing/field.h"

///***********************************************************
/// Implementations of communication routines.
///


template <typename T>
void lattice_struct::reduce_node_sum(T value, bool distribute){}

template <typename T>
void lattice_struct::reduce_node_product(T value, bool distribute){}



#endif //COMM_MPI_H