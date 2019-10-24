#include "../plumbing/globals.h"
#include "../plumbing/field.h"

///***********************************************************
/// Vanilla (non-mpi) implementations of communication routines.
///

template <typename T>
void lattice_struct::reduce_node_sum(T value, bool distribute){}

template <typename T>
void lattice_struct::reduce_node_product(T value, bool distribute){}

//template<typename T>
//void field<T>::start_move(direction d, parity p){}

