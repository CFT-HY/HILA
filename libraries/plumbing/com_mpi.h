#ifndef COMM_MPI_H
#define COMM_MPI_H

#include "plumbing/globals.h"

#ifdef USE_MPI

#include "plumbing/lattice.h"

extern timer start_send_timer, wait_send_timer, post_receive_timer, wait_receive_timer,
    synchronize_timer, reduction_timer, broadcast_timer, send_timer, cancel_send_timer,
    cancel_receive_timer, sublattice_sync_timer;

///***********************************************************
/// Implementations of communication routines.
///

/// Broadcast template for standard type
template <typename T> void broadcast(T &var) {
    static_assert(std::is_trivial<T>::value, "broadcast(var) must use trivial type");
    if (hila::check_input) return;

    broadcast_timer.start();
    MPI_Bcast(&var, sizeof(T), MPI_BYTE, 0, lattice->mpi_comm_lat);
    broadcast_timer.stop();
}

/// Broadcast for std::vector
template <typename T> void broadcast(std::vector<T> &list) {

    static_assert(std::is_trivial<T>::value, "broadcast(vector<T>) must have trivial T");

    if (hila::check_input) return;

    broadcast_timer.start();

    int size = list.size();
    MPI_Bcast(&size, sizeof(int), MPI_BYTE, 0, lattice->mpi_comm_lat);
    if (hila::myrank() != 0) {
        list.resize(size);
    }

    // move vectors directly to the storage
    MPI_Bcast((void *)list.data(), sizeof(T) * size, MPI_BYTE, 0, lattice->mpi_comm_lat);

    broadcast_timer.stop();
}

template <typename T> void broadcast(T *var) {
    static_assert(sizeof(T) > 0 && "Do not use pointers to broadcast()-function");
}

/// Broadcast for arrays where size must be known and same for all nodes
template <typename T> void broadcast_array(T *var, int n) {

    if (hila::check_input) return;

    broadcast_timer.start();
    MPI_Bcast((void *)var, sizeof(T) * n, MPI_BYTE, 0, lattice->mpi_comm_lat);
    broadcast_timer.stop();
}

// DO string bcasts separately
void broadcast(std::string &r);
void broadcast(std::vector<std::string> &l);

/// and broadcast with two values
template <typename T, typename U>
void broadcast(T & t, U & u) {

    if (hila::check_input) return;

    struct { 
        T tv; 
        U uv;
    } s = {t,u};
    
    broadcast(s);
    t = s.tv;
    u = s.uv;
}

// try to get the basic data type of the message
// this is just to enable a bit larger messages
template <typename T> MPI_Datatype get_MPI_number_type(int &size) {

    if (std::is_same<number_type<T>, int>::value) {
        size = sizeof(int);
        return MPI_INT;
    } else if (std::is_same<number_type<T>, long>::value) {
        size = sizeof(long);
        return MPI_LONG;
    } else if (std::is_same<number_type<T>, float>::value) {
        size = sizeof(float);
        return MPI_FLOAT;
    } else if (std::is_same<number_type<T>, double>::value) {
        size = sizeof(double);
        return MPI_DOUBLE;
    }
    size = 1;
    return MPI_BYTE;
}

// Reduction templates
// TODO: implement using custom MPI Ops!  Then no ambiguity

template <typename T>
void lattice_struct::reduce_node_sum(T *value, int N, bool distribute) {
    T work[N];
    MPI_Datatype dtype;

    if (std::is_same<number_type<T>, int>::value) {
        dtype = MPI_INT;
    } else if (std::is_same<number_type<T>, float>::value) {
        dtype = MPI_FLOAT;
    } else if (std::is_same<number_type<T>, double>::value) {
        dtype = MPI_DOUBLE;
    } else {
        static_assert(sizeof(T) > 0, "Unknown number_type in reduce_node_sum");
    }

    reduction_timer.start();
    if (distribute) {
        MPI_Allreduce((void *)value, (void *)work, N * sizeof(T) / sizeof(number_type<T>),
                      dtype, MPI_SUM, lattice->mpi_comm_lat);
        for (int i = 0; i < N; i++)
            value[i] = work[i];
    } else {
        MPI_Reduce((void *)value, (void *)work, N * sizeof(T) / sizeof(number_type<T>),
                   dtype, MPI_SUM, 0, lattice->mpi_comm_lat);
        if (hila::myrank() == 0)
            for (int i = 0; i < N; i++)
                value[i] = work[i];
    }
    reduction_timer.stop();
}

// Product reduction template - so far only for int, float, dbl

template <typename T>
void lattice_struct::reduce_node_product(T *value, int N, bool distribute) {
    T work[N];
    MPI_Datatype dtype;

    if (std::is_same<T, int>::value) {
        dtype = MPI_INT;
    } else if (std::is_same<T, float>::value) {
        dtype = MPI_FLOAT;
    } else if (std::is_same<T, double>::value) {
        dtype = MPI_DOUBLE;
    } else {
        static_assert(sizeof(T) > 0, "Unknown number_type in reduce_node_product");
    }

    reduction_timer.start();
    if (distribute) {
        MPI_Allreduce((void *)value, (void *)work, N, dtype, MPI_PROD,
                      lattice->mpi_comm_lat);
        for (int i = 0; i < N; i++)
            value[i] = work[i];
    } else {
        MPI_Reduce((void *)value, (void *)work, N, dtype, MPI_PROD, 0,
                   lattice->mpi_comm_lat);
        if (hila::myrank() == 0)
            for (int i = 0; i < N; i++)
                value[i] = work[i];
    }
    reduction_timer.stop();
}

#endif // USE_MPI

#endif // COMM_MPI_H
