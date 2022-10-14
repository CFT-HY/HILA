#ifndef COM_MPI_H
#define COM_MPI_H

#include "plumbing/defs.h"

#ifdef USE_MPI

#include "plumbing/lattice.h"

/// let us house the partitions-struct here
namespace hila {
class partitions_struct {
  public:
    unsigned _number, _mylattice;
    bool _sync;

    unsigned number() {
        return _number;
    }
    unsigned mylattice() {
        return _mylattice;
    }
    bool sync() {
        return _sync;
    }
};

extern partitions_struct partitions;

} // namespace hila

// Pile of timers associated with MPI calls
// clang-format off
extern hila::timer 
        start_send_timer, 
        wait_send_timer,
        post_receive_timer,
        wait_receive_timer,
        synchronize_timer,
        reduction_timer,
        reduction_wait_timer,
        broadcast_timer,
        send_timer,
        cancel_send_timer,
        cancel_receive_timer,
        partition_sync_timer;
// clang-format on

///***********************************************************
/// Implementations of communication routines.
///

// The MPI tag generator
int get_next_msg_tag();

// Obtain the MPI data type (MPI_XXX) for a particular type of native numbers.
// Boolean flag "with_int" is used to return "type + int" - types of MPI, used
// in maxloc/minloc reductions

template <typename T>
MPI_Datatype get_MPI_number_type(int &size, bool with_int = false) {

    if (std::is_same<hila::number_type<T>, int>::value) {
        size = sizeof(int);
        return with_int ? MPI_2INT : MPI_INT;
    } else if (std::is_same<hila::number_type<T>, unsigned>::value) {
        size = sizeof(unsigned);
        return with_int ? MPI_2INT : MPI_UNSIGNED; // MPI does not contain MPI_UNSIGNED_INT
    } else if (std::is_same<hila::number_type<T>, long>::value) {
        size = sizeof(long);
        return with_int ? MPI_LONG_INT : MPI_LONG;
    } else if (std::is_same<hila::number_type<T>, int64_t>::value) {
        size = sizeof(int64_t);
        return with_int ? MPI_LONG_INT : MPI_INT64_T; // need to use LONG_INT
    } else if (std::is_same<hila::number_type<T>, uint64_t>::value) {
        size = sizeof(uint64_t);
        return with_int ? MPI_LONG_INT : MPI_UINT64_T; // ditto
    } else if (std::is_same<hila::number_type<T>, float>::value) {
        size = sizeof(float);
        return with_int ? MPI_FLOAT_INT : MPI_FLOAT;
    } else if (std::is_same<hila::number_type<T>, double>::value) {
        size = sizeof(double);
        return with_int ? MPI_DOUBLE_INT : MPI_DOUBLE;
    } else if (std::is_same<hila::number_type<T>, long double>::value) {
        size = sizeof(long double);
        return with_int ? MPI_LONG_DOUBLE_INT : MPI_LONG_DOUBLE;
    }

    size = 1;
    return MPI_BYTE;
}


namespace hila {

///
/// Broadcast the value of _var_ to all nodes from node _rank_ (default=0).
/// Var must be trivial, i.e. plain data. 
/// Returns the broadcast value
/// If var is modifiable, it is changed to the broadcast value
///

template <typename T>
T broadcast(T &var, int rank = 0) {
    static_assert(std::is_trivial<T>::value, "broadcast(var) must use trivial type");
    if (hila::check_input)
        return var;

    assert(0 <= rank && rank < hila::number_of_nodes() && "Invalid sender rank in broadcast()");

    broadcast_timer.start();
    MPI_Bcast(&var, sizeof(T), MPI_BYTE, 0, lattice->mpi_comm_lat);
    broadcast_timer.stop();
    return var;
}

/// Version of broadcast with non-modifiable var
template <typename T>
T broadcast(const T &var, int rank = 0) {
    T tmp = var;
    return broadcast(tmp, rank);
}

/// Broadcast for std::vector
template <typename T>
void broadcast(std::vector<T> &list, int rank = 0) {

    static_assert(std::is_trivial<T>::value, "broadcast(vector<T>) must have trivial T");

    if (hila::check_input)
        return;

    broadcast_timer.start();

    int size = list.size();
    MPI_Bcast(&size, sizeof(int), MPI_BYTE, rank, lattice->mpi_comm_lat);
    if (hila::myrank() != rank) {
        list.resize(size);
    }

    // move vectors directly to the storage
    MPI_Bcast((void *)list.data(), sizeof(T) * size, MPI_BYTE, rank, lattice->mpi_comm_lat);

    broadcast_timer.stop();
}

///
/// Bare pointers cannot be broadcast

template <typename T>
void broadcast(T *var) {
    static_assert(sizeof(T) > 0 &&
                  "Do not use pointers to broadcast()-function. Use 'broadcast_array(T* arr, "
                  "int size)' to broadcast an array");
}

///
/// Broadcast for arrays where size must be known and same for all nodes

template <typename T>
void broadcast_array(T *var, int n, int rank = 0) {

    if (hila::check_input)
        return;

    broadcast_timer.start();
    MPI_Bcast((void *)var, sizeof(T) * n, MPI_BYTE, rank, lattice->mpi_comm_lat);
    broadcast_timer.stop();
}

// DO string bcasts separately
void broadcast(std::string &r, int rank = 0);
void broadcast(std::vector<std::string> &l, int rank = 0);

/// and broadcast with two values
template <typename T, typename U>
void broadcast2(T &t, U &u, int rank = 0) {

    if (hila::check_input)
        return;

    struct {
        T tv;
        U uv;
    } s = {t, u};

    hila::broadcast(s, rank);
    t = s.tv;
    u = s.uv;
}


///
/// Reduce an array across nodes

template <typename T>
void reduce_node_sum(T *value, int send_count, bool allreduce = true) {
    T recv_data[send_count];
    MPI_Datatype dtype;
    int size;

    if (hila::check_input)
        return;

    dtype = get_MPI_number_type<T>(size);

    reduction_timer.start();
    if (allreduce) {
        MPI_Allreduce((void *)value, (void *)recv_data,
                      send_count * sizeof(T) / sizeof(hila::number_type<T>), dtype, MPI_SUM,
                      lattice->mpi_comm_lat);
        for (int i = 0; i < send_count; i++)
            value[i] = recv_data[i];
    } else {
        MPI_Reduce((void *)value, (void *)recv_data,
                   send_count * sizeof(T) / sizeof(hila::number_type<T>), dtype, MPI_SUM, 0,
                   lattice->mpi_comm_lat);
        if (hila::myrank() == 0)
            for (int i = 0; i < send_count; i++)
                value[i] = recv_data[i];
    }
    reduction_timer.stop();
}

///
/// Reduce single variable across nodes.

template <typename T>
T reduce_node_sum(T &var, bool allreduce = true) {
    hila::reduce_node_sum(&var, 1, allreduce);
    return var;
}

// Product reduction template - so far only for int, float, dbl

template <typename T>
void reduce_node_product(T *send_data, int send_count, bool allreduce = true) {
    T recv_data[send_count];
    MPI_Datatype dtype;
    int size;

    if (hila::check_input)
        return;

    dtype = get_MPI_number_type<T>(size);

    reduction_timer.start();
    if (allreduce) {
        MPI_Allreduce((void *)send_data, (void *)recv_data, send_count, dtype, MPI_PROD,
                      lattice->mpi_comm_lat);
        for (int i = 0; i < send_count; i++)
            send_data[i] = recv_data[i];
    } else {
        MPI_Reduce((void *)send_data, (void *)recv_data, send_count, dtype, MPI_PROD, 0,
                   lattice->mpi_comm_lat);
        if (hila::myrank() == 0)
            for (int i = 0; i < send_count; i++)
                send_data[i] = recv_data[i];
    }
    reduction_timer.stop();
}

template <typename T>
T reduce_node_product(T &var, bool allreduce = true) {
    reduce_node_product(&var, 1, allreduce);
    return var;
}

void set_allreduce(bool on = true);
bool get_allreduce();


} // namespace hila


void hila_reduce_double_setup(double *d, int n);
void hila_reduce_float_setup(float *d, int n);
void hila_reduce_sums();


template <typename T>
void hila_reduce_sum_setup(T *value) {

    using b_t = hila::number_type<T>;
    if (std::is_same<b_t, double>::value) {
        hila_reduce_double_setup((double *)value, sizeof(T) / sizeof(double));
    } else if (std::is_same<b_t, float>::value) {
        hila_reduce_float_setup((float *)value, sizeof(T) / sizeof(float));
    } else {
        hila::reduce_node_sum(value, 1, hila::get_allreduce());
    }
}


#endif // USE_MPI

#endif // COMM_MPI_H
