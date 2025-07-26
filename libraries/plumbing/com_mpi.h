#ifndef COM_MPI_H
#define COM_MPI_H

#include "plumbing/defs.h"

#include "plumbing/lattice.h"

#include "datatypes/extended.h"


// let us house the partitions-struct here
namespace hila {
class partitions_struct {
  private:
    unsigned _number, _mylattice;
    bool _sync;

    public: 

    unsigned number() const {
        return _number;
    }

    void set_number(const unsigned u) {
        _number = u;
    }

    unsigned mylattice() const {
        return _mylattice;
    }

    void set_mylattice(const unsigned l) {
        _mylattice = l;
    }

    bool sync() const {
        return _sync;
    }

    void set_sync(bool s) {
        _sync = s;
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
        drop_comms_timer,
        partition_sync_timer;
// clang-format on

///***********************************************************
/// Implementations of communication routines.
///

// The MPI tag generator
int get_next_msg_tag();

/// Obtain the MPI data type (MPI_XXX) for a particular type of native numbers.
///
/// @brief Return MPI data type compatible with native number type
///
/// Boolean flag "with_int" is used to return "type + int" - types of MPI, used
/// in maxloc/minloc reductions
///
/// @tparam T  hila type, hila::arithmetic_type<T> is converted into MPI_type
/// @param size  (optional) size of the MPI_type in bytes
/// @param with_int
/// @return  MPI_datatype (e.g. MPI_INT, MPI_DOUBLE etc.)
template <typename T>
MPI_Datatype get_MPI_number_type(size_t &size, bool with_int = false) {

    if (std::is_same<hila::arithmetic_type<T>, int>::value) {
        size = sizeof(int);
        return with_int ? MPI_2INT : MPI_INT;
    } else if (std::is_same<hila::arithmetic_type<T>, unsigned>::value) {
        size = sizeof(unsigned);
        return with_int ? MPI_2INT : MPI_UNSIGNED; // MPI does not contain MPI_UNSIGNED_INT
    } else if (std::is_same<hila::arithmetic_type<T>, long>::value) {
        size = sizeof(long);
        return with_int ? MPI_LONG_INT : MPI_LONG;
    } else if (std::is_same<hila::arithmetic_type<T>, int64_t>::value) {
        size = sizeof(int64_t);
        return with_int ? MPI_LONG_INT : MPI_INT64_T; // need to use LONG_INT
    } else if (std::is_same<hila::arithmetic_type<T>, uint64_t>::value) {
        size = sizeof(uint64_t);
        return with_int ? MPI_LONG_INT : MPI_UINT64_T; // ditto
    } else if (std::is_same<hila::arithmetic_type<T>, float>::value) {
        size = sizeof(float);
        return with_int ? MPI_FLOAT_INT : MPI_FLOAT;
    } else if (std::is_same<hila::arithmetic_type<T>, double>::value) {
        size = sizeof(double);
        return with_int ? MPI_DOUBLE_INT : MPI_DOUBLE;
    } else if (std::is_same<hila::arithmetic_type<T>, long double>::value) {
        size = sizeof(long double);
        return with_int ? MPI_LONG_DOUBLE_INT : MPI_LONG_DOUBLE;
    }

    size = 1;
    return MPI_BYTE;
}


template <typename T>
MPI_Datatype get_MPI_number_type() {
    size_t s;
    return get_MPI_number_type<T>(s);
}


/// @brief Return MPI complex type equivalent to hila type
///
/// For example, if T is Complex<double>, get_MPI_complex_double<T>() returns MPI_C_DOUBLE_COMPLEX
/// @tparam T  type to be converted
/// @param siz  size of the type in bytes
/// @return MPI_Datatype
template <typename T>
MPI_Datatype get_MPI_complex_type(size_t &siz) {
    if constexpr (std::is_same<T, Complex<double>>::value) {
        siz = sizeof(Complex<double>);
        return MPI_C_DOUBLE_COMPLEX;
    } else if constexpr (std::is_same<T, Complex<float>>::value) {
        siz = sizeof(Complex<float>);
        return MPI_C_FLOAT_COMPLEX;
    } else {
        static_assert(sizeof(T) > 0,
                      "get_MPI_complex_type<T>() called without T being a complex type");
        return MPI_BYTE;
    }
}


namespace hila {

/**
 * @brief Broadcast the value of _var_ to all MPI ranks from _rank_ (default=0).
 *
 * NOTE: the function must be called by all MPI ranks, otherwise the program will deadlock.
 *
 * The type of the variable _var_ can be any standard plain datatype (trivial type),
 * std::string, std::vector or std::array
 *
 * For trivial types, the input _var_ can be non-modifiable value.  In this case
 * the broadcast value is obtained from the broadcast return value.
 *
 * Example:
 * @code{.cpp}
 *     auto rnd = hila::broadcast(hila::random());    // all MPI ranks get the same random value
 * @endcode
 *
 * @param var    variable to be synchronized across the full
 * @param rank   MPI rank from which the
 * @return template <typename T>
 */


template <typename T>
T broadcast(T &var, int rank = 0) {
    static_assert(std::is_trivial<T>::value, "broadcast(var) must use trivial type");
    if (hila::check_input)
        return var;

    assert(0 <= rank && rank < hila::number_of_nodes() && "Invalid sender rank in broadcast()");

    broadcast_timer.start();
    MPI_Bcast(&var, sizeof(T), MPI_BYTE, rank, lattice.mpi_comm_lat);
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

    static_assert(std::is_trivial<T>::value, "broadcast(std::vector<T>) must have trivial T");

    if (hila::check_input)
        return;

    broadcast_timer.start();

    int size = list.size();
    MPI_Bcast(&size, sizeof(int), MPI_BYTE, rank, lattice.mpi_comm_lat);
    if (hila::myrank() != rank) {
        list.resize(size);
    }

    // move vectors directly to the storage
    MPI_Bcast((void *)list.data(), sizeof(T) * size, MPI_BYTE, rank, lattice.mpi_comm_lat);

    broadcast_timer.stop();
}

/// And broadcast for std::array
template <typename T,int n>
void broadcast(std::array<T,n> &arr, int rank = 0) {

    static_assert(std::is_trivial<T>::value, "broadcast(std::array<T>) must have trivial T");

    if (hila::check_input)
        return;

    broadcast_timer.start();

    // move vectors directly to the storage
    MPI_Bcast((void *)arr.data(), sizeof(T) * n, MPI_BYTE, rank, lattice.mpi_comm_lat);

    broadcast_timer.stop();
}




///
/// Bare pointers cannot be broadcast

template <typename T>
void broadcast(T *var, int rank = 0) {
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
    MPI_Bcast((void *)var, sizeof(T) * n, MPI_BYTE, rank, lattice.mpi_comm_lat);
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


template <typename T>
void send_to(int to_rank, const T &data) {
    if (hila::check_input)
        return;

    send_timer.start();
    MPI_Send(&data, sizeof(T), MPI_BYTE, to_rank, hila::myrank(), lattice.mpi_comm_lat);
    send_timer.stop();
}

template <typename T>
void receive_from(int from_rank, T &data) {
    if (hila::check_input)
        return;

    send_timer.start();
    MPI_Recv(&data, sizeof(T), MPI_BYTE, from_rank, from_rank, lattice.mpi_comm_lat,
             MPI_STATUS_IGNORE);
    send_timer.stop();
}

template <typename T>
void send_to(int to_rank, const std::vector<T> &data) {
    if (hila::check_input)
        return;

    send_timer.start();
    size_t s = data.size();
    MPI_Send(&s, sizeof(size_t), MPI_BYTE, to_rank, hila::myrank(), lattice.mpi_comm_lat);

    MPI_Send(data.data(), sizeof(T) * s, MPI_BYTE, to_rank, hila::myrank(), lattice.mpi_comm_lat);
    send_timer.stop();
}

template <typename T>
void receive_from(int from_rank, std::vector<T> &data) {
    if (hila::check_input)
        return;

    send_timer.start();
    size_t s;
    MPI_Recv(&s, sizeof(size_t), MPI_BYTE, from_rank, from_rank, lattice.mpi_comm_lat,
             MPI_STATUS_IGNORE);
    data.resize(s);

    MPI_Recv(data.data(), sizeof(T) * s, MPI_BYTE, from_rank, from_rank, lattice.mpi_comm_lat,
             MPI_STATUS_IGNORE);
    send_timer.stop();
}


///
/// Reduce an array across nodes

template <typename T>
void reduce_node_sum(T *value, int send_count, bool allreduce = true) {

    if (hila::check_input)
        return;

    std::vector<T> recv_data(send_count);
    MPI_Datatype dtype;
    dtype = get_MPI_number_type<T>();
    reduction_timer.start();
    if (allreduce) {
        MPI_Allreduce((void *)value, (void *)recv_data.data(),
                      send_count * (sizeof(T) / sizeof(hila::arithmetic_type<T>)), dtype, MPI_SUM,
                      lattice.mpi_comm_lat);
        for (int i = 0; i < send_count; i++)
            value[i] = recv_data[i];
    } else {
        MPI_Reduce((void *)value, (void *)recv_data.data(),
                   send_count * (sizeof(T) / sizeof(hila::arithmetic_type<T>)), dtype, MPI_SUM, 0,
                   lattice.mpi_comm_lat);
        if (hila::myrank() == 0)
            for (int i = 0; i < send_count; i++)
                value[i] = recv_data[i];
    }
    reduction_timer.stop();
}

///
/// Reduce single variable across nodes.
/// A bit suboptimal, because uses std::vector

template <typename T>
T reduce_node_sum(T &var, bool allreduce = true) {
    hila::reduce_node_sum(&var, 1, allreduce);
    return var;
}

// Product reduction template - so far only for int, float, dbl

template <typename T>
void reduce_node_product(T *send_data, int send_count, bool allreduce = true) {
    std::vector<T> recv_data(send_count);
    MPI_Datatype dtype;

    if (hila::check_input)
        return;

    dtype = get_MPI_number_type<T>();

    reduction_timer.start();
    if (allreduce) {
        MPI_Allreduce((void *)send_data, (void *)recv_data.data(), send_count, dtype, MPI_PROD,
                      lattice.mpi_comm_lat);
        for (int i = 0; i < send_count; i++)
            send_data[i] = recv_data[i];
    } else {
        MPI_Reduce((void *)send_data, (void *)recv_data.data(), send_count, dtype, MPI_PROD, 0,
                   lattice.mpi_comm_lat);
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
void reduce_node_sum_extended(ExtendedPrecision *value, int send_count, bool allreduce = true);

extern MPI_Datatype MPI_ExtendedPrecision_type;
extern MPI_Op MPI_ExtendedPrecision_sum_op;

void create_extended_MPI_type();
void create_extended_MPI_operation();
void extended_sum_op(void *in, void *inout, int *len, MPI_Datatype *datatype);



template <typename T>
void hila_reduce_sum_setup(T *value) {

    using b_t = hila::arithmetic_type<T>;
    if constexpr (std::is_same<T, ExtendedPrecision>::value) {
        reduce_node_sum_extended(value, 1, hila::get_allreduce());
    } else if (std::is_same<b_t, double>::value) {
        hila_reduce_double_setup((double *)value, sizeof(T) / sizeof(double));
    } else if (std::is_same<b_t, float>::value) {
        hila_reduce_float_setup((float *)value, sizeof(T) / sizeof(float));
    } else {
        hila::reduce_node_sum(value, 1, hila::get_allreduce());
    }
}

#endif // COMM_MPI_H
