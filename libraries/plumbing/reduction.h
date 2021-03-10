#ifndef REDUCTION_H_
#define REDUCTION_H_

#include "hila.h"

//#define USE_MPI
#ifdef USE_MPI
#include <mpi.h>
#include "com_mpi.h"

template <typename T, class Enable = void> class Reduction {};

/// Start MPI reduction operations.  Number_type<T> must be defined
/// for the type T.  Put this operation outside of class, so that
/// it is accessible from all Reductions

template <typename T> void do_reduce_operation(MPI_Op operation, Reduction<T> &r) {

    // if for some reason reduction is going on unfinished, wait.
    r.wait();

    if (r.is_nonblocking())
        r.set_comm_on();

    MPI_Datatype dtype;
    int size;

    dtype = get_MPI_number_type<T>(size);

    if (dtype == MPI_BYTE) {
        assert(sizeof(T) < 0 && "Unknown number_type in reduce_node_sum");
    }

    void *ptr = r.get_ptr();
    MPI_Request *rp = r.get_request();

    reduction_timer.start();
    if (r.is_allreduce()) {
        if (r.is_nonblocking()) {
            MPI_Iallreduce(MPI_IN_PLACE, ptr, sizeof(T) / sizeof(number_type<T>), dtype,
                           operation, lattice->mpi_comm_lat, rp);
        } else {
            MPI_Allreduce(MPI_IN_PLACE, ptr, sizeof(T) / sizeof(number_type<T>), dtype,
                          operation, lattice->mpi_comm_lat);
        }
    } else {
        if (hila::myrank() == 0) {
            if (r.is_nonblocking()) {
                MPI_Ireduce(MPI_IN_PLACE, ptr, sizeof(T) / sizeof(number_type<T>),
                            dtype, operation, 0, lattice->mpi_comm_lat, rp);
            } else {
                MPI_Reduce(MPI_IN_PLACE, ptr, sizeof(T) / sizeof(number_type<T>), dtype,
                           operation, 0, lattice->mpi_comm_lat);
            }
        } else {
            if (r.is_nonblocking()) {
                MPI_Ireduce(ptr, ptr, sizeof(T) / sizeof(number_type<T>), dtype,
                            operation, 0, lattice->mpi_comm_lat, rp);
            } else {
                MPI_Reduce(ptr, ptr, sizeof(T) / sizeof(number_type<T>), dtype,
                           operation, 0, lattice->mpi_comm_lat);
            }
        }
    }
    reduction_timer.stop();
}

//////////////////////////////////////////////////////////////////////////////////
/// Special reduction class: declare a reduction variable which
/// can be used in site loops
/// Reductions use non-blocking reduction operations.  The reduction value is
/// calculated when the variable is used for the first time after summing it up
/// E.g.
///    Reduction<double> rv = 0;
///    onsites(ALL) rv += fv[X];
///
///    ...  // possibly do something not involving rv
///
///    // reduction is finalized when rv is first used
///    double average = rv/lattice->volume();
///
/// Reduction<T> r;        : value is broadcast to all nodes
/// Reduction<T> r(false); : value is reduced only to rank==0.
///                          value undefined on other nodes.
///
/// The same reduction variable can be used again
///

template <typename T>
class Reduction<T, typename std::enable_if<!std::is_arithmetic<T>::value>::type>
    : public T {

  private:
    // Value is derived from class T

    /// comm_is_on is true if non-blocking MPI communications are under way.
    bool comm_is_on = false;

    /// Reduction status : is this allreduce, nonblocking, or delayed
    bool is_allreduce_ = true;
    bool is_nonblocking_ = false;
    bool is_delayed_ = false;

    MPI_Request request;

  public:
    /// Initialize to zero by default (? exception to other variables)
    /// allreduce = true by default
    Reduction() {
        (T &)*this = 0;
        comm_is_on = false;
    }
    Reduction(const T &v) {
        (T &)*this = v;
        comm_is_on = false;
    }
    // explicit Reduction(bool all) : val(0), is_allreduce(all) {}
    // explicit Reduction(const T &v, bool all) : val(v), is_allreduce(all) {}


    /// Delete copy construct from another reduction just in case
    /// don't make Reduction temporaries
    Reduction(const Reduction<T> &r) = delete;

    /// Destructor cleans up communications if they are in progress
    ~Reduction() {
        if (comm_is_on) {
            MPI_Cancel(&request);
        }
    }

    /// allreduce(bool) turns allreduce on or off.  By default on.
    Reduction &allreduce(bool b = true) {
        is_allreduce_ = b;
        return *this;
    }
    bool is_allreduce() { return is_allreduce_; }

    /// nonblocking(bool) turns allreduce on or off.  By default on.
    Reduction &nonblocking(bool b = true) {
        is_nonblocking_ = b;
        return *this;
    }
    bool is_nonblocking() { return is_nonblocking_; }

    /// deferred(bool) turns deferred on or off.  By default turns on.
    Reduction &delayed(bool b = true) {
        is_delayed_ = b;
        return *this;
    }
    bool is_delayed() { return is_delayed_; }

    void set_comm_on() { comm_is_on = true; }
    MPI_Request *get_request() { return &request; }

    /// Return value of the reduction variable.  Wait for the comms if needed.
    T &value() { return (T) * this; }

    // operator T() const {
    //     wait();
    //     return *this;
    // }

    /// Cast to root type
    operator T &() { return *this; }

    /// Method set is the same as assignment, but without return value
    /// No need to complete comms if going on
    template <typename S, std::enable_if_t<std::is_assignable<T &, S>::value, int> = 0>
    void set(S &rhs) {
        if (comm_is_on)
            MPI_Cancel(&request);

        comm_is_on = false;
        (T &)*this = rhs;
    }

    /// Assignment is used only outside site loops - drop comms if on, no need to wait
    template <typename S, std::enable_if_t<std::is_assignable<T &, S>::value, int> = 0>
    T &operator=(const S &rhs) {
        set(rhs);
        return *this;
    }

    ///  Make compound ops return void, unconventionally:
    /// these are used within site loops but their value is not well-defined.
    /// Thus  a = (r += b); is disallowed
    /// Within site loops hilapp will change operations
    template <typename S, std::enable_if_t<
                              std::is_assignable<T &, type_plus<T, S>>::value, int> = 0>
    void operator+=(const S &rhs) {
        (T &)*this += rhs;
    }

    template <
        typename S,
        std::enable_if_t<std::is_assignable<T &, type_minus<T, S>>::value, int> = 0>
    void operator-=(const S &rhs) {
        (T &)*this -= rhs;
    }

    template <typename S,
              std::enable_if_t<std::is_assignable<T &, type_mul<T, S>>::value, int> = 0>
    void operator*=(const S &rhs) {
        (T &)*this *= rhs;
    }

    template <typename S,
              std::enable_if_t<std::is_assignable<T &, type_div<T, S>>::value, int> = 0>
    void operator/=(const S &rhs) {
        (T &)*this /= rhs;
    }

    /// Start sum reduction -- works only if the type T addition == element-wise
    /// addition. This is true for all hila predefined data types
    void reduce_sum() { do_reduce_operation(MPI_SUM, *this); }

    /// Product reduction -- currently works only for scalar data types.
    /// For Complex, Matrix and Vector data product is not element-wise.
    /// TODO: Array or std::array ?
    /// TODO: implement using custom MPI ops (if needed)
    void reduce_product() {
        static_assert(std::is_arithmetic<T>::value,
                      "Type not implemented for product reduction");
        do_reduce_operation(MPI_PROD, *this);
    }

    /// Wait for MPI to complete, if it is currently going on
    /// This must be called for non-blocking reduce before use!

    void wait() {
        if (comm_is_on) {
            reduction_wait_timer.start();
            MPI_Status status;
            MPI_Wait(&request, &status);
            reduction_wait_timer.stop();
            comm_is_on = false;
        }
    }

    /// get_ptr gives the pointer to the data held by the variable
    void *get_ptr() { return (void *)this; }
};

/////////////////////////////////////////////////////////////////////////////////
/// Specialization of Reduction to an arithmetic (non-class) type

template <typename T>
class Reduction<T, typename std::enable_if<std::is_arithmetic<T>::value>::type> {

  private:
    T val;

    /// comm_is_on is true if MPI communications are under way.
    bool comm_is_on = false;

    /// Reduction status : is this allreduce, nonblocking
    bool is_allreduce_ = true;
    bool is_nonblocking_ = false;
    MPI_Request request;

  public:
    /// Initialize to zero by default (? exception to other variables)
    /// allreduce = true by default
    Reduction() {
        val = 0;
        comm_is_on = false;
    }
    Reduction(const T &v) {
        val = v;
        comm_is_on = false;
    }
    // explicit Reduction(bool all) : val(0), is_allreduce(all) {}
    // explicit Reduction(const T &v, bool all) : val(v), is_allreduce(all) {}

    /// Assign from assignable type too
    template <typename S, std::enable_if_t<std::is_assignable<T &, S>::value, int> = 0>
    Reduction(const S &s) {
        val = s;
    }

    /// Delete copy construct from another reduction just in case
    /// don't make Reduction temporaries
    Reduction(const Reduction<T> &r) = delete;

    /// Destructor cleans up communications if they are in progress
    ~Reduction() {
        if (comm_is_on) {
            MPI_Cancel(&request);
        }
    }

    /// allreduce(bool) turns allreduce on or off.  By default on.
    void allreduce(bool b = true) { is_allreduce_ = b; }
    bool is_allreduce() { return is_allreduce_; }

    /// nonblocking(bool) turns allreduce on or off.  By default on.
    void nonblocking(bool b = true) { is_nonblocking_ = b; }
    bool is_nonblocking(bool b = true) { return is_nonblocking_; }

    void set_comm_on() { comm_is_on = true; }

    MPI_Request *get_request() { return &request; }

    /// Return value of the reduction variable.  Wait for the comms if needed.
    T &value() { return val; }

    // operator T() const {
    //     wait();
    //     return *this;
    // }

    /// cast to T, with lvalue
    operator T &() { return val; }

    /// Cast to any compatible value
    template <typename S, std::enable_if_t<std::is_assignable<S &, T>::value, int> = 0>
    operator S() {
        return (S)val;
    }

    /// Method set is the same as assignment, but without return value
    /// No need to complete comms if going on
    template <typename S, std::enable_if_t<std::is_assignable<T &, S>::value, int> = 0>
    void set(S &rhs) {
        if (comm_is_on)
            MPI_Cancel(&request);

        comm_is_on = false;
        val = rhs;
    }

    /// Assignment is used only outside site loops - drop comms if on, no need to wait
    template <typename S, std::enable_if_t<std::is_assignable<T &, S>::value, int> = 0>
    T &operator=(const S &rhs) {
        set(rhs);
        return val;
    }

    /// Make compound ops return void, unconventionally:
    /// these are used within site loops but their value is not well-defined.
    /// Thus  a = (r += b); is disallowed
    /// Within site loops hilapp will change operations
    template <typename S, std::enable_if_t<
                              std::is_assignable<T &, type_plus<T, S>>::value, int> = 0>
    void operator+=(const S &rhs) {
        val += rhs;
    }

    template <
        typename S,
        std::enable_if_t<std::is_assignable<T &, type_minus<T, S>>::value, int> = 0>
    void operator-=(const S &rhs) {
        val -= rhs;
    }

    template <typename S,
              std::enable_if_t<std::is_assignable<T &, type_mul<T, S>>::value, int> = 0>
    void operator*=(const S &rhs) {
        val *= rhs;
    }

    template <typename S,
              std::enable_if_t<std::is_assignable<T &, type_div<T, S>>::value, int> = 0>
    void operator/=(const S &rhs) {
        val /= rhs;
    }

    /// Start sum reduction -- works only if the type T addition == element-wise
    /// addition. This is true for all hila predefined data types
    void reduce_sum() { do_reduce_operation(MPI_SUM, *this); }

    /// Product reduction -- currently works only for scalar data types.
    /// For Complex, Matrix and Vector data product is not element-wise.
    /// TODO: Array or std::array ?
    /// TODO: implement using custom MPI ops (if needed)
    void reduce_product() { do_reduce_operation(MPI_PROD, *this); }

    /// Wait for MPI to complete, if it is currently going on
    /// This must be called for non-blocking reduce before use!

    void wait() {
        if (comm_is_on) {
            reduction_wait_timer.start();
            MPI_Status status;
            MPI_Wait(&request, &status);
            reduction_wait_timer.stop();
            comm_is_on = false;
        }
    }

    /// get_ptr gives the pointer to the data held by the variable
    void *get_ptr() { return (void *)&val; }
};

/// Standard output stream function for arithmetic class types
template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
std::ostream &operator<<(std::ostream &strm, Reduction<T> &r) {
    return strm << r.value();
}

#endif // USE_MPI

// TODO - define for !USE_MPI

#if 0


/// Standard arithmetic operations for reductions
template <typename T, typename S>
inline auto operator+(Reduction<T> & a, const S & b) -> type_plus<T, S> {
    return  a.value() + b;
}

template <typename T, typename S>
inline auto operator+(const S & b, Reduction<T> & a) -> type_plus<S, T> {
    return b + a.value();
}

template <typename T, typename S>
inline auto operator-(Reduction<T> & a, const S & b) -> type_minus<T, S> {
    return  a.value() - b;
}

template <typename T, typename S>
inline auto operator-(const S & b, Reduction<T> & a) -> type_minus<S, T> {
    return b - a.value();
}

template <typename T, typename S>
inline auto operator*(Reduction<T> & a, const S & b) -> type_mul<T, S> {
    return  a.value() * b;
}

template <typename T, typename S>
inline auto operator*(const S & b, Reduction<T> & a) -> type_mul<S, T> {
    return b * a.value();
}

template <typename T, typename S>
inline auto operator/(Reduction<T> & a, const S & b) -> type_div<T, S> {
    return  a.value() / b;
}

template <typename T, typename S>
inline auto operator/(const S & b, Reduction<T> & a) -> type_div<S, T> {
    return b / a.value();
}

#endif

#endif
