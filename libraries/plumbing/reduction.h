#ifndef REDUCTION_H_
#define REDUCTION_H_

#include "hila.h"

//#define USE_MPI
#ifdef USE_MPI
#include <mpi.h>
#include "com_mpi.h"

template <typename T, class Enable = void>
class Reduction {};

/// Start MPI reduction operations.  Number_type<T> must be defined
/// for the type T.  Put this operation outside of class, so that
/// it is accessible from all Reductions

template <typename T>
void do_reduce_operation(MPI_Op operation, Reduction<T> &r) {

    // if for some reason reduction is going on unfinished, wait.
    r.wait();

    if (r.is_nonblocking()) r.set_comm_on();

    MPI_Datatype dtype;
    int size;

    dtype = get_MPI_number_type<T>(size);

    assert(dtype != MPI_BYTE && "Unknown number_type in reduction");

    void *ptr = r.get_ptr();
    MPI_Request *rp = r.get_request();

    reduction_timer.start();
    if (r.is_allreduce()) {
        if (r.is_nonblocking()) {
            MPI_Iallreduce(MPI_IN_PLACE, ptr, sizeof(T) / sizeof(hila::number_type<T>),
                           dtype, operation, lattice->mpi_comm_lat, rp);
        } else {
            MPI_Allreduce(MPI_IN_PLACE, ptr, sizeof(T) / sizeof(hila::number_type<T>),
                          dtype, operation, lattice->mpi_comm_lat);
        }
    } else {
        if (hila::myrank() == 0) {
            if (r.is_nonblocking()) {
                MPI_Ireduce(MPI_IN_PLACE, ptr, sizeof(T) / sizeof(hila::number_type<T>),
                            dtype, operation, 0, lattice->mpi_comm_lat, rp);
            } else {
                MPI_Reduce(MPI_IN_PLACE, ptr, sizeof(T) / sizeof(hila::number_type<T>),
                           dtype, operation, 0, lattice->mpi_comm_lat);
            }
        } else {
            if (r.is_nonblocking()) {
                MPI_Ireduce(ptr, ptr, sizeof(T) / sizeof(hila::number_type<T>), dtype,
                            operation, 0, lattice->mpi_comm_lat, rp);
            } else {
                MPI_Reduce(ptr, ptr, sizeof(T) / sizeof(hila::number_type<T>), dtype,
                           operation, 0, lattice->mpi_comm_lat);
            }
        }
    }
    reduction_timer.stop();
}

//////////////////////////////////////////////////////////////////////////////////
/// Special reduction class: declare a reduction variable which
/// can be used in site loops
/// Example:
///    Reduction<Complex<double>> rv = 0;
///    onsites(ALL) rv += fv[X];
///
///    ...  // possibly do something not involving rv
///
/// Reduction can be modified to on/off: allreduce(), nonblocking(), delayed()
/// Example:
///    rv.allreduce(false).delayed();
///
/// Delayed reduction can be used e.g. in expressions like
///
///   Reduction<double> rv = 0;
///   rv.allreduce(false).delayed();
///   foralldir(d) {
///       onsites(ALL) rv += a[X+d]*a[X];
///   }
///   rv.reduce();
///
/// This does only one reduction operation, not for every onsites() -loop.
/// Result is the same
///
/// Reduction variable can be used again
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

    bool delay_is_on = false;   // status of the delayed reduction
    bool is_delayed_sum = true; // sum/product

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
    bool is_allreduce() {
        return is_allreduce_;
    }

    /// nonblocking(bool) turns allreduce on or off.  By default on.
    Reduction &nonblocking(bool b = true) {
        is_nonblocking_ = b;
        return *this;
    }
    bool is_nonblocking() {
        return is_nonblocking_;
    }

    /// deferred(bool) turns deferred on or off.  By default turns on.
    Reduction &delayed(bool b = true) {
        is_delayed_ = b;
        return *this;
    }
    bool is_delayed() {
        return is_delayed_;
    }

    void set_comm_on() {
        comm_is_on = true;
    }

    MPI_Request *get_request() {
        return &request;
    }

    /// Return value of the reduction variable.  Wait for the comms if needed.
    T &value() {
        reduce();
        return (T) * this;
    }

    // operator T() const {
    //     wait();
    //     return *this;
    // }

    /// Cast to root type
    operator T &() {
        return *this;
    }

    /// Method set is the same as assignment, but without return value
    /// No need to complete comms if going on
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    void set(S &rhs) {
        if (comm_is_on) MPI_Cancel(&request);

        comm_is_on = false;
        (T &)*this = rhs;
    }

    /// Assignment is used only outside site loops - drop comms if on, no need to wait
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    T &operator=(const S &rhs) {
        set(rhs);
        return *this;
    }

    ///  Make compound ops return void, unconventionally:
    /// these are used within site loops but their value is not well-defined.
    /// Thus  a = (r += b); is disallowed
    /// Within site loops hilapp will change operations
    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_plus<T, S>>::value,
                               int> = 0>
    void operator+=(const S &rhs) {
        (T &)*this += rhs;
    }

    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_minus<T, S>>::value,
                               int> = 0>
    void operator-=(const S &rhs) {
        (T &)*this -= rhs;
    }

    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_mul<T, S>>::value,
                               int> = 0>
    void operator*=(const S &rhs) {
        (T &)*this *= rhs;
    }

    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_div<T, S>>::value,
                               int> = 0>
    void operator/=(const S &rhs) {
        (T &)*this /= rhs;
    }

    /// Start sum reduction -- works only if the type T addition == element-wise
    /// addition. This is true for all hila predefined data types
    void reduce_sum_node(const T &v) {

        wait(); // wait for possible ongoing

        // add the node values to reduction var
        if (hila::myrank() == 0 || is_delayed_)
            (T &)*this += v;
        else
            (T &)*this = v;

        if (is_delayed_) {
            if (delay_is_on && is_delayed_sum == false) {
                assert(0 && "Cannot mix sum and product reductions!");
            }
            delay_is_on = true;
            is_delayed_sum = true;
        } else {
            do_reduce_operation(MPI_SUM, *this);
        }
    }

    /// Product reduction -- currently works only for scalar data types.
    /// For Complex, Matrix and Vector data product is not element-wise.
    /// TODO: Array or std::array ?
    /// TODO: implement using custom MPI ops (if needed)
    void reduce_product_node(const T &v) {
        static_assert(std::is_arithmetic<T>::value,
                      "Type not implemented for product reduction");

        wait();

        if (hila::myrank() == 0 || is_delayed_)
            (T &)*this *= v;
        else
            (T &)*this = v;

        if (is_delayed_) {
            if (delay_is_on && is_delayed_sum == true) {
                assert(0 && "Cannot mix sum and product reductions!");
            }
            delay_is_on = true;
            is_delayed_sum = false;
        } else {
            do_reduce_operation(MPI_PROD, *this);
        }
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

    /// For delayed reduction, reduce starts or completes the reduction operation
    void reduce() {
        if (delay_is_on) {
            delay_is_on = false;

            if (is_delayed_sum)
                do_reduce_operation(MPI_SUM, *this);
            else
                do_reduce_operation(MPI_PROD, *this);
        }
    }

    /// get_ptr gives the pointer to the data held by the variable
    void *get_ptr() {
        return (void *)this;
    }
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
    bool is_delayed_ = false;

    bool delay_is_on = false;   // status of the delayed reduction
    bool is_delayed_sum = true; // sum/product

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
    Reduction(std::nullptr_t np) {
        val = 0;
        comm_is_on = false;
    }

    // explicit Reduction(bool all) : val(0), is_allreduce(all) {}
    // explicit Reduction(const T &v, bool all) : val(v), is_allreduce(all) {}

    /// Assign from assignable type too
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
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
    Reduction &allreduce(bool b = true) {
        is_allreduce_ = b;
        return *this;
    }
    bool is_allreduce() {
        return is_allreduce_;
    }

    /// nonblocking(bool) turns allreduce on or off.  By default on.
    Reduction &nonblocking(bool b = true) {
        is_nonblocking_ = b;
        return *this;
    }
    bool is_nonblocking() {
        return is_nonblocking_;
    }

    /// deferred(bool) turns deferred on or off.  By default turns on.
    Reduction &delayed(bool b = true) {
        is_delayed_ = b;
        return *this;
    }
    bool is_delayed() {
        return is_delayed_;
    }

    void set_comm_on() {
        comm_is_on = true;
    }

    MPI_Request *get_request() {
        return &request;
    }

    /// Return value of the reduction variable.  Wait for the comms if needed.
    T &value() {
        return val;
    }

    // operator T() const {
    //     wait();
    //     return *this;
    // }

    /// cast to T, with lvalue
    operator T &() {
        return val;
    }
    operator T() const {
        return val;
    }

    /// Cast to any compatible value
    // template <typename S, std::enable_if_t<hila::is_assignable<S &, T>::value, int> =
    // 0> operator S() {
    //    return (S)val;
    //}

    /// Method set is the same as assignment, but without return value
    /// No need to complete comms if going on
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    void set(S &rhs) {
        if (comm_is_on) MPI_Cancel(&request);

        comm_is_on = false;
        val = rhs;
    }

    /// Assignment is used only outside site loops - drop comms if on, no need to wait
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    T &operator=(const S &rhs) {
        set(rhs);
        return val;
    }

    /// Assignment from 0
#pragma hila loop_function
    T &operator=(const std::nullptr_t n) {
        val = 0;
        return val;
    }

    /// Make compound ops return void, unconventionally:
    /// these are used within site loops but their value is not well-defined.
    /// Thus  a = (r += b); is disallowed
    /// Within site loops hilapp will change operations
    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_plus<T, S>>::value,
                               int> = 0>
    void operator+=(const S &rhs) {
        val += rhs;
    }

    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_minus<T, S>>::value,
                               int> = 0>
    void operator-=(const S &rhs) {
        val -= rhs;
    }

    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_mul<T, S>>::value,
                               int> = 0>
    void operator*=(const S &rhs) {
        val *= rhs;
    }

    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_div<T, S>>::value,
                               int> = 0>
    void operator/=(const S &rhs) {
        val /= rhs;
    }

    /// Start sum reduction -- works only if the type T addition == element-wise
    /// addition. This is true for all hila predefined data types
    void reduce_sum_node(const T &v) {

        // if there is ongoing reduction wait for it
        wait();

        // add the node values to reduction var
        if (hila::myrank() == 0 || is_delayed_)
            val += v;
        else
            val = v;

        if (is_delayed_) {
            if (delay_is_on && is_delayed_sum == false) {
                assert(0 && "Cannot mix sum and product reductions!");
            }
            delay_is_on = true;
            is_delayed_sum = true;
        } else {
            do_reduce_operation(MPI_SUM, *this);
        }
    }

    /// Product reduction -- currently works only for scalar data types.
    /// For Complex, Matrix and Vector data product is not element-wise.
    /// TODO: Array or std::array ?
    /// TODO: implement using custom MPI ops (if needed)
    void reduce_product_node(const T &v) {
        static_assert(std::is_arithmetic<T>::value,
                      "Type not implemented for product reduction");

        wait();

        if (hila::myrank() == 0 || is_delayed_)
            val *= v;
        else
            val = v;

        if (is_delayed_) {
            if (delay_is_on && is_delayed_sum == true) {
                assert(0 && "Cannot mix sum and product reductions!");
            }
            delay_is_on = true;
            is_delayed_sum = false;
        } else {
            do_reduce_operation(MPI_PROD, *this);
        }
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

    /// For delayed reduction, reduce starts or completes the reduction operation
    void reduce() {
        if (delay_is_on) {
            delay_is_on = false;

            if (is_delayed_sum)
                do_reduce_operation(MPI_SUM, *this);
            else
                do_reduce_operation(MPI_PROD, *this);
        }
    }

    /// get_ptr gives the pointer to the data held by the variable
    void *get_ptr() {
        return (void *)&val;
    }
};

/// Standard output stream function for arithmetic class types
template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
std::ostream &operator<<(std::ostream &strm, Reduction<T> &r) {
    return strm << r.value();
}

#if defined(CUDA) || defined(HIP)
#include "backend_cuda/gpu_reduction.h"

template <typename T>
T Field<T>::reduce_sum(bool allreduce) const {
    return gpu_reduce_field(*this, allreduce);
}

#else

template <typename T>
T Field<T>::reduce_sum(bool allreduce) const {

    Reduction<T> result;
    result.allreduce(allreduce);
    onsites (ALL)
        result += (*this)[X];
    return result.value();
}

#endif

#endif // USE_MPI

// TODO - define for !USE_MPI

#endif
