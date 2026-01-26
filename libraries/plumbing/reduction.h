#ifndef HILA_REDUCTION_H_
#define HILA_REDUCTION_H_

#include "hila.h"


//////////////////////////////////////////////////////////////////////////////////
/// @brief Special reduction class: enables delayed and non-blocking reductions, which
/// are not possible with the standard reduction. See [user guide](@ref reduction_type_guide)
/// for details
///

template <typename T>
class Reduction {

  private:
    // value holder
    T val;

    // comm_is_on is true if non-blocking MPI communications are under way.
    bool comm_is_on = false;

    // Reduction status : is this allreduce, nonblocking, or delayed
    bool is_allreduce_ = true;
    bool is_nonblocking_ = false;
    bool is_delayed_ = false;

    bool delay_is_on = false;   // status of the delayed reduction
    bool is_delayed_sum = true; // sum/product

    MPI_Request request;

    // start the actual reduction

    void do_reduce_operation(MPI_Op operation) {

        // if for some reason reduction is going on unfinished, wait.
        wait();

        if (is_nonblocking())
            comm_is_on = true;

        MPI_Datatype dtype;

        dtype = get_MPI_number_type<T>();

        assert(dtype != MPI_BYTE && "Unknown number_type in reduction");

        void *ptr = &val;

        reduction_timer.start();
        if (is_allreduce()) {
            if (is_nonblocking()) {
                MPI_Iallreduce(MPI_IN_PLACE, ptr, sizeof(T) / sizeof(hila::arithmetic_type<T>),
                               dtype, operation, lattice->mpi_comm_lat, &request);
            } else {
                MPI_Allreduce(MPI_IN_PLACE, ptr, sizeof(T) / sizeof(hila::arithmetic_type<T>),
                              dtype, operation, lattice->mpi_comm_lat);
            }
        } else {
            if (hila::myrank() == 0) {
                if (is_nonblocking()) {
                    MPI_Ireduce(MPI_IN_PLACE, ptr, sizeof(T) / sizeof(hila::arithmetic_type<T>),
                                dtype, operation, 0, lattice->mpi_comm_lat, &request);
                } else {
                    MPI_Reduce(MPI_IN_PLACE, ptr, sizeof(T) / sizeof(hila::arithmetic_type<T>),
                               dtype, operation, 0, lattice->mpi_comm_lat);
                }
            } else {
                if (is_nonblocking()) {
                    MPI_Ireduce(ptr, ptr, sizeof(T) / sizeof(hila::arithmetic_type<T>), dtype,
                                operation, 0, lattice->mpi_comm_lat, &request);
                } else {
                    MPI_Reduce(ptr, ptr, sizeof(T) / sizeof(hila::arithmetic_type<T>), dtype,
                               operation, 0, lattice->mpi_comm_lat);
                }
            }
        }
        reduction_timer.stop();
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

  public:
    /// Initialize to zero by default (? exception to other variables)
    /// allreduce = true by default
    Reduction() {
        val = 0;
        comm_is_on = false;
    }
    /// Initialize to value only on rank 0 - ensures expected result for delayed reduction
    ///
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    Reduction(const S &v) {
        if (hila::myrank() == 0) {
            val = v;
        } else {
            val = 0;
        }
        comm_is_on = false;
    }

    /// Delete copy construct from another reduction just in case
    /// don't make Reduction temporaries
    Reduction(const Reduction<T> &r) = delete;

    /// Destructor cleans up communications if they are in progress
    ~Reduction() {
        wait();
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

    /// Return value of the reduction variable.  Wait for the comms if needed.
    const T value() {
        reduce();
        return val;
    }


    /// Method set is the same as assignment, but without return value
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    void set(const S &rhs) {
        *this = rhs;
    }

    /// Assignment is used only outside site loops - drop comms if on, no need to wait
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    T operator=(const S &rhs) {
        wait();

        comm_is_on = false;
        T ret = rhs;
        if (hila::myrank() == 0) {
            val = ret;
        } else {
            val = 0;
        }
        // all ranks return the same value
        return ret;
    }

    /// Compound operator += is used in reduction but can be used outside onsites too.
    /// returns void, unconventionally:
    /// these are used within site loops but their value is not well-defined.
    /// Thus  a = (r += b); is disallowed
    template <typename S,
              std::enable_if_t<hila::is_assignable<T &, hila::type_plus<T, S>>::value, int> = 0>
    void operator+=(const S &rhs) {
        val += rhs;
    }

    // template <typename S,
    //           std::enable_if_t<hila::is_assignable<T &, hila::type_mul<T, S>>::value,
    //                            int> = 0>
    // void operator*=(const S &rhs) {
    //     val *= rhs;
    // }

    // template <typename S,
    //           std::enable_if_t<hila::is_assignable<T &, hila::type_div<T, S>>::value,
    //                            int> = 0>
    // void operator/=(const S &rhs) {
    //     val /= rhs;
    // }

    // Start sum reduction -- works only if the type T addition == element-wise
    // addition. This is true for all hila predefined data types
    void reduce_sum_node(const T &v) {

        wait(); // wait for possible ongoing

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
            do_reduce_operation(MPI_SUM);
        }
    }

    /// Product reduction -- currently works only for scalar data types.
    /// For Complex, Matrix and Vector data product is not element-wise.
    /// TODO: Array or std::array ?
    /// TODO: implement using custom MPI ops (if needed)
    // void reduce_product_node(const T &v) {

    //     reduce();

    //     if (hila::myrank() == 0 || is_delayed_)
    //         val *= v;
    //     else
    //         val = v;

    //     if (is_delayed_) {
    //         if (delay_is_on && is_delayed_sum == true) {
    //             assert(0 && "Cannot mix sum and product reductions!");
    //         }
    //         delay_is_on = true;
    //         is_delayed_sum = false;
    //     } else {
    //         do_reduce_operation(MPI_PROD);
    //     }
    // }


    /// For delayed reduction, start_reduce starts or completes the reduction operation
    void start_reduce() {
        if (!comm_is_on) {
            if (delay_is_on) {
                delay_is_on = false;

                if (is_delayed_sum)
                    do_reduce_operation(MPI_SUM);
                else
                    do_reduce_operation(MPI_PROD);
            }
        }
    }

    /// Complete the reduction - start if not done, and wait if ongoing
    void reduce() {
        start_reduce();
        wait();
    }
};


////////////////////////////////////////////////////////////////////////////////////

// #if defined(CUDA) || defined(HIP)
// #include "backend_gpu/gpu_reduction.h"

// template <typename T>
// T Field<T>::sum(Parity par, bool allreduce) const {
//     return gpu_reduce_sum(allreduce, par, false);
// }

// #else
// // This for not-gpu branch

template <typename T>
T Field<T>::sum(Parity par, bool allreduce) const {

    Reduction<T> result;
    result.allreduce(allreduce);
    onsites (par)
        result += (*this)[X];
    return result.value();
}

template <typename T>
T Field<T>::product(Parity par, bool allreduce) const {
    static_assert(std::is_arithmetic<T>::value,
                  ".product() reduction only for integer or floating point types");
    Reduction<T> result;
    result = 1;
    result.allreduce(allreduce);
    onsites (par)
        result *= (*this)[X];
    return result.value();
}


// get global minimum/maximums - meant to be used through .min() and .max()

#if defined(OPENMP) && !defined(HILAPP) 
#include <omp.h>
#endif

#if defined(CUDA) || defined(HIP)
// #include "backend_gpu/gpu_reduction.h"
#include "backend_gpu/gpu_minmax.h"
#endif

template <typename T>
T Field<T>::minmax(bool is_min, Parity par, CoordinateVector &loc) const {

    static_assert(std::is_same<T, int>::value || std::is_same<T, long>::value ||
                      std::is_same<T, float>::value || std::is_same<T, double>::value ||
                      std::is_same<T, long double>::value,
                  "In Field .min() and .max() methods the Field element type must be one of "
                  "(int/long/float/double/long double)");

#if defined(CUDA) || defined(HIP)
    T val = gpu_minmax(is_min, par, loc);
#else
    // get suitable initial value
    T val = is_min ? std::numeric_limits<T>::max() : std::numeric_limits<T>::min();

    // write the loop with explicit OpenMP parallel region.  It has negligible effect
    // on non-OpenMP code, and the pragmas are ignored.
    #pragma omp parallel shared(val, loc, is_min)
    {
        CoordinateVector loc_th(0);
        T val_th = is_min ? std::numeric_limits<T>::max() : std::numeric_limits<T>::min();

        // Pragma "hila omp_parallel_region" is necessary here, because this is within
        // omp parallel
        #pragma hila novector omp_parallel_region direct_access(loc_th, val_th)
        onsites (par) {
            if ((is_min && (*this)[X] < val_th) || (!is_min && (*this)[X] > val_th)) {
                val_th = (*this)[X];
                loc_th = X.coordinates();
            }
        }
        #pragma omp critical
        if ((is_min && val_th < val) || (!is_min && val_th > val)) {
            val = val_th;
            loc = loc_th;
        }
    }
    
#endif

    if (hila::number_of_nodes() > 1) {
        size_t size;
        // this returns MPI+int -type, for example MPI_DOUBLE_INT
        MPI_Datatype dtype = get_MPI_number_type<T>(size, true);

        struct {
            T v;
            int rank;
        } rdata;

        static_assert(sizeof(T) % sizeof(int) == 0,
                      "min/max reduction: datatype struct not packed!");

        rdata.v = val;
        rdata.rank = hila::myrank();

        // after allreduce rdata contains the min value and rank where it is
        if (is_min) {
            MPI_Allreduce(MPI_IN_PLACE, &rdata, 1, dtype, MPI_MINLOC, lattice->mpi_comm_lat);
        } else {
            MPI_Allreduce(MPI_IN_PLACE, &rdata, 1, dtype, MPI_MAXLOC, lattice->mpi_comm_lat);
        }
        val = rdata.v;

        // send the coordinatevector of the minloc to all nodes
        MPI_Bcast(&loc, sizeof(CoordinateVector), MPI_BYTE, rdata.rank, lattice->mpi_comm_lat);
    }

    return val;
}


/// Find minimum value from Field
template <typename T>
T Field<T>::min(Parity par) const {
    CoordinateVector loc;
    return minmax(true, par, loc);
}

/// Find minimum value and location from Field
template <typename T>
T Field<T>::min(CoordinateVector &loc) const {
    return minmax(true, ALL, loc);
}

/// Find minimum value and location from Field
template <typename T>
T Field<T>::min(Parity par, CoordinateVector &loc) const {
    return minmax(true, par, loc);
}


/// Find maximum value from Field
template <typename T>
T Field<T>::max(Parity par) const {
    CoordinateVector loc;
    return minmax(false, par, loc);
}

/// Find maximum value and location from Field
template <typename T>
T Field<T>::max(CoordinateVector &loc) const {
    return minmax(false, ALL, loc);
}

/// Find maximum value and location from Field
template <typename T>
T Field<T>::max(Parity par, CoordinateVector &loc) const {
    return minmax(false, par, loc);
}


#endif
