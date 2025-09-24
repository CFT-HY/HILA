#ifndef HILA_FFT_H
#define HILA_FFT_H

/** @file fft.h */

#if !defined(CUDA) && !defined(HIP)
#define USE_FFTW
#endif

#include "plumbing/defs.h"
#include "datatypes/cmplx.h"
#include "plumbing/coordinates.h"
#include "plumbing/field.h"
#include "plumbing/timing.h"

#ifdef USE_FFTW
#include <fftw3.h>
#endif

// just some values here, make less than 100 just in case
#define WRK_GATHER_TAG 42
#define WRK_SCATTER_TAG 43


/**
 * @brief Convert momentum space CoordinateVector to wave number k, where -pi/2 < k_i <= pi/2
 * 
 * CoordinateVector is (periodically) modded to valid lattice coordinate,
 * and folded so that if n_i > lattice.size(i), n_i = n_i - lattice.size(i)
 * now  k_i = 2 pi n_i / lattice.size(i)
 */


#pragma hila novector
template <typename T>
inline Vector<NDIM, double> CoordinateVector_t<T>::convert_to_k() const {
    Vector<NDIM, double> k;
    foralldir(d) {
        int n = pmod((*this).e(d), lattice.size(d));
        if (n > lattice.size(d) / 2)
            n -= lattice.size(d);

        k[d] = n * 2.0 * M_PI / lattice.size(d);
    }
    return k;
}

inline Vector<NDIM, double> convert_to_k(const CoordinateVector &cv) {
    return cv.convert_to_k();
}

// hold static fft node data structures
struct pencil_struct {
    int node;               // node rank to send stuff for fft:ing
    unsigned size_to_dir;   // size of "node" to fft-dir
    unsigned column_offset; // first perp-plane column to be handled by "node"
    unsigned column_number; // and number of columns to be sent
    size_t recv_buf_size;   // size of my fft collect buffer (in units of sizeof(T)
                            // for stuff received from / returned to "node"
};

//
extern std::vector<pencil_struct> hila_pencil_comms[NDIM];

/// Build offsets to buffer arrays:
///   Fastest Direction = dir, offset 1
///   next fastest index of complex_t elements in T, offset elem_offset
///   and then other directions, in order
/// Returns element_offset and sets offset and nmin vectors

size_t pencil_get_buffer_offsets(const Direction dir, const size_t elements,
                                 CoordinateVector &offset, CoordinateVector &nmin);

/// Initialize fft direction - defined in fft.cpp
void init_pencil_direction(Direction d);

// Helper class to transform data
template <typename T, typename cmplx_t>
union T_union {
    T val;
    cmplx_t c[sizeof(T) / sizeof(cmplx_t)];
};

/// Class to hold fft relevant variables - note: fft.cpp holds static info, which is not
/// here

template <typename cmplx_t>
class hila_fft {
  public:
    Direction dir;
    int elements;
    fft_direction fftdir;
    size_t buf_size;
    size_t local_volume;

    bool only_reflect;

    cmplx_t *send_buf;
    cmplx_t *receive_buf;

    // data structures which point to to-be-copied buffers
    std::vector<cmplx_t *> rec_p;
    std::vector<int> rec_size;

    // initialize fft, allocate buffers
    hila_fft(int _elements, fft_direction _fftdir, bool _reflect = false) {
        extern size_t pencil_recv_buf_size[NDIM];

        elements = _elements;
        fftdir = _fftdir;
        only_reflect = _reflect;

        local_volume = lattice->mynode.volume;

        // init dirs here at one go
        foralldir(d) init_pencil_direction(d);

        buf_size = 1;
        foralldir(d) {
            if (pencil_recv_buf_size[d] > buf_size)
                buf_size = pencil_recv_buf_size[d];
        }
        if (buf_size < local_volume)
            buf_size = local_volume;

        // get fully aligned buffer space
        send_buf = (cmplx_t *)d_malloc(buf_size * sizeof(cmplx_t) * elements);
        receive_buf = (cmplx_t *)d_malloc(buf_size * sizeof(cmplx_t) * elements);
        //        if (buf_size > 0)
        //            fft_wrk_buf = (cmplx_t *)d_malloc(buf_size * sizeof(cmplx_t) *
        //            elements);
    }

    ~hila_fft() {
        d_free(send_buf);
        d_free(receive_buf);
        // if (buf_size > 0)
        //  d_free(fft_wrk_buf);
    }

    // make_plan does the fft plan, as appropriate
    void make_plan();
    // the actual transform is done here.  Custom for fftw and others
    void transform();

    // reflection using special call
    void reflect();

    /////////////////////////////////////////////////////////////////////////
    /// Initialize fft to Direction dir.

    void setup_direction(Direction _dir) {

        dir = _dir;

        // now in transform itself
        // make_fft_plan();

        rec_p.resize(hila_pencil_comms[dir].size());
        rec_size.resize(hila_pencil_comms[dir].size());

        cmplx_t *p = receive_buf;
        int i = 0;
        for (pencil_struct &fn : hila_pencil_comms[dir]) {

            if (fn.node != hila::myrank()) {

                // usually, out/in buffer is the same
                rec_p[i] = p;
                rec_size[i] = fn.size_to_dir;
                p += fn.recv_buf_size * elements;

            } else {

                // for local node, point directly to send_buf arrays

                rec_p[i] = send_buf + fn.column_offset * elements;
                rec_size[i] = fn.size_to_dir;
            }
            i++;
        }
    }

    /// Collect the data from field to send_buf for sending or fft'ing.
    /// Order: Direction dir goes fastest, then the index to complex data in T,
    /// and other directions are slowest.

    template <typename T>
    void collect_data(const Field<T> &f) {

        extern hila::timer pencil_collect_timer;
        pencil_collect_timer.start();

        constexpr int elements = sizeof(T) / sizeof(cmplx_t);

        // Build vector offset, which encodes where the data should be written
        // elem_offset is the same for the offset of the elements of T
        CoordinateVector offset, nmin;

        const size_t elem_offset =
            pencil_get_buffer_offsets(dir, sizeof(T) / sizeof(cmplx_t), offset, nmin);

        cmplx_t *sb = send_buf;

        // and collect the data
#pragma hila novector direct_access(sb)
        onsites(ALL) {

            T_union<T, cmplx_t> v;
            v.val = f[X];
            int off = offset.dot(X.coordinates() - nmin);
            for (int i = 0; i < elements; i++) {
                sb[off + i * elem_offset] = v.c[i];
            }
        }

        pencil_collect_timer.stop();
    }

    /// Inverse of the fft_collect_data: write fft'd data from receive_buf to field.

    template <typename T>
    void save_result(Field<T> &f) {

        constexpr int elements = sizeof(T) / sizeof(cmplx_t);

        extern hila::timer pencil_save_timer;
        pencil_save_timer.start();

        // Build vector offset, which encodes where the data should be written
        CoordinateVector offset, nmin;

        const size_t elem_offset = pencil_get_buffer_offsets(dir, elements, offset, nmin);

        cmplx_t *rb = receive_buf;

// and collect the data from buffers
#pragma hila novector direct_access(rb)
        onsites(ALL) {

            T_union<T, cmplx_t> v;

            size_t off = offset.dot(X.coordinates() - nmin);
            for (int i = 0; i < elements; i++) {
                v.c[i] = rb[off + i * elem_offset];
            }
            f[X] = v.val;
        }

        pencil_save_timer.stop();
    }

    /////////////////////////////////////////////////////////////////////////////
    ///  Reshuffle data, given that previous fft dir was to prev_dir and now to dir
    ///  Assuming here that the data is in receive_buf after fft and copy to send_buf
    ///  This requires swapping send_buf and receive_buf ptrs after 1 fft

    void reshuffle_data(Direction prev_dir) {

        extern hila::timer pencil_reshuffle_timer;
        pencil_reshuffle_timer.start();

        int elem = elements;

        CoordinateVector offset_in, offset_out, nmin;

        const size_t e_offset_in = pencil_get_buffer_offsets(prev_dir, elements, offset_in, nmin);
        const size_t e_offset_out = pencil_get_buffer_offsets(dir, elements, offset_out, nmin);

        cmplx_t *sb = send_buf;
        cmplx_t *rb = receive_buf;

#pragma hila novector direct_access(sb, rb)
        onsites(ALL) {
            CoordinateVector v = X.coordinates() - nmin;
            size_t off_in = offset_in.dot(v);
            size_t off_out = offset_out.dot(v);
            for (int e = 0; e < elem; e++) {
                sb[off_out + e * e_offset_out] = rb[off_in + e * e_offset_in];
            }
        }

        pencil_reshuffle_timer.stop();
    }

    // free the work buffers
    void cleanup() {}

    // just swap the buf pointers
    void swap_buffers() {
        std::swap(send_buf, receive_buf);
    }

    // communication functions for slicing the lattice
    void scatter_data();
    void gather_data();

    ////////////////////////////////////////////////////////////////////////
    /// Do the transform itself (fft or reflect only)

    template <typename T>
    void full_transform(const Field<T> &input, Field<T> &result,
                        const CoordinateVector &directions) {

        // Make sure the result is allocated and mark it changed
        result.check_alloc();

        bool first_dir = true;
        Direction prev_dir;

        foralldir(dir) {
            if (directions[dir]) {

                setup_direction(dir);

                if (first_dir) {
                    collect_data(input);
                    // in_p = &result;
                } else {
                    reshuffle_data(prev_dir);
                }

                gather_data();

                if (!only_reflect)
                    transform();
                else
                    reflect();

                scatter_data();

                cleanup();

                // fft_save_result( result, dir, receive_buf );

                prev_dir = dir;
                first_dir = false;

                // swap the pointers
                swap_buffers();
            }
        }

        save_result(result);

        result.mark_changed(ALL);
    }
};

// prototype for plan deletion
void FFT_delete_plans();


// Implementation dependent core fft collect and transforms are defined here

#if defined(USE_FFTW)

#include "plumbing/fft_fftw_transform.h"

#elif defined(CUDA) || defined(HIP)

#include "plumbing/backend_gpu/fft_gpu_transform.h"

#endif

/////////////////////////////////////////////////////////////////////////////////////////
/// Complex-to-complex FFT transform of a field input, result in result.
/// input and result can be same, "in-place".
///
/// Both input and output are of type Field<T>, where T must contain complex type, 
/// either Complex<float> or Complex<double>. 
/// As an example, if T is Vector<Complex<double>,3> the result has the same type
/// and the 3 components of the vector will contain FFTs of the input components.
/// 
/// directions: if directions[dir] == false (or 0), transform is not done to
/// coordinate direction dir. 
/// 
/// fftdir: direction of the transform itself:
///     fft_direction::forward (default)  x -> k
///     fft_direction::inverse  k-> x
///
/// FFT is unnormalized: transform + inverse transform yields source multiplied
///     by the product of the size of the lattice to active directions
///     If all directions are active, result = source * lattice.volume()
///
/// See also Field<T>::FFT()
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
inline void FFT_field(const Field<T> &input, Field<T> &result, const CoordinateVector &directions,
                      fft_direction fftdir = fft_direction::forward) {

    static_assert(hila::contains_complex<T>::value,
                  "FFT_field argument fields must contain complex type");

    // get the type of the complex number here
    using cmplx_t = Complex<hila::arithmetic_type<T>>;
    constexpr size_t elements = sizeof(T) / sizeof(cmplx_t);

    extern hila::timer fft_timer;
    fft_timer.start();

    hila_fft<cmplx_t> fft(elements, fftdir);

    fft.full_transform(input, result, directions);

    fft_timer.stop();
}

//////////////////////////////////////////////////////////////////////////////////
///
/// Complex-to-complex FFT transform of a field input, result in result.
/// Same as FFT_field(input,result,directions,fftdir)
/// with all directions active.
///
/// See also Field<T>::FFT() 
//////////////////////////////////////////////////////////////////////////////////

template <typename T>
inline void FFT_field(const Field<T> &input, Field<T> &result,
                      fft_direction fftdir = fft_direction::forward) {

    CoordinateVector dirs;
    dirs.fill(true); // set all directions OK

    FFT_field(input, result, dirs, fftdir);
}


/**
 * @brief Field method for performing FFT
 * @details
 * Used as 
 *   res = f.FFT(<args>);
 * where f and res are of type Field<T>, where T must contain complex type, either 
 * Complex<float> or Complex<double>. 
 * 
 * As an example, if T is Vector<Complex<double>,3> the result has the same type
 * and the 3 components of the vector will contain FFTs of the input components.
 * 
 * parameter dirs: if dirs[dir] == false (or 0), transform is not done to coordinate 
 * direction dir. By default calling without arguments will execute FFT in all directions.
 * 
 * parameter fftdir: direction of the transform itself:
 *      fft_direction::forward (default)  x -> k 
 *      fft_direction::inverse  k-> x
 * 
 * NOTE: the transform is unnormalized, i.e. forward + inverse transform yields
 * the original Field multiplied by the product of the size of the lattice to active directions.
 *
 * @code{.cpp}
 * .
 * . // Field f is defined
 * .
 * auto res = f.FFT() //Forward transform
 * auto res_2 = res.FFT(fft_direction::back) // res_2 == f * lattice.volume()
 * @endcode
 *
 * One can also specify the direction of the FFT with a coordinate vector:
 * @code{.cpp}
 * .
 * . // Field f is defined
 * .
 * auto res = f.FFT(e_x) //Forward transform in x-direction
 * auto res_2 = res.FFT(e_X,fft_direction::back) // res_2 == f * lattice.size(e_x)
 * @endcode
 *
 * With this in mind `f.FFT(e_x+e_y+e_z) = f.FFT()`
 *
 * @tparam T
 * @param dirs Direction to perform FFT in, default is all directions
 * @param fftdir fft_direction::forward (default)  or fft_direction::back
 * @return Field<T> Transformed field
 */
template <typename T>
Field<T> Field<T>::FFT(const CoordinateVector &dirs, fft_direction fftdir) const {
    Field<T> res;
    FFT_field(*this, res, dirs, fftdir);
    return res;
}

template <typename T>
Field<T> Field<T>::FFT(fft_direction fftdir) const {
    CoordinateVector cv;
    cv.fill(true);
    Field<T> res;
    FFT_field(*this, res, cv, fftdir);
    return res;
}


//////////////////////////////////////////////////////////////////////////////////
/// FFT_real_to_complex:
/// Field must be a real-valued field, result is a complex-valued field of the same type
/// Implemented just by doing a FFT with a complex field with im=0;
/// fft_direction::back gives a complex conjugate of the forward transform
/// Result is  f(-x) = f(L - x) = f(x)^*
//////////////////////////////////////////////////////////////////////////////////

template <typename T>
Field<Complex<hila::arithmetic_type<T>>> Field<T>::FFT_real_to_complex(fft_direction fftdir) const {

    static_assert(hila::is_arithmetic<T>::value,
                  "FFT_real_to_complex can be applied only to Field<real-type> variable");

    Field<Complex<T>> cf;
    cf[ALL] = Complex<T>((*this)[X], 0.0);
    return cf.FFT(fftdir);
}

//////////////////////////////////////////////////////////////////////////////////
/// @internal Helper function for complex to real FFT, used by FFT_complex_to_real()
//////////////////////////////////////////////////////////////////////////////////
namespace hila {
inline int FFT_complex_to_real_site(const CoordinateVector &cv) {

    // foralldir continues only if cv[d] == 0 or cv[d] == size(d)/2
    foralldir(d) {
        if (cv[d] > 0 && cv[d] < lattice.size(d) / 2)
            return 1;
        if (cv[d] > lattice.size(d) / 2)
            return -1;
    }
    // we get here only if all coords are 0 or size(d)/2
    return 0;
}

} // namespace hila



//////////////////////////////////////////////////////////////////////////////////
/// FFT_complex_to_real;
/// Field must be a complex-valued field, result is a real field of the same number type
/// Not optimized, should not be used on a hot path
///
/// Because the complex field must have the property f(-x) = f(L-x) = f(x)^*, only
/// half of the values in input field are significant, the routine does the appropriate
/// symmetrization.
///
/// Routine hila::FFT_complex_to_real_site(CoordinateVector cv) gives the significant values at
/// location cv:
///   = +1  significant complex value,
///   =  0  significant real part, imag ignored
///   = -1  value ignored here
/// Example: in 2d 8x8 lattice the sites are:  (* = (0,0), value 0)
///
///   - + + + - - - -                          - - - 0 + + + 0
///   - + + + - - - -                          - - - + + + + +
///   - + + + - - - -    after centering       - - - + + + + +
///   0 + + + 0 - - -    (0,0) to center       - - - + + + + +
///   + + + + + - - -    ----------------->    - - - * + + + 0
///   + + + + + - - -                          - - - - + + + -
///   + + + + + - - -                          - - - - + + + -
///   * + + + 0 - - -                          - - - - + + + -
///
//////////////////////////////////////////////////////////////////////////////////


template <typename T>
Field<hila::arithmetic_type<T>> Field<T>::FFT_complex_to_real(fft_direction fftdir) const {

    static_assert(hila::is_complex<T>::value,
                  "FFT_complex_to_real can be applied only to Field<Complex<>> type variable");

    foralldir(d) {
        assert(lattice.size(d) % 2 == 0 &&
               "FFT_complex_to_real works only with even lattice size to all directions");
    }

    // first, do a full reflection of the field, giving rf(x) = f(L-x) = "f(-x)"
    auto rf = this->reflect();
    // And symmetrize the field appropriately - can use rf
    onsites(ALL) {
        int type = hila::FFT_complex_to_real_site(X.coordinates());
        if (type == 1) {
            rf[X] = (*this)[X];
        } else if (type == -1) {
            rf[X] = rf[X].conj();
        } else {
            rf[X].real() = (*this)[X].real();
            rf[X].imag() = 0;
        }
    }

    FFT_field(rf, rf, fftdir);

    double ims = 0;
    double rss = 0;
    onsites(ALL) {
        ims += ::squarenorm(rf[X].imag());
        rss += ::squarenorm(rf[X].real());
    }

    Field<hila::arithmetic_type<T>> res;
    onsites(ALL) res[X] = rf[X].real();
    return res;
}


//////////////////////////////////////////////////////////////////////////////////
/// Field<T>::reflect() reflects the field around the desired axis
/// This is here because it uses similar communications as fft
/// TODO: refactorise so that there is separate "make columns" class!

/**
 * @brief  Reflect the Field around the desired axis
 * @details Can be called in the following ways:
 *
 * __Reflect on all axes:__
 *
 * \code {.cpp}
 * Field<MyType> f;
 * .
 * .
 * .
 * f.reflect()
 * \endcode
 *
 * @todo refactorise so that there is separate "make columns" class!
 * @tparam T
 * @param dirs
 * @return Field<T>
 */

template <typename T>
Field<T> Field<T>::reflect(const CoordinateVector &dirs) const {

    constexpr int elements = 1;

    Field<T> result;

    hila_fft<T> refl(elements, fft_direction::forward, true);
    refl.full_transform(*this, result, dirs);

    return result;
}

template <typename T>
Field<T> Field<T>::reflect() const {

    CoordinateVector c;
    c.fill(true);
    return reflect(c);
}

template <typename T>
Field<T> Field<T>::reflect(Direction dir) const {

    CoordinateVector c;
    c.fill(false);
    c[dir] = true;
    return reflect(c);
}


#endif
