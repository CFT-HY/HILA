#ifndef HILA_FFT_H
#define HILA_FFT_H

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

// just some values here
#define WRK_GATHER_TAG 42
#define WRK_SCATTER_TAG 43

// hold static fft node data structures
struct fftnode_struct {
    int node;               // node rank to send stuff for fft:ing
    unsigned size_to_dir;   // size of "node" to fft-dir
    unsigned column_offset; // first perp-plane column to be handled by "node"
    unsigned column_number; // and number of columns to be sent
    size_t recv_buf_size;   // size of my fft collect buffer (in units of
                          // elements*cmplx_size) for stuff received from / returned to
                          // "node"
};

//
extern std::vector<fftnode_struct> hila_fft_comms[NDIM];

/// Build offsets to buffer arrays:
///   Fastest Direction = dir, offset 1
///   next fastest index of complex_t elements in T, offset elem_offset
///   and then other directions, in order
/// Returns element_offset and sets offset and nmin vectors

size_t fft_get_buffer_offsets(const Direction dir, const size_t elements,
                              CoordinateVector &offset, CoordinateVector &nmin);

/// Initialize fft direction - defined in fft.cpp
void init_fft_direction(Direction d);

// Helper class to transform data
template <typename T, typename cmplx_t> union T_union {
    T val;
    cmplx_t c[sizeof(T) / sizeof(cmplx_t)];
};

/// Class to hold fft relevant variables - note: fft.cpp holds static info, which is not
/// here

template <typename cmplx_t> class hila_fft {
  public:
    Direction dir;
    int elements;
    int cmplx_size;
    fft_direction fftdir;
    size_t buf_size;
    size_t local_volume;

    cmplx_t *RESTRICT send_buf;
    cmplx_t *RESTRICT receive_buf;

    // data structures which point to to-be-copied buffers
    std::vector<cmplx_t *> rec_p;
    std::vector<int> rec_size;

#ifdef USE_MPI
    MPI_Datatype mpi_cmplx_t;
#endif

    // initialize fft, allocate buffers
    hila_fft(int _elements, fft_direction _fftdir) {
        extern size_t fft_recv_buf_size[NDIM];

        elements = _elements;
        fftdir = _fftdir;
        cmplx_size = sizeof(cmplx_t);

        local_volume = lattice->mynode.volume();

        // init dirs here at one go
        foralldir(d) init_fft_direction(d);

        buf_size = 1;
        foralldir(d) {
            if (fft_recv_buf_size[d] > buf_size)
                buf_size = fft_recv_buf_size[d];
        }
        if (buf_size < local_volume)
            buf_size = local_volume;

        // get fully aligned buffer space
        send_buf = (cmplx_t *)d_malloc(buf_size * sizeof(cmplx_t) * elements);
        receive_buf = (cmplx_t *)d_malloc(buf_size * sizeof(cmplx_t) * elements);
        //        if (buf_size > 0)
        //            fft_wrk_buf = (cmplx_t *)d_malloc(buf_size * sizeof(cmplx_t) *
        //            elements);

#ifdef USE_MPI
        mpi_cmplx_t = (sizeof(cmplx_t) == sizeof(Complex<double>))
                          ? MPI_C_DOUBLE_COMPLEX
                          : MPI_C_FLOAT_COMPLEX;

#endif
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

    /////////////////////////////////////////////////////////////////////////
    /// Initialize fft to Direction dir.

    void setup_direction(Direction _dir) {

        dir = _dir;

        // now in transform itself
        // make_fft_plan();

        rec_p.resize(hila_fft_comms[dir].size());
        rec_size.resize(hila_fft_comms[dir].size());

        cmplx_t *p = receive_buf;
        int i = 0;
        for (fftnode_struct &fn : hila_fft_comms[dir]) {

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

    template <typename T> void collect_data(const Field<T> &f) {

        extern hila::timer fft_collect_timer;
        fft_collect_timer.start();

        constexpr int elements = sizeof(T) / sizeof(cmplx_t);

        // Build vector offset, which encodes where the data should be written
        // elem_offset is the same for the offset of the elements of T
        CoordinateVector offset, nmin;

        const size_t elem_offset =
            fft_get_buffer_offsets(dir, sizeof(T) / sizeof(cmplx_t), offset, nmin);

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

        fft_collect_timer.stop();
    }

    /// Inverse of the fft_collect_data: write fft'd data from receive_buf to field.

    template <typename T> void save_result(Field<T> &f) {

        constexpr int elements = sizeof(T) / sizeof(cmplx_t);

        extern hila::timer fft_save_timer;
        fft_save_timer.start();

        // Build vector offset, which encodes where the data should be written
        CoordinateVector offset, nmin;

        const size_t elem_offset = fft_get_buffer_offsets(dir, elements, offset, nmin);

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

        fft_save_timer.stop();
    }

    /////////////////////////////////////////////////////////////////////////////
    ///  Reshuffle data, given that previous fft dir was to prev_dir and now to dir
    ///  Assuming here that the data is in receive_buf after fft and copy to send_buf
    ///  This requires swapping send_buf and receive_buf ptrs after 1 fft

    void reshuffle_data(Direction prev_dir) {

        extern hila::timer fft_reshuffle_timer;
        fft_reshuffle_timer.start();

        int elem = elements;

        CoordinateVector offset_in, offset_out, nmin;

        const size_t e_offset_in =
            fft_get_buffer_offsets(prev_dir, elements, offset_in, nmin);
        const size_t e_offset_out =
            fft_get_buffer_offsets(dir, elements, offset_out, nmin);

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

        fft_reshuffle_timer.stop();
    }

    // free the work buffers
    void cleanup() {}

    // just swap the buf pointers
    void swap_buffers() { std::swap(send_buf, receive_buf); }

    // communication functions for slicing the lattice
    void scatter_data();
    void gather_data();
};

// Implementation dependent core fft collect and transforms are defined here

#if defined(USE_FFTW)

#include "plumbing/fft_fftw_transform.h"

#elif defined(CUDA) || defined(HIP)

#include "plumbing/backend_cuda/fft_hip_transform.h"

#endif

/////////////////////////////////////////////////////////////////////////////////////////
/// Complex-to-complex FFT transform of a field input, result in result.
/// input and result can be same, "in-place".
/// Both input and output are of type Field<T>, where T must contain complex type,
/// Complex<float> or Complex<double>.
/// directions: if directions[dir] == false (or 0), transform is not done to direction
/// dir. 
/// fftdir: direction of the transform itself:
///     fft_direction::forward (default)  x -> k
///     fft_direction::inverse  k-> x
/// FFT is unnormalized: transform + inverse transform yields source multiplied
///     by the product of the size of the lattice to active directions
///     If all directions are active, result = source * lattice->volume():
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
inline void FFT_field(const Field<T> &input, Field<T> &result,
                      const CoordinateVector &directions,
                      fft_direction fftdir = fft_direction::forward) {

    static_assert(hila::contains_complex<T>::value,
                  "FFT_field argument fields must contain complex type");

    // get the type of the complex number here
    using cmplx_t = Complex<hila::number_type<T>>;
    constexpr size_t elements = sizeof(T) / sizeof(cmplx_t);

    assert(lattice == input.fs->lattice && "Default lattice mismatch in fft");

    extern hila::timer fft_timer;
    fft_timer.start();

    // Make sure the result is allocated and mark it changed
    result.check_alloc();

    hila_fft<cmplx_t> fft(elements, fftdir);

    bool first_dir = true;
    Direction prev_dir;

    foralldir(dir) if (directions[dir]) {

        fft.setup_direction(dir);

        if (first_dir) {
            fft.collect_data(input);
            // in_p = &result;
        } else {
            fft.reshuffle_data(prev_dir);
        }

        fft.gather_data();

        fft.transform();

        fft.scatter_data();

        fft.cleanup();

        // fft_save_result( result, dir, receive_buf );

        prev_dir = dir;
        first_dir = false;

        // swap the pointers
        fft.swap_buffers();
    }

    fft.save_result(result);

    fft_timer.stop();

    result.mark_changed(ALL);
}

//////////////////////////////////////////////////////////////////////////////////
///
/// Complex-to-complex FFT transform of a field input, result in result.
/// Same as FFT_field(input,result,directions,fftdir)
/// with all directions active.
///
//////////////////////////////////////////////////////////////////////////////////

template <typename T>
inline void FFT_field(const Field<T> &input, Field<T> &result,
                      fft_direction fftdir = fft_direction::forward) {

    CoordinateVector dirs;
    dirs.asArray() = true; // set all directions OK

    FFT_field(input, result, dirs, fftdir);
}

//////////////////////////////////////////////////////////////////////////////////
///
/// Field method for performing FFT
///   a.FFT();   does in-place transform of field a
/// fftdir:  fft_direction::forward (default)  or fft_direction::inverse
//////////////////////////////////////////////////////////////////////////////////

template <typename T> void Field<T>::FFT(fft_direction fftdir) {
    FFT_field(*this, *this, fftdir);
}



#endif
