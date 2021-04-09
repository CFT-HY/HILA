#ifndef FFT_FIELD_H
#define FFT_FIELD_H

#include "datatypes/cmplx.h"
#include "fftw3.h"
#include "plumbing/coordinates.h"
#include "plumbing/field.h"
#include "plumbing/timing.h"

/// Initialize fft to Direction dir.  _elements: number of complex values in
/// field, T_size the size of the field variable.  fftdir is
/// fft_direction::FORWARD or fft

void init_fft_direction(Direction dir, size_t _elements, size_t T_size,
                        fft_direction fftdir, void *const buffer1, void *const buffer2);

/// Build offsets to buffer arrays:
///   Fastest Direction = dir, offset 1
///   next fastest index of complex_t elements in T, offset elem_offset
///   and then other directions, in order
/// Returns element_offset and sets offset and nmin vectors

size_t fft_get_buffer_offsets(const Direction dir, const size_t elements,
                              CoordinateVector &offset, CoordinateVector &nmin);

/// FFT execute does the actual fft.  It should be called as
/// fft_execute<type>();  where type is the cmplx type
template <typename cmplx_t> void fft_execute() { assert("Don't call this!"); }

template <> void fft_execute<Complex<double>>();
template <> void fft_execute<Complex<float>>();

// Define type to map from T to complex_t
template <typename T, typename cmplx_t> union T_union {
    T val;
    cmplx_t c[sizeof(T) / sizeof(cmplx_t)];
};

/// Collect the data from field to buffer for sending or fft'ing.
/// Order: Direction dir goes fastest, then the index to complex data in T,
/// and other directions are slowest.

template <typename T, typename cmplx_t>
inline void fft_collect_data(const Field<T> &f, const Direction dir,
                             cmplx_t *const RESTRICT buffer) {

    constexpr int elements = sizeof(T) / sizeof(cmplx_t);

    extern hila::timer fft_collect_timer;
    fft_collect_timer.start();

    // Build vector offset, which encodes where the data should be written
    // elem_offset is the same for the offset of the elements of T
    CoordinateVector offset, nmin;

    const size_t elem_offset =
        fft_get_buffer_offsets(dir, sizeof(T) / sizeof(cmplx_t), offset, nmin);

// and collect the data
#pragma hila novector direct_access(buffer)
    onsites(ALL) {

        T_union<T, cmplx_t> v;
        v.val = f[X];
        size_t off = offset.dot(X.coordinates() - nmin);
        for (int i = 0; i < elements; i++) {

            buffer[off + i * elem_offset] = v.c[i];
        }
    }

    fft_collect_timer.stop();
}

/// Inverse of the fft_collect_data: writeh fft'd data to field.

template <typename T, typename cmplx_t>
inline void fft_save_result(Field<T> &f, const Direction dir,
                            const cmplx_t *const RESTRICT buffer) {

    constexpr size_t elements = sizeof(T) / sizeof(cmplx_t); // cmplx elements in T

    extern hila::timer fft_save_timer;
    fft_save_timer.start();

    // Build vector offset, which encodes where the data should be written
    CoordinateVector offset, nmin;

    const size_t elem_offset = fft_get_buffer_offsets(dir, elements, offset, nmin);

// and collect the data from buffers
#pragma hila novector direct_access(buffer)
    onsites(ALL) {

        T_union<T, cmplx_t> v;

        size_t off = offset.dot(X.coordinates() - nmin);
        for (size_t i = 0; i < elements; i++) {
            v.c[i] = buffer[off + i * elem_offset];
        }
        f[X] = v.val;
    }

    fft_save_timer.stop();
}

/// Increment the coordinates in order of odirs-array

inline void increment_current_coord(CoordinateVector &current,
                                    Direction odirs[NDIM - 1]) {

    int i = 0;
    while (i < NDIM - 1 && ++current[odirs[i]] >= lattice->mynode.size[odirs[i]]) {
        current[odirs[i]] = 0;
        ++i;
    }
}

///  Reshuffle data, given that previous fft dir was to prev_dir and now to dir
///  We order this so that the in buffer is accessed in order, out in ~random
///  order - theoretically random writes are faster than reads

template <typename cmplx_t>
inline void fft_reshuffle_data(const Direction fft_dir, cmplx_t *const RESTRICT out,
                               const Direction prev_dir,
                               const cmplx_t *const RESTRICT in,
                               const size_t elements) {

    extern hila::timer fft_reshuffle_timer;
    fft_reshuffle_timer.start();

    CoordinateVector offset_in, offset_out, nmin;

    const size_t e_offset_in =
        fft_get_buffer_offsets(prev_dir, elements, offset_in, nmin);
    const size_t e_offset_out =
        fft_get_buffer_offsets(fft_dir, elements, offset_out, nmin);

#if 1
    // Don't understand why this below seems to be faster than the more "direct"
    // construct below. This is random in - random out, whereas the below has
    // one of them in order

#pragma hila novector direct_access(out, in)
    onsites(ALL) {
        CoordinateVector v = X.coordinates() - nmin;
        size_t off_in = offset_in.dot(v);
        size_t off_out = offset_out.dot(v);
        for (size_t e = 0; e < elements; e++) {
            out[off_out + e * e_offset_out] = in[off_in + e * e_offset_in];
        }
    }

#else

    // what order directions?  Tally these up to odirs-array (other directions)
    Direction odirs[NDIM - 1];
    int i = 0;
    size_t n_columns = 1;
    foralldir(d) if (d != fft_dir) {
        odirs[i++] = d;
        n_columns *= lattice->mynode.size[d];
    }

    CoordinateVector current(0);
    const CoordinateVector nodesize = lattice->mynode.size;
    const size_t ns = nodesize[(int)fft_dir];

    for (size_t col = 0; col < n_columns; col++) {

        current[fft_dir] = 0;
        size_t off_out = offset_out.dot(current);

        for (size_t e = 0; e < elements; e++) {
            for (size_t d0 = 0; d0 < ns; d0++) {
                current[fft_dir] = d0;

                size_t off_in = offset_in.dot(current);
                /// and copy cmplx number at a time
                out[d0 + off_out + e * e_offset_out] = in[e * e_offset_in + off_in];
                // memcpy(out + (off_out + e*e_offset_out),  in + d0 +
                // e*e_offset_in , sizeof(cmplx_t));
            }
        }
        increment_current_coord(current, odirs);
    }

#endif
    //   increment_current_coord(current,odirs);

    fft_reshuffle_timer.stop();
}

void fft_post_gather();
void fft_start_gather(void *buffer);
void fft_wait_send();
void fft_wait_receive();
void fft_post_scatter(void *buffer);
void fft_start_scatter();
void fft_cleanup();

template <typename T>
inline void FFT_field(const Field<T> &input, Field<T> &result,
                      fft_direction fftdir = fft_direction::forward) {

    static_assert(hila::contains_complex<T>::value,
                  "FFT_field argument fields must contain complex type");

    // get the type of the complex number here
    using cmplx_t = Complex<hila::number_type<T>>;
    constexpr size_t elements = sizeof(T) / sizeof(cmplx_t);

    assert(lattice == input.fs->lattice && "Default lattice mismatch in fft");

    size_t local_volume = lattice->mynode.volume();

    size_t max_volume = 1;
    foralldir(d) max_volume *= lattice->nodes.max_size[d];

    extern hila::timer fft_timer;
    fft_timer.start();

    // Make sure the result is allocated and mark it changed
    result.check_alloc();
    result.mark_changed(ALL);

    cmplx_t *const RESTRICT collect_buffer =
        (cmplx_t *)d_malloc(local_volume * sizeof(T));
    cmplx_t *const RESTRICT receive_buffer =
        (cmplx_t *)d_malloc(local_volume * sizeof(T));

    bool first_dir = true;
    Direction prev_dir;

    foralldir(dir) {

        init_fft_direction(dir, elements, sizeof(T), fftdir, collect_buffer,
                           receive_buffer);

        fft_post_gather();

        if (first_dir) {
            fft_collect_data(input, dir, collect_buffer);
            // in_p = &result;
        } else {
            fft_reshuffle_data(dir, collect_buffer, prev_dir, receive_buffer, elements);
        }

        fft_start_gather(collect_buffer);

        fft_wait_receive(); // possibility to interleave later?
        fft_post_scatter(receive_buffer);
        fft_wait_send();

        fft_execute<cmplx_t>();

        fft_start_scatter();
        fft_wait_receive();
        fft_wait_send();

        fft_cleanup();

        // fft_save_result( result, dir, receive_buffer );

        prev_dir = dir;
        first_dir = false;
    }

    fft_save_result(result, prev_dir, receive_buffer);

    d_free(collect_buffer);
    d_free(receive_buffer);

    fft_timer.stop();
}

#endif
