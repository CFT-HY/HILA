#ifndef FFT_FIELD_H
#define FFT_FIELD_H

#include "datatypes/cmplx.h"
#include "fftw3.h"
#include "plumbing/coordinates.h"
#include "plumbing/field.h"
#include "plumbing/timing.h"

/// Initialize fft to direction dir.  _elements: number of complex values in
/// field, T_size the size of the field variable.  fftdir is
/// fft_direction::FORWARD or fft

void init_fft_direction(direction dir, size_t _elements, size_t T_size,
                        fft_direction fftdir, void *const buffer1, void *const buffer2);

/// Build offsets to buffer arrays:
///   Fastest direction = dir, offset 1
///   next fastest index of complex_t elements in T, offset elem_offset
///   and then other directions, in order
/// Returns element_offset and sets offset and nmin vectors

size_t fft_get_buffer_offsets(const direction dir, const size_t elements,
                              CoordinateVector &offset, CoordinateVector &nmin);

/// FFT execute does the actual fft.  It should be called as
/// fft_execute<type>();  where type is the cmplx type
template <typename cmplx_t> void fft_execute() { assert("Don't call this!"); }

template <> void fft_execute<Cmplx<double>>();
template <> void fft_execute<Cmplx<float>>();

/// Collect the data from field to buffer for sending or fft'ing.
/// Order: Direction dir goes fastest, then the index to complex data in T,
/// and other directions are slowest.

template <typename T, typename cmplx_t>
inline void fft_collect_data(const Field<T> &f, const direction dir,
                             cmplx_t *const RESTRICT buffer) {

    constexpr size_t elements = sizeof(T) / sizeof(cmplx_t); // cmplx elements in T

    extern hila::timer fft_collect_timer;
    fft_collect_timer.start();

    // Use this union to convert from T to cmplx numbers
    union T_union {
        T val;
        cmplx_t c[elements];
    };

    // Build vector offset, which encodes where the data should be written
    // elem_offset is the same for the offset of the elements of T
    CoordinateVector offset, nmin;

    const size_t elem_offset = fft_get_buffer_offsets(dir, elements, offset, nmin);

// and collect the data
#pragma hila novector direct_access(buffer)
    onsites(ALL) {
        T_union v;
        v.val = f[X];
        size_t off = offset.dot(X.coordinates() - nmin);
        for (size_t i = 0; i < elements; i++) {

            buffer[off + i * elem_offset] = v.c[i];
        }
    }

    fft_collect_timer.stop();
}

/// Inverse of the fft_collect_data: writeh fft'd data to field.

template <typename T, typename cmplx_t>
inline void fft_save_result(Field<T> &f, const direction dir,
                            const cmplx_t *const RESTRICT buffer) {

    constexpr size_t elements = sizeof(T) / sizeof(cmplx_t); // cmplx elements in T

    extern hila::timer fft_save_timer;
    fft_save_timer.start();

    // Use this union to convert from T to cmplx numbers
    union T_union {
        T val;
        cmplx_t c[elements];
    };

    // Build vector offset, which encodes where the data should be written
    CoordinateVector offset, nmin;

    const size_t elem_offset = fft_get_buffer_offsets(dir, elements, offset, nmin);

// and collect the data from buffers
#pragma hila novector direct_access(buffer)
    onsites(ALL) {
        T_union v;

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
                                    direction odirs[NDIM - 1]) {

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
inline void fft_reshuffle_data(const direction fft_dir, cmplx_t *const RESTRICT out,
                               const direction prev_dir, const cmplx_t *const RESTRICT in,
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
    direction odirs[NDIM - 1];
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

template <typename T, typename cmplx_t>
inline void FFT_field_complex(const Field<T> &input, Field<T> &result,
                              fft_direction fftdir = fft_direction::forward) {

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
        (cmplx_t *)memalloc(local_volume * elements * sizeof(cmplx_t));
    cmplx_t *const RESTRICT receive_buffer =
        (cmplx_t *)memalloc(local_volume * elements * sizeof(cmplx_t));

    bool first_dir = true;
    direction prev_dir;

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

    free(collect_buffer);
    free(receive_buffer);

    fft_timer.stop();
}

/// Match a given type T to it's underlying complex type
template <typename T, class Enable = void> struct complex_base {};

/// Match to a complex type
template <> struct complex_base<Cmplx<float>> { using type = Cmplx<float>; };

/// Match to a complex type
template <> struct complex_base<Cmplx<double>> { using type = Cmplx<double>; };

/// Match templated class B to it's underlying complex type
template <template <typename B> class C, typename B> struct complex_base<C<B>> {
    using type = typename complex_base<B>::type;
};

/// Match templated class B to it's underlying complex type
template <template <int a, typename B> class C, int a, typename B>
struct complex_base<C<a, B>> {
    using type = typename complex_base<B>::type;
};

/// Match templated class B to it's underlying complex type
template <template <int a, int b, typename B> class C, int a, int b, typename B>
struct complex_base<C<a, b, B>> {
    using type = typename complex_base<B>::type;
};

/// Run fourier transform on a complex field
// Called with any type T with a Cmplx type nested in the lowest level
template <typename T>
void FFT_field(const Field<T> &input, Field<T> &result,
               fft_direction fdir = fft_direction::forward) {
    FFT_field_complex<T, typename complex_base<T>::type>(input, result, fdir);
}

#endif
