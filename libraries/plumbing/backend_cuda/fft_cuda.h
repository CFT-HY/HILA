#ifndef FFT_CUDA_H
#define FFT_CUDA_H

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


template <typename complex_type>
__global__ void gather_column(cufftDoubleComplex *data, char *field_elem, int elements,
                              int cpn, int nnodes, int node_column_size, int column_size,
                              int block_size) {
    int ind = threadIdx.x + blockIdx.x * blockDim.x;
    if (ind < node_column_size * cpn) {
        int t = ind % node_column_size;
        int l = ind / node_column_size;
        cufftDoubleComplex *d = data + l * column_size * elements;
        complex_type *f = (complex_type *)(field_elem) + node_column_size * elements * l;
        for (int s = 0; s < nnodes; s++) {
            complex_type *f_sl = f + cpn * node_column_size * elements * s;
            for (int e = 0; e < elements; e++) { // Complex elements / field element
                d[t + node_column_size * s + e * column_size].x =
                    f_sl[e + elements * t].re;
                d[t + node_column_size * s + e * column_size].y =
                    f_sl[e + elements * t].im;
            }
        }
    }
}




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

    // alloc work buffers
    char * collect_buffer, * receive_buffer;
    cudaMalloc((void **)&collect_buffer, local_volume * elements * sizeof(cmplx_t));
    cudaMalloc((void **)&receive_buffer, local_volume * elements * sizeof(cmplx_t));

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




#endif