#ifndef GPU_MINMAX_H_
#define GPU_MINMAX_H_

#if !defined(HILAPP)
// hilapp does not have to read through this, nothing to convert
// and hilapp does not understand .cuh -files

#include "hila.h"

#if defined(CUDA)
#include <cub/cub.cuh>
namespace gpucub = cub;
using keyvalueindexT = int;
#endif

#if defined(HIP)
#include <hipcub/hipcub.hpp>
namespace gpucub = hipcub;
#if defined(HIPCUB_VERSION) && HIPCUB_VERSION >= 300200
// For some strange reason hipcub 3.3 used long in KeyValuePair ?? bug in library?
using keyvalueindexT = long;
#else
using keyvalueindexT = int;
#endif
#endif // HIP


template <typename T>
gpucub::KeyValuePair<keyvalueindexT, T> gpu_minmax_argreduce_range(bool is_min, const T *data_in,
                                                                    int64_t num_items) {
    gpucub::KeyValuePair<keyvalueindexT, T> *result_p, result;

    gpuMalloc(&result_p, sizeof(gpucub::KeyValuePair<keyvalueindexT, T>));

    // Determine temporary device storage requirements
    void *d_temp_storage = nullptr;
    size_t temp_storage_bytes = 0;

    if (is_min) {
        GPU_CHECK(gpucub::DeviceReduce::ArgMin(d_temp_storage, temp_storage_bytes, data_in,
                                               result_p, num_items));
    } else {
        GPU_CHECK(gpucub::DeviceReduce::ArgMax(d_temp_storage, temp_storage_bytes, data_in,
                                               result_p, num_items));
    }

    // Allocate temporary storage
    gpuMalloc(&d_temp_storage, temp_storage_bytes);

    // Run argmin-reduction
    if (is_min) {
        GPU_CHECK(gpucub::DeviceReduce::ArgMin(d_temp_storage, temp_storage_bytes, data_in,
                                               result_p, num_items));
    } else {
        GPU_CHECK(gpucub::DeviceReduce::ArgMax(d_temp_storage, temp_storage_bytes, data_in,
                                               result_p, num_items));
    }

    gpuFree(d_temp_storage);

    gpuMemcpy(&result, result_p, sizeof(gpucub::KeyValuePair<keyvalueindexT, T>),
              gpuMemcpyDeviceToHost);

    gpuFree(result_p);

    return result;
}

template <typename T>
T Field<T>::gpu_minmax(bool is_min, Parity par, CoordinateVector &loc) const {


    const lattice_struct &mylat = fs->mylattice.ref();

#ifdef GPU_OVERLAP_COMM
    hila::iter_range_t _hila_ranges;
    int _hila_loops = mylat.loop_ranges(par, 0, _hila_ranges);

    T *base = this->field_buffer();
    auto result =
        gpu_minmax_argreduce_range(is_min, base + _hila_ranges.min[0],
                                   _hila_ranges.max[0] - _hila_ranges.min[0]);
    int loop_begin = _hila_ranges.min[0];

    if (_hila_loops == 2) {
        auto result2 =
            gpu_minmax_argreduce_range(is_min, base + _hila_ranges.min[1],
                                       _hila_ranges.max[1] - _hila_ranges.min[1]);
        bool take_2 = is_min ? (result2.value < result.value) : (result2.value > result.value);
        if (take_2) {
            result = result2;
            loop_begin = _hila_ranges.min[1];
        }
    }

    loc = mylat.coordinates(result.key + loop_begin);
    return result.value;
#else
    int _hila_loop_begin = mylat.loop_begin(par);
    int _hila_loop_end = mylat.loop_end(par);
    int64_t num_items = _hila_loop_end - _hila_loop_begin;

    T *data_in = this->field_buffer() + _hila_loop_begin; // ptr to data
    auto result = gpu_minmax_argreduce_range(is_min, data_in, num_items);

    loc = mylat.coordinates(result.key + _hila_loop_begin);
    return result.value;
#endif
}

#endif // ! HILAPP


#endif