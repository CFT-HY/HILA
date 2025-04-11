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
#if defined(HIPCUB_VERSION) && HIPCUB_VERSION >= 300300
// For some strange reason hipcub 3.3 used long in KeyValuePair ?? bug in library?
using keyvalueindexT = long;
#else
using keyvalueindexT = int;
#endif
#endif  // HIP


template <typename T>
T Field<T>::gpu_minmax(bool is_min, Parity par, CoordinateVector &loc) const {

    // skip the cub/hipcub bits in hilapp, not needed

    int64_t num_items = lattice.loop_end(par) - lattice.loop_begin(par);

    // Declare, allocate, and initialize device-accessible pointers
    // for input and output
    T *data_in = this->field_buffer() + lattice.loop_begin(par); // ptr to data
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
    // hila::out0 << "gpu_minmax: alloc " << temp_storage_bytes << " bytes\n";
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

    gpuMemcpy(&result, result_p, sizeof(gpucub::KeyValuePair<keyvalueindexT, T>), gpuMemcpyDeviceToHost);

    gpuFree(result_p);

    loc = lattice.coordinates(result.key + lattice.loop_begin(par));
    return result.value;
}

#endif  // ! HILAPP


#endif