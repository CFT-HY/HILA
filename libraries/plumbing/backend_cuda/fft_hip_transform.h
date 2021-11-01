#ifndef FFT_HIP_TRANSFORM_H
#define FFT_HIP_TRANSFORM_H

#ifndef HILAPP

#if defined(CUDA)

#include <cufft.h>

using gpufftComplex = cufftComplex;
using gpufftDoubleComplex = cufftDoubleComplex;
using gpufftHandle = cufftHandle;
#define gpufftExecC2C cufftExecC2C
#define gpufftExecZ2Z cufftExecZ2Z
#define gpufftPlan1d cufftPlan1d
#define gpufftDestroy cufftDestroy

#define GPUFFT_FORWARD CUFFT_FORWARD
#define GPUFFT_INVERSE CUFFT_INVERSE

#define GPUFFT_C2C CUFFT_C2C
#define GPUFFT_Z2Z CUFFT_Z2Z

#else


#include "hip/hip_runtime.h"
#include <hipfft.h>

using gpufftComplex = hipfftComplex;
using gpufftDoubleComplex = hipfftDoubleComplex;
using gpufftHandle = hipfftHandle;
#define gpufftExecC2C hipfftExecC2C
#define gpufftExecZ2Z hipfftExecZ2Z
#define gpufftPlan1d hipfftPlan1d
#define gpufftDestroy hipfftDestroy

#define GPUFFT_FORWARD HIPFFT_FORWARD
#define GPUFFT_INVERSE HIPFFT_BACKWARD

#define GPUFFT_C2C HIPFFT_C2C
#define GPUFFT_Z2Z HIPFFT_Z2Z

#endif

#ifndef HILAPP

/// Gather one element column from the mpi buffer
template <typename cmplx_t>
__global__ void gather_column(cmplx_t *RESTRICT data, cmplx_t *RESTRICT *d_ptr,
                              int *RESTRICT d_size, int n, int colsize, int columns) {

    int ind = threadIdx.x + blockIdx.x * blockDim.x;
    if (ind < columns) {
        int s = colsize * ind;

        int k = s;
        for (int i = 0; i < n; i++) {
            int offset = ind * d_size[i];
            for (int j = 0; j < d_size[i]; j++, k++) {
                data[k] = d_ptr[i][j + offset];
            }
        }
    }
}

/// Gather one element column from the mpi buffer
template <typename cmplx_t>
__global__ void scatter_column(cmplx_t *RESTRICT data, cmplx_t *RESTRICT *d_ptr,
                               int *RESTRICT d_size, int n, int colsize, int columns) {

    int ind = threadIdx.x + blockIdx.x * blockDim.x;
    if (ind < columns) {
        int s = colsize * ind;

        int k = s;
        for (int i = 0; i < n; i++) {
            int offset = ind * d_size[i];
            for (int j = 0; j < d_size[i]; j++, k++) {
                d_ptr[i][j + offset] = data[k];
            }
        }
    }
}

/// Define appropriate cufft-type depending on the cmplx_t -type
/// these types are 1-1 compatible anyway
/// use as cufft_cmplx_t<T>::type
template <typename cmplx_t>
using fft_cmplx_t = typename std::conditional<sizeof(gpufftComplex) == sizeof(cmplx_t),
                                              gpufftComplex, gpufftDoubleComplex>::type;

/// Templates for cufftExec float and double complex

template <typename cmplx_t,
          std::enable_if_t<sizeof(cmplx_t) == sizeof(gpufftComplex), int> = 0>
inline void gpufft_execute(gpufftHandle plan, cmplx_t *buf, int direction) {
    gpufftExecC2C(plan, (gpufftComplex *)buf, (gpufftComplex *)buf, direction);
}

template <typename cmplx_t,
          std::enable_if_t<sizeof(cmplx_t) == sizeof(gpufftDoubleComplex), int> = 0>
inline void gpufft_execute(gpufftHandle plan, cmplx_t *buf, int direction) {
    gpufftExecZ2Z(plan, (gpufftDoubleComplex *)buf, (gpufftDoubleComplex *)buf,
                  direction);
}

template <typename cmplx_t>
void hila_fft<cmplx_t>::transform() {

    extern unsigned hila_fft_my_columns[NDIM];
    extern hila::timer fft_plan_timer, fft_execute_timer, fft_buffer_timer;

    constexpr bool is_float = (sizeof(cmplx_t) == sizeof(Complex<float>));

    int n_columns = hila_fft_my_columns[dir] * elements;

    int direction =
        (fftdir == fft_direction::forward) ? GPUFFT_FORWARD : GPUFFT_INVERSE;


    // allocate here fftw plans.  TODO: perhaps store, if plans take appreciable time?
    // Timer will tell the proportional timing

    gpufftHandle plan;
    int batch = hila_fft_my_columns[dir];
    int n_fft = 1;
    // reduce very large batch to smaller, avoid large buffer space

    fft_plan_timer.start();


    bool is_divisible = true;
    while (batch > GPUFFT_BATCH_SIZE && is_divisible) {
        is_divisible = false;
        for (int div : {2, 3, 5, 7}) {
            if (batch % div == 0) {
                batch /= div;
                n_fft *= div;
                is_divisible = true;
                break;
            }
        }
    }

    // output0 << " Batch " << batch << " nfft " << n_fft << '\n';

    // HIPFFT_C2C for float transform, Z2Z for double
    gpufftPlan1d(&plan, lattice->size(dir), is_float ? GPUFFT_C2C : GPUFFT_Z2Z, batch);
    check_device_error("FFT plan");

    fft_plan_timer.stop();

    // alloc work array
    cmplx_t *fft_wrk = (cmplx_t *)d_malloc(buf_size * sizeof(cmplx_t) * elements);

    // Reorganize the data to form columns of a single element
    // move from receive_buf to fft_wrk
    // first need to copy index arrays to device

    fft_buffer_timer.start();

    cmplx_t **d_ptr = (cmplx_t **)d_malloc(sizeof(cmplx_t *) * rec_p.size());
    int *d_size = (int *)d_malloc(sizeof(int) * rec_p.size());

    gpuMemcpy(d_ptr, rec_p.data(), rec_p.size() * sizeof(cmplx_t *),
              gpuMemcpyHostToDevice);
    gpuMemcpy(d_size, rec_size.data(), rec_size.size() * sizeof(int),
              gpuMemcpyHostToDevice);

    int N_blocks = (n_columns + N_threads - 1) / N_threads;

#if defined(CUDA)
    gather_column<cmplx_t><<<N_blocks, N_threads>>>(
        fft_wrk, d_ptr, d_size, rec_p.size(), lattice->size(dir), n_columns);
#else
    hipLaunchKernelGGL(HIP_KERNEL_NAME(gather_column<cmplx_t>), dim3(N_blocks),
                       dim3(N_threads), 0, 0, fft_wrk, d_ptr, d_size, rec_p.size(),
                       lattice->size(dir), n_columns);
#endif

    fft_buffer_timer.stop();

    // do the fft
    fft_execute_timer.start();

    for (int i = 0; i < n_fft; i++) {

        cmplx_t *cp = fft_wrk + i * (batch * lattice->size(dir));

        gpufft_execute(plan, cp, direction);
        check_device_error("FFT execute");
    }

    fft_execute_timer.stop();

    fft_buffer_timer.start();

#if defined(CUDA)
    scatter_column<cmplx_t><<<N_blocks, N_threads>>>(
        fft_wrk, d_ptr, d_size, rec_p.size(), lattice->size(dir), n_columns);
#else
    hipLaunchKernelGGL(HIP_KERNEL_NAME(scatter_column<cmplx_t>), dim3(N_blocks),
                       dim3(N_threads), 0, 0, fft_wrk, d_ptr, d_size, rec_p.size(),
                       lattice->size(dir), n_columns);
#endif

    fft_buffer_timer.stop();

    d_free(d_size);
    d_free(d_ptr);
    d_free(fft_wrk);

    gpufftDestroy(plan);
}

#endif // HILAPP

////////////////////////////////////////////////////////////////////
/// send column data to nodes

template <typename cmplx_t>
void hila_fft<cmplx_t>::gather_data() {

#ifdef USE_MPI

    extern hila::timer fft_MPI_timer;
    fft_MPI_timer.start();

    // post receive and send
    int n_comms = hila_fft_comms[dir].size() - 1;

    MPI_Request sendreq[n_comms], recreq[n_comms];
    MPI_Status stat[n_comms];

#ifndef CUDA_AWARE_MPI
    cmplx_t *send_p[n_comms];
    cmplx_t *receive_p[n_comms];
#endif

    int i = 0;
    int j = 0;
    for (auto &fn : hila_fft_comms[dir]) {
        if (fn.node != hila::myrank()) {

            size_t siz = fn.recv_buf_size * elements;
            if (siz >= (1ULL << 30)) {
                hila::output << "Too large MPI message in FFT! Size " << siz
                             << " complex numbers\n";
                hila::terminate(1);
            }

#ifndef CUDA_AWARE_MPI
            cmplx_t *p = receive_p[i] = (cmplx_t *)memalloc(sizeof(cmplx_t) * siz);
#else
            cmplx_t *p = rec_p[j];
#endif

            MPI_Irecv(p, (int)siz, mpi_cmplx_t, fn.node, WRK_GATHER_TAG,
                      lattice->mpi_comm_lat, &recreq[i]);

            i++;
        }
        j++;
    }

    i = 0;
    for (auto &fn : hila_fft_comms[dir]) {
        if (fn.node != hila::myrank()) {

            cmplx_t *p = send_buf + fn.column_offset * elements;
            int n = fn.column_number * elements * lattice->mynode.size[dir];

#ifndef CUDA_AWARE_MPI
            send_p[i] = (cmplx_t *)memalloc(sizeof(cmplx_t) * n);
            gpuMemcpy(send_p[i], p, sizeof(cmplx_t) * n, gpuMemcpyDeviceToHost);
            p = send_p[i];
#endif

            MPI_Isend(p, n, mpi_cmplx_t, fn.node, WRK_GATHER_TAG, lattice->mpi_comm_lat,
                      &sendreq[i]);
            i++;
        }
    }

    // and wait for the send and receive to complete
    if (n_comms > 0) {
        MPI_Waitall(n_comms, recreq, stat);
        MPI_Waitall(n_comms, sendreq, stat);

#ifndef CUDA_AWARE_MPI
        i = j = 0;

        for (auto &fn : hila_fft_comms[dir]) {
            if (fn.node != hila::myrank()) {

                size_t siz = fn.recv_buf_size * elements;

                gpuMemcpy(rec_p[j], receive_p[i], siz * sizeof(cmplx_t),
                          gpuMemcpyHostToDevice);
                i++;
            }
            j++;
        }

        for (i = 0; i < n_comms; i++) {
            free(receive_p[i]);
            free(send_p[i]);
        }
#endif
    }

    fft_MPI_timer.stop();

#endif
}

//////////////////////////////////////////////////////////////////////////////////////
/// inverse of gather_data

template <typename cmplx_t>
void hila_fft<cmplx_t>::scatter_data() {

#ifdef USE_MPI

    extern hila::timer fft_MPI_timer;
    fft_MPI_timer.start();

    int n_comms = hila_fft_comms[dir].size() - 1;

    MPI_Request sendreq[n_comms], recreq[n_comms];
    MPI_Status stat[n_comms];

#ifndef CUDA_AWARE_MPI
    cmplx_t *send_p[n_comms];
    cmplx_t *receive_p[n_comms];
#endif

    int i = 0;

    for (auto &fn : hila_fft_comms[dir]) {
        if (fn.node != hila::myrank()) {

            int n = fn.column_number * elements * lattice->mynode.size[dir];
#ifdef CUDA_AWARE_MPI
            cmplx_t *p = send_buf + fn.column_offset * elements;
#else
            cmplx_t *p = receive_p[i] = (cmplx_t *)memalloc(sizeof(cmplx_t) * n);
#endif

            MPI_Irecv(p, n, mpi_cmplx_t, fn.node, WRK_SCATTER_TAG,
                      lattice->mpi_comm_lat, &recreq[i]);

            i++;
        }
    }

    i = 0;
    int j = 0;
    for (auto &fn : hila_fft_comms[dir]) {
        if (fn.node != hila::myrank()) {

            int n = fn.recv_buf_size * elements;
#ifndef CUDA_AWARE_MPI
            cmplx_t *p = send_p[i] = (cmplx_t *)memalloc(sizeof(cmplx_t) * n);
            gpuMemcpy(p, rec_p[j], sizeof(cmplx_t) * n, gpuMemcpyDeviceToHost);
#else
            cmplx_t *p = rec_p[j];
#endif
            MPI_Isend(p, n, mpi_cmplx_t, fn.node, WRK_SCATTER_TAG,
                      lattice->mpi_comm_lat, &sendreq[i]);

            i++;
        }
        j++;
    }

    // and wait for the send and receive to complete
    if (n_comms > 0) {
        MPI_Waitall(n_comms, recreq, stat);
        MPI_Waitall(n_comms, sendreq, stat);

#ifndef CUDA_AWARE_MPI
        i = 0;
        for (auto &fn : hila_fft_comms[dir]) {
            if (fn.node != hila::myrank()) {

                int n = fn.column_number * elements * lattice->mynode.size[dir];
                cmplx_t *p = send_buf + fn.column_offset * elements;

                gpuMemcpy(p, receive_p[i], n * sizeof(cmplx_t), gpuMemcpyHostToDevice);
                i++;
            }
        }

        for (i = 0; i < n_comms; i++) {
            free(receive_p[i]);
            free(send_p[i]);
        }
#endif
    }

    fft_MPI_timer.stop();
#endif
}

#endif // HILAPP

#endif