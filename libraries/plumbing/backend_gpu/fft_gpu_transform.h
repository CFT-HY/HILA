#ifndef FFT_GPU_TRANSFORM_H
#define FFT_GPU_TRANSFORM_H

#ifndef HILAPP

#include "plumbing/fft_structs.h"

#if defined(CUDA)

#include <cufft.h>

using gpufftComplex = cufftComplex;
using gpufftDoubleComplex = cufftDoubleComplex;
using gpufftHandle = cufftHandle;
using gpufftResult = cufftResult;
#define gpufftExecC2C cufftExecC2C
#define gpufftExecZ2Z cufftExecZ2Z
#define gpufftPlan1d cufftPlan1d
#define gpufftDestroy cufftDestroy
#define gpufftCreate cufftCreate
#define gpufftSetAutoAllocation cufftSetAutoAllocation
#define gpufftSetWorkArea cufftSetWorkArea
#define gpufftMakePlan1d cufftMakePlan1d
#define gpufftGetSize1d cufftGetSize1d
#define gpufftEstimate1d cufftEstimate1d

#define GPUFFT_FORWARD CUFFT_FORWARD
#define GPUFFT_INVERSE CUFFT_INVERSE

#define GPUFFT_C2C CUFFT_C2C
#define GPUFFT_Z2Z CUFFT_Z2Z

#else

#include "hip/hip_runtime.h"
#include <hipfft/hipfft.h>

using gpufftComplex = hipfftComplex;
using gpufftDoubleComplex = hipfftDoubleComplex;
using gpufftHandle = hipfftHandle;
#define gpufftExecC2C hipfftExecC2C
#define gpufftExecZ2Z hipfftExecZ2Z
#define gpufftPlan1d hipfftPlan1d
#define gpufftDestroy hipfftDestroy
#define gpufftCreate hipfftCreate
#define gpufftSetAutoAllocation hipfftSetAutoAllocation
#define gpufftSetWorkArea hipfftSetWorkArea
#define gpufftMakePlan1d hipfftMakePlan1d
#define gpufftGetSize1d hipfftGetSize1d
#define gpufftEstimate1d hipfftEstimate1d

#define GPUFFT_FORWARD HIPFFT_FORWARD
#define GPUFFT_INVERSE HIPFFT_BACKWARD

#define GPUFFT_C2C HIPFFT_C2C
#define GPUFFT_Z2Z HIPFFT_Z2Z

#endif


/// Gather one element column from the mpi buffer
template <typename cmplx_t>
__global__ void hila_fft_gather_column(cmplx_t *RESTRICT data, cmplx_t *RESTRICT *d_ptr,
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
__global__ void hila_fft_scatter_column(cmplx_t *RESTRICT data, cmplx_t *RESTRICT *d_ptr,
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

// Define datatype for saved plans

#define GPUFFT_SHARE_PLAN_MEMORY

class hila_saved_fftplan_t {
  public:
    struct plan_d {
        gpufftHandle plan;
        int size;
        int batch;
        bool is_float;
    };

    std::vector<plan_d> plans;

#if defined(GPUFFT_SHARE_PLAN_MEMORY)
    void *work_area = nullptr;
    size_t work_area_size = 0;
#endif

    unsigned long seq = 0;


    hila_saved_fftplan_t() {}

    ~hila_saved_fftplan_t() {
        delete_plans();
    }

    void delete_plans() {
        static bool deleted = false;

        if (deleted)
            return;
        deleted = true;

        for (auto &p : plans) {
            gpufftDestroy(p.plan);
        }
#if defined(GPUFFT_SHARE_PLAN_MEMORY)
        if (work_area != nullptr)
            d_free(work_area);
        hila::out0 << "GPU FFT Work area size on node 0: " << work_area_size << '\n';
#endif
        hila::out0 << "GPU FFT Plans: " << plans.size() << " plans used\n";
        plans.clear();
    }

    // get cached plan or create new.  If the saved plan is incompatible with
    // the one required, destroy plans
    gpufftHandle get_plan(int size, int batch, bool is_float) {

        extern hila::timer fft_plan_timer;

        // do we have saved plan of the same type?

        seq++;

        for (auto &p : plans) {
            if (p.size == size && p.batch == batch && p.is_float == is_float) {
                // Now we got it!
                return p.plan;
            }
        }

        // not cached, make new

        fft_plan_timer.start();

        // new plan hook
        plans.emplace_back();
        plan_d *pp = &plans.back();

        pp->size = size;
        pp->batch = batch;
        pp->is_float = is_float;

#if defined(GPUFFT_SHARE_PLAN_MEMORY)
        
        gpufftCreate(&pp->plan);
        // turn off auto allocate
        gpufftSetAutoAllocation(pp->plan, 0);
        // and get the required scratch size
        size_t worksize;
        gpufftEstimate1d(size, is_float ? GPUFFT_C2C : GPUFFT_Z2Z, batch, &worksize);
        // gpufftGetSize1d(pp->plan, size, is_float ? GPUFFT_C2C : GPUFFT_Z2Z, batch, &worksize);
        check_device_error("FFT Plan size");

        if (work_area_size >= worksize) {
            // use existing area
            gpufftSetWorkArea(pp->plan, work_area);
            check_device_error("FFT Plan work area set");

        } else {
            // (re-)allocate the work area
            if (work_area != nullptr) {
                d_free(work_area);
            }

            hila::out0 << "FFT: allocate new work area of size " << worksize << '\n';

            work_area_size = worksize;
            work_area = (void *)d_malloc(work_area_size);

            for (auto &pref : plans) {
                gpufftSetWorkArea(pref.plan, work_area);
                check_device_error("FFT Plan work area reset");
            }
        }

        auto siz = work_area_size;
        gpufftMakePlan1d(pp->plan, size, is_float ? GPUFFT_C2C : GPUFFT_Z2Z, batch,
                         &siz);
        check_device_error("FFT Plan make1d");
        assert(work_area_size >= siz && "GPU work area mismatch");

#else

        // _C2C for float transform, Z2Z for double

        gpufftPlan1d(&(pp->plan), size, is_float ? GPUFFT_C2C : GPUFFT_Z2Z, batch);
        check_device_error("FFT plan");

#endif

        fft_plan_timer.stop();

        return pp->plan;
    }
};

/// Define appropriate cufft-type depending on the cmplx_t -type
/// these types are 1-1 compatible anyway
/// use as cufft_cmplx_t<T>::type
template <typename cmplx_t>
using fft_cmplx_t = typename std::conditional<sizeof(gpufftComplex) == sizeof(cmplx_t),
                                              gpufftComplex, gpufftDoubleComplex>::type;

/// Templates for cufftExec float and double complex

template <typename cmplx_t, std::enable_if_t<sizeof(cmplx_t) == sizeof(gpufftComplex), int> = 0>
inline void hila_gpufft_execute(gpufftHandle plan, cmplx_t *buf, int direction) {
    gpufftExecC2C(plan, (gpufftComplex *)buf, (gpufftComplex *)buf, direction);
}

template <typename cmplx_t,
          std::enable_if_t<sizeof(cmplx_t) == sizeof(gpufftDoubleComplex), int> = 0>
inline void hila_gpufft_execute(gpufftHandle plan, cmplx_t *buf, int direction) {
    gpufftExecZ2Z(plan, (gpufftDoubleComplex *)buf, (gpufftDoubleComplex *)buf, direction);
}

template <typename cmplx_t>
void hila_fft<cmplx_t>::transform() {

    // these externs defined in fft.cpp
    extern hila::timer fft_execute_timer, fft_buffer_timer;
    extern hila_saved_fftplan_t hila_saved_fftplan;

    constexpr bool is_float = (sizeof(cmplx_t) == sizeof(Complex<float>));

    const hila::fftdata_struct &fft = *(lattice->fftdata);

    int n_columns = fft.hila_fft_my_columns[dir] * elements;

    int direction = (fftdir == fft_direction::forward) ? GPUFFT_FORWARD : GPUFFT_INVERSE;

    // allocate here fftw plans.  TODO: perhaps store, if plans take appreciable time?
    // Timer will tell the proportional timing

    int batch = fft.hila_fft_my_columns[dir];
    int n_fft = elements;
    // reduce very large batch to smaller, avoid large buffer space

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

    gpufftHandle plan = hila_saved_fftplan.get_plan(lattice.size(dir), batch, is_float);
    // hila::out0 << " Batch " << batch << " nfft " << n_fft << '\n';

    // alloc work array
    cmplx_t *fft_wrk = (cmplx_t *)d_malloc(buf_size * sizeof(cmplx_t) * elements);

    // Reorganize the data to form columns of a single element
    // move from receive_buf to fft_wrk
    // first need to copy index arrays to device

    fft_buffer_timer.start();

    cmplx_t **d_ptr = (cmplx_t **)d_malloc(sizeof(cmplx_t *) * rec_p.size());
    int *d_size = (int *)d_malloc(sizeof(int) * rec_p.size());

    gpuMemcpy(d_ptr, rec_p.data(), rec_p.size() * sizeof(cmplx_t *), gpuMemcpyHostToDevice);
    gpuMemcpy(d_size, rec_size.data(), rec_size.size() * sizeof(int), gpuMemcpyHostToDevice);

    int N_blocks = (n_columns + N_threads - 1) / N_threads;

#if defined(CUDA)
    hila_fft_gather_column<cmplx_t><<<N_blocks, N_threads>>>(fft_wrk, d_ptr, d_size, rec_p.size(),
                                                             lattice.size(dir), n_columns);
#else
    hipLaunchKernelGGL(HIP_KERNEL_NAME(hila_fft_gather_column<cmplx_t>), dim3(N_blocks),
                       dim3(N_threads), 0, 0, fft_wrk, d_ptr, d_size, rec_p.size(),
                       lattice.size(dir), n_columns);
#endif

    fft_buffer_timer.stop();

    // do the fft
    fft_execute_timer.start();

    for (int i = 0; i < n_fft; i++) {

        cmplx_t *cp = fft_wrk + i * (batch * lattice.size(dir));

        hila_gpufft_execute(plan, cp, direction);
        check_device_error("FFT execute");
    }

    fft_execute_timer.stop();

    fft_buffer_timer.start();

#if defined(CUDA)
    hila_fft_scatter_column<cmplx_t><<<N_blocks, N_threads>>>(fft_wrk, d_ptr, d_size, rec_p.size(),
                                                              lattice.size(dir), n_columns);
#else
    hipLaunchKernelGGL(HIP_KERNEL_NAME(hila_fft_scatter_column<cmplx_t>), dim3(N_blocks),
                       dim3(N_threads), 0, 0, fft_wrk, d_ptr, d_size, rec_p.size(),
                       lattice.size(dir), n_columns);
#endif

    fft_buffer_timer.stop();

    d_free(d_size);
    d_free(d_ptr);
    d_free(fft_wrk);
}


////////////////////////////////////////////////////////////////////
/// send column data to nodes

template <typename cmplx_t>
void hila_fft<cmplx_t>::gather_data() {


    extern hila::timer pencil_MPI_timer;
    pencil_MPI_timer.start();

    const hila::fftdata_struct &fft = *(lattice->fftdata);

    // post receive and send
    int n_comms = fft.hila_pencil_comms[dir].size() - 1;

    std::vector<MPI_Request> sendreq(n_comms), recreq(n_comms);
    std::vector<MPI_Status> stat(n_comms);

#ifndef GPU_AWARE_MPI
    std::vector<cmplx_t *> send_p(n_comms);
    std::vector<cmplx_t *> receive_p(n_comms);
#endif

    int i = 0;
    int j = 0;

    // this synchronization should be enough for all MPI's in the
    gpuStreamSynchronize(0);

    size_t mpi_type_size;
    MPI_Datatype mpi_type = get_MPI_complex_type<cmplx_t>(mpi_type_size);

    for (auto &fn : fft.hila_pencil_comms[dir]) {
        if (fn.node != hila::myrank()) {

            size_t siz = fn.recv_buf_size * elements * sizeof(cmplx_t);
            if (siz >= (1ULL << 31) * mpi_type_size) {
                hila::out << "Too large MPI message in pencils! Size " << siz << " bytes ("
                          << siz / mpi_type_size << " elements)\n";
                hila::terminate(1);
            }

#ifndef GPU_AWARE_MPI
            cmplx_t *p = receive_p[i] = (cmplx_t *)memalloc(siz);
#else
            cmplx_t *p = rec_p[j];
#endif

            MPI_Irecv(p, (int)(siz / mpi_type_size), mpi_type, fn.node, WRK_GATHER_TAG,
                      lattice->mpi_comm_lat, &recreq[i]);

            i++;
        }
        j++;
    }

    i = 0;
    for (auto &fn : fft.hila_pencil_comms[dir]) {
        if (fn.node != hila::myrank()) {

            cmplx_t *p = send_buf + fn.column_offset * elements;
            size_t n = fn.column_number * elements * lattice->mynode.size[dir] * sizeof(cmplx_t);

#ifndef GPU_AWARE_MPI
            // now not GPU_AWARE_MPI
            send_p[i] = (cmplx_t *)memalloc(n);
            gpuMemcpy(send_p[i], p, n, gpuMemcpyDeviceToHost);
            p = send_p[i];
#endif

            MPI_Isend(p, (int)(n / mpi_type_size), mpi_type, fn.node, WRK_GATHER_TAG,
                      lattice->mpi_comm_lat, &sendreq[i]);
            i++;
        }
    }

    // and wait for the send and receive to complete
    if (n_comms > 0) {
        MPI_Waitall(n_comms, recreq.data(), stat.data());
        MPI_Waitall(n_comms, sendreq.data(), stat.data());

#ifndef GPU_AWARE_MPI
        i = j = 0;

        for (auto &fn : fft.hila_pencil_comms[dir]) {
            if (fn.node != hila::myrank()) {

                size_t siz = fn.recv_buf_size * elements;

                gpuMemcpy(rec_p[j], receive_p[i], siz * sizeof(cmplx_t), gpuMemcpyHostToDevice);
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

    pencil_MPI_timer.stop();
}

//////////////////////////////////////////////////////////////////////////////////////
/// inverse of gather_data

template <typename cmplx_t>
void hila_fft<cmplx_t>::scatter_data() {


    extern hila::timer pencil_MPI_timer;
    pencil_MPI_timer.start();

    const hila::fftdata_struct &fft = *(lattice->fftdata);

    int n_comms = fft.hila_pencil_comms[dir].size() - 1;

    std::vector<MPI_Request> sendreq(n_comms), recreq(n_comms);
    std::vector<MPI_Status> stat(n_comms);

#ifndef GPU_AWARE_MPI
    std::vector<cmplx_t *> send_p(n_comms);
    std::vector<cmplx_t *> receive_p(n_comms);
#endif

    int i = 0;

    size_t mpi_type_size;
    MPI_Datatype mpi_type = get_MPI_complex_type<cmplx_t>(mpi_type_size);

    gpuStreamSynchronize(0);

    for (auto &fn : fft.hila_pencil_comms[dir]) {
        if (fn.node != hila::myrank()) {

            size_t n = fn.column_number * elements * lattice->mynode.size[dir] * sizeof(cmplx_t);
#ifdef GPU_AWARE_MPI
            cmplx_t *p = send_buf + fn.column_offset * elements;
#else
            cmplx_t *p = receive_p[i] = (cmplx_t *)memalloc(n);
#endif

            MPI_Irecv(p, (int)(n / mpi_type_size), mpi_type, fn.node, WRK_SCATTER_TAG,
                      lattice->mpi_comm_lat, &recreq[i]);

            i++;
        }
    }

    i = 0;
    int j = 0;
    for (auto &fn : fft.hila_pencil_comms[dir]) {
        if (fn.node != hila::myrank()) {

            size_t n = fn.recv_buf_size * elements * sizeof(cmplx_t);
#ifdef GPU_AWARE_MPI
            cmplx_t *p = rec_p[j];
//             gpuStreamSynchronize(0);
#else
            cmplx_t *p = send_p[i] = (cmplx_t *)memalloc(n);
            gpuMemcpy(p, rec_p[j], n, gpuMemcpyDeviceToHost);
#endif
            MPI_Isend(p, (int)(n / mpi_type_size), mpi_type, fn.node, WRK_SCATTER_TAG,
                      lattice->mpi_comm_lat, &sendreq[i]);

            i++;
        }
        j++;
    }

    // and wait for the send and receive to complete
    if (n_comms > 0) {
        MPI_Waitall(n_comms, recreq.data(), stat.data());
        MPI_Waitall(n_comms, sendreq.data(), stat.data());

#ifndef GPU_AWARE_MPI
        i = 0;
        for (auto &fn : fft.hila_pencil_comms[dir]) {
            if (fn.node != hila::myrank()) {

                size_t n =
                    fn.column_number * elements * lattice->mynode.size[dir] * sizeof(cmplx_t);
                cmplx_t *p = send_buf + fn.column_offset * elements;

                gpuMemcpy(p, receive_p[i], n, gpuMemcpyHostToDevice);
                i++;
            }
        }

        for (i = 0; i < n_comms; i++) {
            free(receive_p[i]);
            free(send_p[i]);
        }
#endif
    }

    pencil_MPI_timer.stop();
}

///////////////////////////////////////////////////////////////////////////////////
/// Separate reflect operation
/// Reflect flips the coordinates so that negative direction becomes positive,
/// and x=0 plane remains,
/// r(x) <- f(L - x)  -  note that x == 0 layer is as before
/// r(0) = f(0), r(1) = f(L-1), r(2) = f(L-2)  ...
///////////////////////////////////////////////////////////////////////////////////

/// Reflect data in the array
template <typename cmplx_t>
__global__ void hila_reflect_dir_kernel(cmplx_t *RESTRICT data, const int colsize,
                                        const int columns) {

    int ind = threadIdx.x + blockIdx.x * blockDim.x;
    if (ind < columns) {
        const int s = colsize * ind;

        for (int i = 1; i < colsize / 2; i++) {
            int i1 = s + i;
            int i2 = s + colsize - i;
            cmplx_t tmp = data[i1];
            data[i1] = data[i2];
            data[i2] = tmp;
        }
    }
}


template <typename cmplx_t>
void hila_fft<cmplx_t>::reflect() {


    constexpr bool is_float = (sizeof(cmplx_t) == sizeof(Complex<float>));

    const hila::fftdata_struct &fft = *(lattice->fftdata);

    int n_columns = fft.hila_fft_my_columns[dir] * elements;

    // reduce very large batch to smaller, avoid large buffer space

    // alloc work array
    cmplx_t *fft_wrk = (cmplx_t *)d_malloc(buf_size * sizeof(cmplx_t) * elements);

    // Reorganize the data to form columns of a single element
    // move from receive_buf to fft_wrk
    // first need to copy index arrays to device

    cmplx_t **d_ptr = (cmplx_t **)d_malloc(sizeof(cmplx_t *) * rec_p.size());
    int *d_size = (int *)d_malloc(sizeof(int) * rec_p.size());

    gpuMemcpy(d_ptr, rec_p.data(), rec_p.size() * sizeof(cmplx_t *), gpuMemcpyHostToDevice);
    gpuMemcpy(d_size, rec_size.data(), rec_size.size() * sizeof(int), gpuMemcpyHostToDevice);

    int N_blocks = (n_columns + N_threads - 1) / N_threads;

#if defined(CUDA)
    hila_fft_gather_column<cmplx_t><<<N_blocks, N_threads>>>(fft_wrk, d_ptr, d_size, rec_p.size(),
                                                             lattice.size(dir), n_columns);
#else
    hipLaunchKernelGGL(HIP_KERNEL_NAME(hila_fft_gather_column<cmplx_t>), dim3(N_blocks),
                       dim3(N_threads), 0, 0, fft_wrk, d_ptr, d_size, rec_p.size(),
                       lattice.size(dir), n_columns);
#endif

#if defined(CUDA)
    hila_reflect_dir_kernel<cmplx_t>
        <<<N_blocks, N_threads>>>(fft_wrk, lattice.size(dir), n_columns);
#else
    hipLaunchKernelGGL(HIP_KERNEL_NAME(hila_reflect_dir_kernel<cmplx_t>), dim3(N_blocks),
                       dim3(N_threads), 0, 0, fft_wrk, lattice.size(dir), n_columns);
#endif


#if defined(CUDA)
    hila_fft_scatter_column<cmplx_t><<<N_blocks, N_threads>>>(fft_wrk, d_ptr, d_size, rec_p.size(),
                                                              lattice.size(dir), n_columns);
#else
    hipLaunchKernelGGL(HIP_KERNEL_NAME(hila_fft_scatter_column<cmplx_t>), dim3(N_blocks),
                       dim3(N_threads), 0, 0, fft_wrk, d_ptr, d_size, rec_p.size(),
                       lattice.size(dir), n_columns);
#endif


    d_free(d_size);
    d_free(d_ptr);
    d_free(fft_wrk);
}


#endif // HILAPP

#endif