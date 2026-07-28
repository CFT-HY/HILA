
// Define below to deactivate "extern" in global var defs
#define IN_HILA_GPU

#include "plumbing/defs.h"
#include "plumbing/lattice.h"
#include "plumbing/field.h"
#include "plumbing/backend_gpu/defs.h"

// hilapp needs to transform the include files above, to make them __device__
// callable...

#ifndef HILAPP

#if defined(CUDA)

#include <curand_kernel.h>

//using gpurandState = curandState_t;
using gpurandState = curandStateXORWOW_t;
#define gpurand_init curand_init
#define gpurand_uint32 curand
#define gpurand_uniform curand_uniform
#define gpuGetDeviceCount(a) GPU_CHECK(cudaGetDeviceCount(a))
#define gpuSetDevice(dev) GPU_CHECK(cudaSetDevice(dev))
#define gpuGetLastError cudaGetLastError
#define gpuGetErrorString cudaGetErrorString

#elif defined(HIP)

#include <hip/hip_runtime.h>
#include <hiprand/hiprand_kernel.h>

//using gpurandState = hiprandState_t;
using gpurandState = hiprandStateXORWOW_t;
#define gpurand_init hiprand_init
#define gpurand_uint32 hiprand
#define gpurand_uniform hiprand_uniform
#define gpuGetDeviceCount(a) GPU_CHECK(hipGetDeviceCount(a))
#define gpuSetDevice(dev) GPU_CHECK(hipSetDevice(dev))
#define gpuGetLastError hipGetLastError
#define gpuGetErrorString hipGetErrorString

#endif

#if defined(CUDA) && defined(GPU_OVERLAP_COMM)
// EXPERIMENT (exp-overlap-aware-sync): runtime-selectable overlap stream strategy.
//   HILA_OVERLAP_STREAM_MODE = default | priority | green   (unset -> priority)
//     default  : plain non-blocking streams (original behaviour; the comparison baseline)
//     priority : halo_stream at greatest priority, so the comm kernel (NCCL SendRecv /
//                NVSHMEM proxy) is scheduled onto SMs as compute waves retire
//     green    : split the SMs into two green contexts -- comm on HILA_GREEN_SM SMs
//                (default 32), compute on the rest -- so they run on DISJOINT SMs.
//                Green requires the CUDA driver API and is compiled in only when
//                GPU_GREEN_CTX is defined (pulls in -lcuda); otherwise a green
//                request falls back to priority with a warning.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#ifdef GPU_GREEN_CTX
#include <cuda.h> // CUDA driver API: green contexts (SM partitioning)
#endif
namespace {
enum class OverlapStreamMode { Default, Priority, Green };
OverlapStreamMode overlap_stream_mode() {
    static OverlapStreamMode mode = [] {
        OverlapStreamMode m = OverlapStreamMode::Priority; // default when env unset
        const char *e = getenv("HILA_OVERLAP_STREAM_MODE");
        if (e) {
            if (strcmp(e, "default") == 0)
                m = OverlapStreamMode::Default;
            else if (strcmp(e, "priority") == 0)
                m = OverlapStreamMode::Priority;
            else if (strcmp(e, "green") == 0) {
#ifdef GPU_GREEN_CTX
                m = OverlapStreamMode::Green;
#else
                fprintf(stderr, "HILA_OVERLAP_STREAM_MODE=green requested but not "
                                "compiled in (rebuild with GPU_GREEN_CTX) -> priority\n");
                m = OverlapStreamMode::Priority;
#endif
            } else
                fprintf(stderr, "HILA_OVERLAP_STREAM_MODE: unknown '%s' -> priority\n", e);
        }
        return m;
    }();
    return mode;
}
const char *overlap_stream_mode_name(OverlapStreamMode m) {
    return m == OverlapStreamMode::Default ? "default"
           : m == OverlapStreamMode::Green ? "green"
                                           : "priority";
}
#ifdef GPU_GREEN_CTX
// Requested green-context comm-SM count (HILA_GREEN_SM, default 32). The realized
// split after granularity rounding is reported by green_streams() on first use.
int green_sm_request() {
    int comm_sm = 32;
    if (const char *e = getenv("HILA_GREEN_SM")) {
        int v = atoi(e);
        if (v > 0)
            comm_sm = v;
    }
    return comm_sm;
}

struct GreenStreams {
    gpuStream_t compute;
    gpuStream_t halo;
};
// Single fixed SM split (comm on HILA_GREEN_SM SMs, compute on the rest), built once.
// NOTE: an adaptive per-kernel load balancer over a pool of splits was prototyped and
// dropped -- see agent/HILA.md "Green-context overlap TODO". This is the user-pinned
// split only.
GreenStreams &green_streams() {
    static GreenStreams gs = [] {
        GreenStreams r{};
        int comm_sm = green_sm_request(); // HILA_GREEN_SM overrides (default 32)
#define HILA_DRV(x)                                                                                \
    do {                                                                                           \
        CUresult _e = (x);                                                                         \
        if (_e != CUDA_SUCCESS) {                                                                  \
            const char *_s = nullptr;                                                              \
            cuGetErrorString(_e, &_s);                                                             \
            fprintf(stderr, "green-ctx: %s -> %s\n", #x, _s ? _s : "?");                           \
            abort();                                                                               \
        }                                                                                          \
    } while (0)
        HILA_DRV(cuInit(0));
        int devId;
        GPU_CHECK(cudaGetDevice(&devId));
        CUdevice dev;
        HILA_DRV(cuDeviceGet(&dev, devId));
        CUdevResource sm;
        HILA_DRV(cuDeviceGetDevResource(dev, &sm, CU_DEV_RESOURCE_TYPE_SM));
        CUdevResource comm_res, compute_res;
        unsigned int nb = 1;
        HILA_DRV(cuDevSmResourceSplitByCount(&comm_res, &nb, &sm, &compute_res, 0, comm_sm));
        CUdevResourceDesc comm_desc, compute_desc;
        HILA_DRV(cuDevResourceGenerateDesc(&comm_desc, &comm_res, 1));
        HILA_DRV(cuDevResourceGenerateDesc(&compute_desc, &compute_res, 1));
        CUgreenCtx comm_g, compute_g;
        HILA_DRV(cuGreenCtxCreate(&comm_g, comm_desc, dev, CU_GREEN_CTX_DEFAULT_STREAM));
        HILA_DRV(cuGreenCtxCreate(&compute_g, compute_desc, dev, CU_GREEN_CTX_DEFAULT_STREAM));
        int lo, hi;
        cudaDeviceGetStreamPriorityRange(&lo, &hi);
        CUstream cs, hs;
        HILA_DRV(cuGreenCtxStreamCreate(&cs, compute_g, CU_STREAM_NON_BLOCKING, 0));
        HILA_DRV(cuGreenCtxStreamCreate(&hs, comm_g, CU_STREAM_NON_BLOCKING, hi));
        r.compute = (gpuStream_t)cs;
        r.halo = (gpuStream_t)hs;
        if (hila::myrank() == 0)
            fprintf(stderr, "green-ctx: comm=%u SMs  compute=%u SMs\n", comm_res.sm.smCount,
                    compute_res.sm.smCount);
#undef HILA_DRV
        return r;
    }();
    return gs;
}
#endif // GPU_GREEN_CTX
} // namespace
#endif

void hila::report_overlap_config(std::ostream &out) {
#if defined(GPU_OVERLAP_COMM)
#if defined(CUDA)
    OverlapStreamMode m = overlap_stream_mode();
    out << "Overlap stream mode   : " << overlap_stream_mode_name(m)
#ifdef GPU_GREEN_CTX
        << "  (HILA_OVERLAP_STREAM_MODE=default|priority|green)\n";
    if (m == OverlapStreamMode::Green)
        out << "Green-ctx comm SMs    : " << green_sm_request()
            << " requested  (HILA_GREEN_SM; realized split logged at first use)\n";
#else
        << "  (HILA_OVERLAP_STREAM_MODE=default|priority; green needs GPU_GREEN_CTX)\n";
#endif
#elif defined(HIP)
    out << "Overlap stream mode   : priority (fixed; high-priority halo_stream)\n";
#endif
#endif
}

gpuStreamPool &hila::stream_pool() {
    static gpuStreamPool instance;
    return instance;
}

gpuStream_t &hila::halo_stream() {
#if defined(CUDA) && defined(GPU_OVERLAP_COMM)
    // runtime lever: green -> dedicated green-context stream; priority -> greatest
    // priority; default -> plain non-blocking (the comparison baseline)
#ifdef GPU_GREEN_CTX
    if (overlap_stream_mode() == OverlapStreamMode::Green)
        return green_streams().halo;
#endif
    static gpuStream_t instance = [] {
        gpuStream_t stream;
        if (overlap_stream_mode() == OverlapStreamMode::Priority) {
            int lo, hi;
            cudaDeviceGetStreamPriorityRange(&lo, &hi);
            GPU_CHECK(cudaStreamCreateWithPriority(&stream, cudaStreamNonBlocking, hi));
        } else { // Default
            GPU_CHECK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
        }
        return stream;
    }();
    return instance;
#elif defined(HIP) && defined(GPU_OVERLAP_COMM)
    static gpuStream_t instance = [] {
        gpuStream_t stream;
        int lo, hi;
        hipDeviceGetStreamPriorityRange(&lo, &hi);
        GPU_CHECK(hipStreamCreateWithPriority(&stream, hipStreamNonBlocking, hi));
        return stream;
    }();
    return instance;
#else
    static gpuStream_t instance = [] {
        gpuStream_t stream;
        gpuStreamCreateWithFlags(&stream, gpuStreamNonBlocking);
        return stream;
    }();
    return instance;
#endif
}

gpuEvent_t &hila::halo_event() {
    static gpuEvent_t instance = [] {
        gpuEvent_t e;
        gpuEventCreate(&e);
        return e;
    }();
    return instance;
}

gpuStream_t &hila::compute_stream() {
#if defined(GPU_OVERLAP_COMM)
#if defined(CUDA) && defined(GPU_GREEN_CTX)
    // green mode runs the bulk compute on the (larger) compute green context;
    // default/priority use a normal full-GPU non-blocking stream
    if (overlap_stream_mode() == OverlapStreamMode::Green)
        return green_streams().compute;
#endif
    static gpuStream_t instance = [] {
        gpuStream_t stream;
        gpuStreamCreateWithFlags(&stream, gpuStreamNonBlocking);
        return stream;
    }();
    return instance;
#else
    static gpuStream_t instance = 0;
    return instance;
#endif
}

gpuEvent_t &hila::compute_event() {
    static gpuEvent_t instance = [] {
        gpuEvent_t e;
        gpuEventCreate(&e);
        return e;
    }();
    return instance;
}


// // Save "constants" lattice size and volume here
// __constant__ int64_t _d_volume;
// // __constant__ int _d_size[NDIM];
// __constant__ CoordinateVector _d_size;
// #ifndef EVEN_SITES_FIRST
// __constant__ int _d_nodesize[NDIM];
// __constant__ int _d_nodemin[NDIM];
// __constant__ int _d_nodefactor[NDIM];
// #endif

/* Random number generator */
static gpurandState *gpurandstateptr = nullptr;
__constant__ gpurandState *d_gpurandstateptr;

// check if rng on device is OK

bool hila::is_device_rng_on() {
    return gpurandstateptr != nullptr;
}

/* Set seed on device */
__global__ void seed_random_kernel(unsigned long long seed) {
    unsigned x = threadIdx.x + blockIdx.x * blockDim.x;
    //  d_gpurandstateptr set now using memcpyToSymbol
    //  d_gpurandstateptr = state;
    gpurand_init(seed + x, x, 0, &d_gpurandstateptr[x]);
}

/* Set seed on device and host */
void hila::initialize_device_rng(uint64_t seed) {

#if defined(GPU_RNG_THREAD_BLOCKS) && GPU_RNG_THREAD_BLOCKS < 0
    hila::out0 << "GPU RANDOM NUMBERS DISABLED, GPU_RNG_THREAD_BLOCKS < 0\n";
    return;

#else

    if (is_device_rng_on()) {
        hila::out0 << "Reseeding GPU random numbers, new seed " << seed << '\n';
    }

    unsigned long n_blocks = (lattice->mynode.volume + N_threads - 1) / N_threads;

#if defined(GPU_RNG_THREAD_BLOCKS) && GPU_RNG_THREAD_BLOCKS > 0
    // If we have limited rng block number
    if (GPU_RNG_THREAD_BLOCKS < n_blocks) {
        n_blocks = GPU_RNG_THREAD_BLOCKS;
    }

    if (!is_device_rng_on()) {
        hila::out0 << "Initializing GPU random number generator\n";
        hila::out0 << "GPU random number thread blocks: " << n_blocks << " of size " << N_threads
                   << " threads\n";
    }
#else
    if (!is_device_rng_on()) {
        hila::out0 << "Initializing GPU random number generator\n";
        hila::out0 << "GPU random numbers: using one generator/site (GPU_RNG_THREAD_BLOCKS = 0 or "
                      "undefined)\n";
    }
#endif

    unsigned long long n_sites = n_blocks * N_threads;
    unsigned long long myseed = seed + hila::myrank() * n_sites;

    if (!is_device_rng_on()) {
        // allocate random state and copy the ptr to d_gpurandstateptr
        gpuMalloc(&gpurandstateptr, n_sites * sizeof(gpurandState));
        gpuMemcpyToSymbol(d_gpurandstateptr, &gpurandstateptr, sizeof(gpurandState *), 0,
                          gpuMemcpyHostToDevice);
    }

#ifdef CUDA
    seed_random_kernel<<<n_blocks, N_threads>>>(myseed);
#else
    hipLaunchKernelGGL(seed_random_kernel, dim3(n_blocks), dim3(N_threads), 0, 0, myseed);
#endif
    check_device_error("seed_random kernel");

#endif

}

void hila::free_device_rng() {
    if (is_device_rng_on()) {
        gpuFree(gpurandstateptr);
        gpurandstateptr = nullptr;
        // set d_gpurandstateptr <- nullptr.
        gpuMemcpyToSymbol(d_gpurandstateptr, &gpurandstateptr, sizeof(gpurandState *), 0,
                          gpuMemcpyHostToDevice);

        // good to purge the memory pool after releasing a large chunk
        gpu_memory_pool_purge();
    }
}

/* Generate random numbers on device or host. Values are in range [0,1) */
__device__ __host__ double hila::random() {
#ifdef _GPU_DEVICE_COMPILE_
    unsigned x = threadIdx.x + blockIdx.x * blockDim.x;

    // ---
    // Method 1: convert 32-bit unsigned int directly to double,
    // i.e. return  uint32 * 2^(-32)
    // This number below is 2^(-32) in double
    // THERE SHOULD BE NO REASON TO USE THIS
    //
    // return gpurand_uint32(&d_gpurandstateptr[x]) * 2.3283064365386963e-10;

    // ---

    // Method 2: use 2 32-bit unsigned random ints to generate a 53-bit random value,
    // which is then converted to double by multiplying by 2^(-53)
    // This should be the preferred method
    unsigned u1 = gpurand_uint32(&d_gpurandstateptr[x]);
    unsigned u2 = gpurand_uint32(&d_gpurandstateptr[x]);
    unsigned long long z = ((unsigned long long)u1) ^ ((unsigned long long)u2 << (53 - 32));
    return z * 1.1102230246251565e-16; // = 2^(-53)

    // ---

    // This below is std cuda/hip uniform rng
    // DO NOT USE, GIVES WRONG RANGE (0,1], AND ACTUALLY RETURNS FLOATS
    // return gpurand_uniform(&d_gpurandstateptr[x]);

#else
    return hila::host_random();
#endif
}


///////////////////////////////////////////////////////////////////////////////////////
// Setup the lattice struct on GPUs:
// allocate neighbour and coordinate arrays
// setup global variables in __constant__ memory

void backend_lattice_struct::setup(lattice_struct &lat) {
    CoordinateVector *tmp;

    /* Setup neighbour fields in all directions */
    for (int d = 0; d < NDIRS; d++) {
        // For normal boundaries
        gpuMalloc(&(d_neighb[d]), lat.mynode.volume * sizeof(unsigned));
        gpuMemcpy(d_neighb[d], lat.neighb[d], lat.mynode.volume * sizeof(unsigned),
                  gpuMemcpyHostToDevice);

#ifdef SPECIAL_BOUNDARY_CONDITIONS
        // For special boundaries
        // TODO: check this really works now!
        const unsigned *special_neighb =
            lat.get_neighbour_array((Direction)d, hila::bc::ANTIPERIODIC);

        if (special_neighb != lat.neighb[d]) {
            gpuMalloc(&(d_neighb_special[d]), lat.mynode.volume * sizeof(unsigned));
            gpuMemcpy(d_neighb_special[d], special_neighb, lat.mynode.volume * sizeof(unsigned),
                      gpuMemcpyHostToDevice);
        } else {
            d_neighb_special[d] = d_neighb[d];
        }
#endif
    }

#ifdef EVEN_SITES_FIRST
    /* Setup the location field */
    gpuMalloc(&(d_coordinates), lat.mynode.volume * sizeof(CoordinateVector));
    tmp = (CoordinateVector *)memalloc(lat.mynode.volume * sizeof(CoordinateVector));
    for (unsigned i = 0; i < lat.mynode.volume; i++)
        tmp[i] = lat.coordinates(i);

    gpuMemcpy(d_coordinates, tmp, lat.mynode.volume * sizeof(CoordinateVector),
              gpuMemcpyHostToDevice);
    free(tmp);


#endif

    // Other backend_lattice parameters
    field_alloc_size = lat.mynode.field_alloc_size;

    set_device_globals(lat);
}

#endif // not HILAPP

    // set some gobal variables, visible on GPUs
    // thus, hilapp needs to see this definition

    void backend_lattice_struct::set_device_globals(const lattice_struct &lat) {

#ifndef HILAPP
#if defined(GPU_CCL) || defined(GPU_SHMEM)
        {
            int n_devices;
            gpuGetDeviceCount(&n_devices);
            check_device_error("Could not get device count");
            gpuSetDevice(lat.mynode.rank % n_devices);
        }
#endif
#endif

#ifdef EVEN_SITES_FIRST

        gpuMemcpyToSymbol(_dev_coordinates, &d_coordinates, sizeof(CoordinateVector *), 0,
                          gpuMemcpyHostToDevice);
#endif

        gpuMemcpyToSymbol(_dev_field_alloc_size, &field_alloc_size, sizeof(unsigned), 0,
                          gpuMemcpyHostToDevice);

        _d_volume = lat.l_volume;
        _d_size = lat.l_size;

#ifndef EVEN_SITES_FIRST

        _d_nodesize = lat.mynode.size;
        _d_nodemin = lat.mynode.min;
        _d_nodefactor = lat.mynode.size_factor;

        // foralldir(d) s[d] = lat.mynode.size[d];
        // gpuMemcpyToSymbol(_d_nodesize, s, sizeof(int) * NDIM, 0, gpuMemcpyHostToDevice);

        // foralldir(d) s[d] = lat.mynode.min[d];
        // gpuMemcpyToSymbol(_d_nodemin, s, sizeof(int) * NDIM, 0, gpuMemcpyHostToDevice);

        // foralldir(d) s[d] = lat.mynode.size_factor[d];
        // gpuMemcpyToSymbol(_d_nodefactor, s, sizeof(int) * NDIM, 0, gpuMemcpyHostToDevice);

#endif
    }

#ifndef HILAPP
    // again, hilapp can skip this part

    void initialize_gpu(int rank, int device) {
        int n_devices, my_device;

        gpuGetDeviceCount(&n_devices);
        check_device_error("Could not get device count");
        // This assumes that each node has the same number of mpi ranks and GPUs
        // TODO:generalize (if needed)
        if (device > 0 && hila::number_of_nodes() == 1) {
            if (device >= n_devices) {
                hila::out0 << "-device " << device << ": too large device number, maximum "
                           << n_devices - 1 << " on this machine\n";
                hila::terminate(0);
            }

            my_device = device;
        } else {
            my_device = rank % n_devices;
        }


        hila::out0 << "GPU devices accessible from node 0: " << n_devices << '\n';

        // TODO: this only for node 0?
        if (n_devices > 1 && rank < 6) {
            hila::out << "GPU: MPI rank " << rank << " choosing device " << my_device << std::endl;
            if (hila::number_of_nodes() > 6) {
                hila::out0 << "  + " << hila::number_of_nodes() - 6 << " more nodes\n";
            }
        }

        gpuSetDevice(my_device);

        // set gpu rng state to "off", to prevent accidental use
        gpurandstateptr = nullptr;
        // set d_gpurandstateptr <- nullptr.
        gpuMemcpyToSymbol(d_gpurandstateptr, &gpurandstateptr, sizeof(gpurandState *), 0,
                          gpuMemcpyHostToDevice);


#if defined(CUDA_MALLOC_ASYNC)
        // set memory pool
        cudaMemPool_t mempool;
        cudaDeviceGetDefaultMemPool(&mempool, my_device);
        uint64_t threshold = UINT64_MAX;
        cudaMemPoolSetAttribute(mempool, cudaMemPoolAttrReleaseThreshold, &threshold);

#endif
    }

    // if using NCCL or RCCL or NVSHMEM
    namespace hila {

#ifdef GPU_CCL
    /**
     * @brief initialize nccl/rccl communicator for GPU communication
     * @details Uses same mapping as MPI communicator (gpu per task)
     *
     */
    void initialize_gccl_communication() {
        int rank = lattice->mynode.rank;
        int size = lattice->nodes.number;
        std::cout << "pre set device " << "rank: " << rank << " num ranks " << size << std::endl;

        {
            int n_devices;
            gpuGetDeviceCount(&n_devices);
            check_device_error("Could not get device count");
            gpuSetDevice(rank % n_devices);
        }
        // gpuSetDevice(rank);
        std::cout << "Post set device " << "rank: " << rank << " num ranks " << size << std::endl;
        gcclComm_t communicator;
        std::cout << "Post gcclComm_t constructor " << "rank: " << rank << " num ranks " << size
                  << std::endl;

        gcclUniqueId unique_id;
        std::cout << "Post gcclUniqueId constructor " << "rank: " << rank << " num ranks " << size
                  << std::endl;

        if (rank == 0) {
            gcclGetUniqueId(&unique_id);
        }
        std::cout << "Post get unique id " << "rank: " << rank << " num ranks " << size
                  << std::endl;

        MPI_Bcast(&unique_id, sizeof(unique_id), MPI_BYTE, 0, lattice->mpi_comm_lat);
        std::cout << "Post Bcast" << "rank: " << rank << " num ranks " << size << std::endl;

        gcclCommInitRank(&communicator, size, unique_id, rank);
        std::cout << "Post Init comm rank" << "rank: " << rank << " num ranks " << size
                  << std::endl;

        lattice.ptr()->gccl_comm_lat = communicator;
        // double *broadcast_val;
        // double *recieve;
        // double recieve_host;
        // gpuMalloc(&recieve, sizeof(double) * 2048 * 200);
        // gpuMalloc(&broadcast_val, sizeof(double) * 2048 * 200);

        // double val = 2.718281828459045;
        // gpuMemcpy(broadcast_val, &val, sizeof(double), gpuMemcpyHostToDevice);

        // gcclGroupStart();
        // gcclAllReduce(broadcast_val, recieve, 2048 * 200, gccl_type<double>::value, ncclSum,
        //               communicator, hila::compute_stream());
        //// gcclSend(broadcast_val, 2048*200, gccl_type<double>::value, (rank+size/2)%size,
        /// communicator, / hila::compute_stream()); gcclRecv(recieve, 2048*200,
        /// gccl_type<double>::value,
        //// (rank-size/2+size)%size, communicator, hila::compute_stream());
        // gcclGroupEnd();
        // gpuStreamSynchronize(hila::compute_stream());
        // gpuMemcpy(&recieve_host, recieve, sizeof(double), gpuMemcpyDeviceToHost);
        // std::cout << "Post Broadcast " << "rank: " << rank << " num " << recieve_host <<
        // std::endl;
    }
#endif // GPU_CCL
#ifdef GPU_SHMEM
#include <nvshmem.h>
#include <nvshmemx.h>
    void initialize_nvshmem_communication() {
        int rank, size;

        {
            int n_devices;
            gpuGetDeviceCount(&n_devices);
            check_device_error("Could not get device count");
            gpuSetDevice(lattice->mynode.rank % n_devices);
        }

        nvshmemx_init_attr_t attr;
        attr.mpi_comm = &lattice.ptr()->mpi_comm_lat;
        nvshmemx_init_attr(NVSHMEMX_INIT_WITH_MPI_COMM, &attr);
        rank = nvshmem_my_pe();
        size = nvshmem_n_pes();
        std::cout << "NVSHMEM initialized. Rank: " << rank << " Size: " << size << std::endl;
    }

    void finalize_nvshmem_communication() {}

#endif // GPU_SHMEM
    } // namespace hila


#ifdef CUDA

#ifdef OPEN_MPI
// here functions to inquire cuda-aware MPI defined
#include "mpi-ext.h"
#endif

    void gpu_device_info() {
        if_rank0 () {
            const int kb = 1024;
            const int mb = kb * kb;

            int driverVersion, rtVersion;
            GPU_CHECK(cudaDriverGetVersion(&driverVersion));
            GPU_CHECK(cudaRuntimeGetVersion(&rtVersion));
            hila::out << "CUDA driver version: " << driverVersion << ", runtime " << rtVersion
                      << '\n';
            hila::out << "CUDART_VERSION " << CUDART_VERSION << '\n';
#if defined(CUDA_MALLOC_ASYNC)
            if (CUDART_VERSION >= 11020) {
                hila::out << "Using cudaMallocAsync() to allocate memory\n";
            }
#endif

            cudaDeviceProp props;
            int my_device;
            GPU_CHECK(cudaGetDevice(&my_device));
            GPU_CHECK(cudaGetDeviceProperties(&props, my_device));
            hila::out << "Device on node rank 0 device " << my_device << ":\n";
            hila::out << "  " << props.name << "  capability: " << props.major << "." << props.minor
                      << '\n';
            hila::out << "  Global memory:   " << props.totalGlobalMem / mb << "MB" << '\n';
            hila::out << "  Shared memory:   " << props.sharedMemPerBlock / kb << "kB" << '\n';
            hila::out << "  Constant memory: " << props.totalConstMem / kb << "kB" << '\n';
            hila::out << "  Block registers: " << props.regsPerBlock << '\n';

            hila::out << "  Warp size:         " << props.warpSize << '\n';
            hila::out << "  Threads per block: " << props.maxThreadsPerBlock << '\n';
            hila::out << "  Max block dimensions: [ " << props.maxThreadsDim[0] << ", "
                      << props.maxThreadsDim[1] << ", " << props.maxThreadsDim[2] << " ]" << '\n';
            hila::out << "  Max grid dimensions:  [ " << props.maxGridSize[0] << ", "
                      << props.maxGridSize[1] << ", " << props.maxGridSize[2] << " ]" << '\n';

            hila::out << "Thread block size used: " << N_threads << '\n';


            hila::out << "WARNING: GPU_BLOCK_REDUCTION_THREADS (" << GPU_BLOCK_REDUCTION_THREADS
                      << ") may exceed available shared memory (" << props.sharedMemPerBlock
                      << " bytes). Consider reducing it.\n";


// Following should be OK in open MPI
#ifdef OPEN_MPI
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
            hila::out << "OpenMPI library supports CUDA-Aware MPI\n";
            if (MPIX_Query_cuda_support() == 1)
                hila::out << "  Runtime library supports CUDA-Aware MPI\n";
            else {
                hila::out << "  Runtime library does not support CUDA-Aware MPI!\n";
#if defined(GPU_AWARE_COMM)
                hila::out << "GPU_AWARE_COMM is defined -- THIS MAY CRASH IN MPI\n";
#endif
            }
#else
            hila::out << "OpenMPI library does not support CUDA-Aware MPI\n";
#if defined(GPU_AWARE_COMM)
            hila::out << "GPU_AWARE_COMM is defined -- THIS MAY CRASH IN MPI\n";
#endif
#endif // MPIX
#endif // OPEN_MPI
        }
    }
#endif

#ifdef HIP

    void gpu_device_info() {
        if_rank0 () {
            const int kb = 1024;
            const int mb = kb * kb;

            int driverVersion, rtVersion;
            GPU_CHECK(hipDriverGetVersion(&driverVersion));
            GPU_CHECK(hipRuntimeGetVersion(&rtVersion));
            hila::out << "HIP driver version: " << driverVersion << ", runtime " << rtVersion
                      << '\n';

            hipDeviceProp_t props;
            int my_device;
            GPU_CHECK(hipGetDevice(&my_device));
            GPU_CHECK(hipGetDeviceProperties(&props, my_device));
            hila::out << "Device on node rank 0 device " << my_device << ":\n";
            hila::out << "  " << props.name << "  capability: " << props.major << "." << props.minor
                      << '\n';
            hila::out << "  Global memory:   " << props.totalGlobalMem / mb << "MB" << '\n';
            hila::out << "  Shared memory:   " << props.sharedMemPerBlock / kb << "kB" << '\n';
            hila::out << "  Constant memory: " << props.totalConstMem / kb << "kB" << '\n';
            hila::out << "  Block registers: " << props.regsPerBlock << '\n';

            hila::out << "  Warp size:         " << props.warpSize << '\n';
            hila::out << "  Threads per block: " << props.maxThreadsPerBlock << '\n';
            hila::out << "  Max block dimensions: [ " << props.maxThreadsDim[0] << ", "
                      << props.maxThreadsDim[1] << ", " << props.maxThreadsDim[2] << " ]" << '\n';
            hila::out << "  Max grid dimensions:  [ " << props.maxGridSize[0] << ", "
                      << props.maxGridSize[1] << ", " << props.maxGridSize[2] << " ]" << '\n';
            hila::out << "Thread block size used: " << N_threads << '\n';

            hila::out << "WARNING: GPU_BLOCK_REDUCTION_THREADS (" << GPU_BLOCK_REDUCTION_THREADS
                      << ") may exceed available shared memory (" << props.sharedMemPerBlock
                      << " bytes). Consider reducing it.\n";
        }
    }

#endif

    void gpu_exit_on_error(const char *msg, const char *file, int line) {
        gpuError code = gpuGetLastError();
        if (gpuSuccess != code) {
            hila::out << GPUTYPESTR << " error: " << msg << " in file " << file << " line " << line
                      << '\n';
            hila::out << GPUTYPESTR << " error string: " << gpuGetErrorString(code) << "\n";

            hila::terminate(0);
        }
    }

    void gpu_exit_on_error(gpuError code, const char *msg, const char *file, int line) {
        if (gpuSuccess != code) {
            hila::out << GPUTYPESTR << " error in command: " << msg << " in file " << file
                      << " line " << line << '\n';
            hila::out << GPUTYPESTR << " error string: " << gpuGetErrorString(code) << "\n";

            hila::terminate(0);
        }
    }

#endif // not HILAPP
