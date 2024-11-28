
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

using gpurandState = curandState_t;
#define gpurand_init curand_init
#define gpurand_uniform curand_uniform
#define gpuGetDeviceCount(a) GPU_CHECK(cudaGetDeviceCount(a))
#define gpuSetDevice(dev) GPU_CHECK(cudaSetDevice(dev))
#define gpuGetLastError cudaGetLastError
#define gpuGetErrorString cudaGetErrorString

#elif defined(HIP)

#include <hip/hip_runtime.h>
#include <hiprand/hiprand_kernel.h>

using gpurandState = hiprandState_t;
#define gpurand_init hiprand_init
#define gpurand_uniform hiprand_uniform
#define gpuGetDeviceCount(a) GPU_CHECK(hipGetDeviceCount(a))
#define gpuSetDevice(dev) GPU_CHECK(hipSetDevice(dev))
#define gpuGetLastError hipGetLastError
#define gpuGetErrorString hipGetErrorString

#endif

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
static gpurandState *gpurandstateptr;
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
    gpurand_init(seed + x, 0, 0, &d_gpurandstateptr[x]);
}

/* Set seed on device and host */
void hila::initialize_device_rng(uint64_t seed) {
    unsigned long n_blocks = (lattice.mynode.volume() + N_threads - 1) / N_threads;

#if defined(GPU_RNG_THREAD_BLOCKS) && GPU_RNG_THREAD_BLOCKS > 0
    // If we have limited rng block number
    if (GPU_RNG_THREAD_BLOCKS < n_blocks) {
        n_blocks = GPU_RNG_THREAD_BLOCKS;
    }

    hila::out0 << "GPU random number generator initialized\n";
    hila::out0 << "GPU random number thread blocks: " << n_blocks << " of size " << N_threads
               << " threads\n";
#elif defined(GPU_RNG_THREAD_BLOCKS) && GPU_RNG_THREAD_BLOCKS < 0
    hila::out0 << "GPU RANDOM NUMBERS DISABLED, GPU_RNG_THREAD_BLOCKS < 0\n";
#else
    hila::out0 << "GPU random number generator initialized\n";
    hila::out0
        << "GPU random numbers: using on generator/site (GPU_RNG_THREAD_BLOCKS = 0 or undefined)\n";
#endif

    unsigned long long n_sites = n_blocks * N_threads;
    unsigned long long myseed = seed + hila::myrank() * n_sites;

    // allocate random state and copy the ptr to d_gpurandstateptr
    gpuMalloc(&gpurandstateptr, n_sites * sizeof(gpurandState));
    gpuMemcpyToSymbol(d_gpurandstateptr, &gpurandstateptr, sizeof(gpurandState *), 0,
                      gpuMemcpyHostToDevice);

#ifdef CUDA
    seed_random_kernel<<<n_blocks, N_threads>>>(myseed);
#else
    hipLaunchKernelGGL(seed_random_kernel, dim3(n_blocks), dim3(N_threads), 0, 0, myseed);
#endif
    check_device_error("seed_random kernel");
}

void hila::free_device_rng() {
    if (gpurandstateptr != nullptr) {
        gpuFree(gpurandstateptr);
        gpurandstateptr = nullptr;
        // set d_gpurandstateptr <- nullptr.
        gpuMemcpyToSymbol(d_gpurandstateptr, &gpurandstateptr, sizeof(gpurandState *), 0,
                          gpuMemcpyHostToDevice);

        // good to purge the memory pool after releasing a large chunk
        gpu_memory_pool_purge();
    }
}

/* Generate random numbers on device or host */
__device__ __host__ double hila::random() {
#ifdef _GPU_DEVICE_COMPILE_
    unsigned x = threadIdx.x + blockIdx.x * blockDim.x;
    return gpurand_uniform(&d_gpurandstateptr[x]);
#else
    return hila::host_random();
#endif
}


///////////////////////////////////////////////////////////////////////////////////////
// Setup the lattice struct on GPUs:
// allocate neighbour and coordinate arrays
// setup global variables in __constant__ memory

void backend_lattice_struct::setup(lattice_struct &lattice) {
    CoordinateVector *tmp;

    /* Setup neighbour fields in all directions */
    for (int d = 0; d < NDIRS; d++) {
        // For normal boundaries
        gpuMalloc(&(d_neighb[d]), lattice.mynode.volume() * sizeof(unsigned));
        gpuMemcpy(d_neighb[d], lattice.neighb[d], lattice.mynode.volume() * sizeof(unsigned),
                  gpuMemcpyHostToDevice);

#ifdef SPECIAL_BOUNDARY_CONDITIONS
        // For special boundaries
        // TODO: check this really works now!
        const unsigned *special_neighb =
            lattice.get_neighbour_array((Direction)d, hila::bc::ANTIPERIODIC);

        if (special_neighb != lattice.neighb[d]) {
            gpuMalloc(&(d_neighb_special[d]), lattice.mynode.volume() * sizeof(unsigned));
            gpuMemcpy(d_neighb_special[d], special_neighb,
                      lattice.mynode.volume() * sizeof(unsigned), gpuMemcpyHostToDevice);
        } else {
            d_neighb_special[d] = d_neighb[d];
        }
#endif
    }

#ifdef EVEN_SITES_FIRST
    /* Setup the location field */
    gpuMalloc(&(d_coordinates), lattice.mynode.volume() * sizeof(CoordinateVector));
    tmp = (CoordinateVector *)memalloc(lattice.mynode.volume() * sizeof(CoordinateVector));
    for (unsigned i = 0; i < lattice.mynode.volume(); i++)
        tmp[i] = lattice.coordinates(i);

    gpuMemcpy(d_coordinates, tmp, lattice.mynode.volume() * sizeof(CoordinateVector),
              gpuMemcpyHostToDevice);
    free(tmp);
#endif

    // Other backend_lattice parameters
    field_alloc_size = lattice.field_alloc_size();

    set_lattice_globals(lattice);

}

#endif // not HILAPP

// set some gobal variables, visible on GPUs
// thus, hilapp needs to see this definition

void backend_lattice_struct::set_lattice_globals(lattice_struct &lattice) {

    _d_volume = lattice.volume();
    _d_size = lattice.size();

#ifndef EVEN_SITES_FIRST

    _d_nodesize = lattice.mynode.size;
    _d_nodemin = lattice.mynode.min;
    _d_nodefactor = lattice.mynode.size_factor;

    // foralldir(d) s[d] = lattice.mynode.size[d];
    // gpuMemcpyToSymbol(_d_nodesize, s, sizeof(int) * NDIM, 0, gpuMemcpyHostToDevice);

    // foralldir(d) s[d] = lattice.mynode.min[d];
    // gpuMemcpyToSymbol(_d_nodemin, s, sizeof(int) * NDIM, 0, gpuMemcpyHostToDevice);

    // foralldir(d) s[d] = lattice.mynode.size_factor[d];
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

#ifdef CUDA

#ifdef OPEN_MPI
// here functions to inquire cuda-aware MPI defined
#include "mpi-ext.h"
#endif

void gpu_device_info() {
    if (hila::myrank() == 0) {
        const int kb = 1024;
        const int mb = kb * kb;

        int driverVersion, rtVersion;
        GPU_CHECK(cudaDriverGetVersion(&driverVersion));
        GPU_CHECK(cudaRuntimeGetVersion(&rtVersion));
        hila::out << "CUDA driver version: " << driverVersion << ", runtime " << rtVersion << '\n';
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

// Following should be OK in open MPI
#ifdef OPEN_MPI
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
        hila::out << "OpenMPI library supports CUDA-Aware MPI\n";
        if (MPIX_Query_cuda_support() == 1)
            hila::out << "  Runtime library supports CUDA-Aware MPI\n";
        else {
            hila::out << "  Runtime library does not support CUDA-Aware MPI!\n";
#if defined(GPU_AWARE_MPI)
            hila::out << "GPU_AWARE_MPI is defined -- THIS MAY CRASH IN MPI\n";
#endif
        }
#else
        hila::out << "OpenMPI library does not support CUDA-Aware MPI\n";
#if defined(GPU_AWARE_MPI)
        hila::out << "GPU_AWARE_MPI is defined -- THIS MAY CRASH IN MPI\n";
#endif
#endif // MPIX
#endif // OPEN_MPI
    }
}
#endif

#ifdef HIP

void gpu_device_info() {
    if (hila::myrank() == 0) {
        const int kb = 1024;
        const int mb = kb * kb;

        int driverVersion, rtVersion;
        GPU_CHECK(hipDriverGetVersion(&driverVersion));
        GPU_CHECK(hipRuntimeGetVersion(&rtVersion));
        hila::out << "HIP driver version: " << driverVersion << ", runtime " << rtVersion << '\n';

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
        hila::out << GPUTYPESTR << " error in command: " << msg << " in file " << file << " line "
                  << line << '\n';
        hila::out << GPUTYPESTR << " error string: " << gpuGetErrorString(code) << "\n";

        hila::terminate(0);
    }
}

#endif  // not HILAPP



