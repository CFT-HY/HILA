
#include "plumbing/defs.h"
#include "plumbing/lattice.h"
#include "plumbing/field.h"
#include "plumbing/backend_cuda/defs.h"

// hilapp needs to transform the include files above, to make them __device__
// callable...

#ifndef HILAPP

#if defined(CUDA)

#include <curand_kernel.h>

using gpurandState = curandState_t;
#define gpurand_init curand_init
#define gpurand_uniform curand_uniform
#define gpuMemcpyToSymbol(a, b, size, c, dir) \
    GPU_CHECK( cudaMemcpyToSymbol(a, b, size, c, dir) )
#define gpuGetDeviceCount(a) GPU_CHECK( cudaGetDeviceCount(a) )
#define gpuSetDevice(dev) GPU_CHECK( cudaSetDevice(dev) )
#define gpuGetLastError cudaGetLastError
#define gpuGetErrorString cudaGetErrorString

#elif defined(HIP)

#include <hip/hip_runtime.h>
#include <hiprand_kernel.h>

using gpurandState = hiprandState_t;
#define gpurand_init hiprand_init
#define gpurand_uniform hiprand_uniform
#define gpuMemcpyToSymbol(a, b, size, c, dir)                                          \
    GPU_CHECK( hipMemcpyToSymbol(HIP_SYMBOL(a), b, size, c, dir) )
#define gpuGetDeviceCount(a) GPU_CHECK( hipGetDeviceCount(a) )
#define gpuSetDevice(dev) GPU_CHECK( hipSetDevice(dev) )
#define gpuGetLastError hipGetLastError
#define gpuGetErrorString hipGetErrorString

#endif

// Save "constants" lattice size and volume here
__constant__ int64_t _d_volume;
__constant__ int _d_size[NDIM];
#ifndef EVEN_SITES_FIRST
__constant__ int _d_nodesize[NDIM];
__constant__ int _d_nodemin[NDIM];
__constant__ int _d_nodefactor[NDIM];
#endif

/* Random number generator */
gpurandState *gpurandstate;
__device__ gpurandState *d_gpurandstate;

/* Set seed on device */
__global__ void seed_random_kernel(gpurandState *state, unsigned long seed,
                                   unsigned int iters_per_kernel, unsigned int stride) {
    int x = threadIdx.x + blockIdx.x * blockDim.x;
    d_gpurandstate = state;
    for (int i = 0; i < iters_per_kernel; i++) {
        gpurand_init(seed + i * stride + x, i * stride + x, 0,
                     &d_gpurandstate[i * stride + x]);
    }
}

/* Set seed on device and host */
void hila::seed_device_rng(unsigned long seed) {
    unsigned int iters_per_kernel = 16;
    unsigned long n_blocks =
        lattice->mynode.volume() / (N_threads * iters_per_kernel) + 1;
    unsigned long n_sites = N_threads * n_blocks * iters_per_kernel;
    unsigned long myseed = seed + hila::myrank() * n_sites;
    gpuMalloc((void **)&gpurandstate, n_sites * sizeof(gpurandState));
#ifdef CUDA
    seed_random_kernel<<<n_blocks, N_threads>>>(gpurandstate, myseed, iters_per_kernel,
                                                n_blocks * N_threads);
#else
    hipLaunchKernelGGL(seed_random_kernel, dim3(n_blocks), dim3(N_threads), 0, 0,
                       gpurandstate, myseed, iters_per_kernel, n_blocks * N_threads);
#endif
    check_device_error("seed_random kernel");
}

/* Generate random numbers on device or host */
__device__ __host__ double hila::random() {
#ifdef __GPU_DEVICE_COMPILE__
    int x = threadIdx.x + blockIdx.x * blockDim.x;
    return gpurand_uniform(&d_gpurandstate[x]);
#else
    return mersenne();
#endif
}

// Then, define global functions loop_lattice_size() and _volume()
__device__ __host__ int loop_lattice_size(Direction dir) {
#ifdef __GPU_DEVICE_COMPILE__
    return _d_size[dir];
#else
    return lattice->size(dir);
#endif
}
__device__ __host__ CoordinateVector loop_lattice_size(void) {
#ifdef __GPU_DEVICE_COMPILE__
    CoordinateVector v;
    foralldir (d)
        v[d] = _d_size[d];
    return v;
#else
    return lattice->size();
#endif
}
__device__ __host__ int64_t loop_lattice_volume(void) {
#ifdef __GPU_DEVICE_COMPILE__
    return _d_volume;
#else
    return lattice->volume();
#endif
}

#ifndef EVEN_SITES_FIRST

__device__ const CoordinateVector backend_lattice_struct::coordinates(unsigned idx) const {
    CoordinateVector c;
    unsigned vdiv,ndiv;

    vdiv = idx;
    for (int d = 0; d < NDIM-1; ++d) {
        ndiv = vdiv / _d_nodesize[d];
        c[d] = vdiv - ndiv * _d_nodesize[d] + _d_nodemin[d];
        vdiv = ndiv;
    }
    c[NDIM-1] = vdiv + _d_nodemin[NDIM-1];

    return c;
}

__device__ int backend_lattice_struct::coordinate(unsigned idx, Direction dir) const {
    return ( idx / _d_nodefactor[dir] ) % _d_nodesize[dir] + _d_nodemin[dir];
}

#endif


void backend_lattice_struct::setup(lattice_struct *lattice) {
    CoordinateVector *tmp;

    /* Setup neighbour fields in all directions */
    for (int d = 0; d < NDIRS; d++) {
        // For normal boundaries
        gpuMalloc((void **)&(d_neighb[d]), lattice->mynode.volume() * sizeof(unsigned));
        gpuMemcpy(d_neighb[d], lattice->neighb[d],
                  lattice->mynode.volume() * sizeof(unsigned), gpuMemcpyHostToDevice);

#ifdef SPECIAL_BOUNDARY_CONDITIONS
        // For special boundaries
        // TODO: check this really works now!
        const unsigned *special_neighb =
            lattice->get_neighbour_array((Direction)d, BoundaryCondition::ANTIPERIODIC);

        if (special_neighb != lattice->neighb[d]) {
            gpuMalloc((void **)&(d_neighb_special[d]),
                      lattice->mynode.volume() * sizeof(unsigned));
            gpuMemcpy(d_neighb_special[d], special_neighb,
                      lattice->mynode.volume() * sizeof(unsigned), gpuMemcpyHostToDevice);
        } else {
            d_neighb_special[d] = d_neighb[d];
        }
#endif
    }

#ifdef EVEN_SITES_FIRST
    /* Setup the location field */
    gpuMalloc((void **)&(d_coordinates),
              lattice->mynode.volume() * sizeof(CoordinateVector));
    tmp = (CoordinateVector *)memalloc(lattice->mynode.volume() *
                                       sizeof(CoordinateVector));
    for (int i = 0; i < lattice->mynode.volume(); i++) tmp[i] = lattice->coordinates(i);

    gpuMemcpy(d_coordinates, tmp, lattice->mynode.volume() * sizeof(CoordinateVector),
              gpuMemcpyHostToDevice);
    free(tmp);
#endif

    // Other backend_lattice parameters
    field_alloc_size = lattice->field_alloc_size();

    int64_t v = lattice->volume();
    gpuMemcpyToSymbol(_d_volume, &v, sizeof(int64_t), 0, gpuMemcpyHostToDevice);
    int s[NDIM];
    foralldir (d)
        s[d] = lattice->size(d);
    gpuMemcpyToSymbol(_d_size, s, sizeof(int) * NDIM, 0, gpuMemcpyHostToDevice);

#ifndef EVEN_SITES_FIRST
    foralldir(d) s[d] = lattice->mynode.size[d];
    gpuMemcpyToSymbol(_d_nodesize, s, sizeof(int) * NDIM, 0, gpuMemcpyHostToDevice);

    foralldir(d) s[d] = lattice->mynode.min[d];
    gpuMemcpyToSymbol(_d_nodemin, s, sizeof(int) * NDIM, 0, gpuMemcpyHostToDevice);

    foralldir(d) s[d] = lattice->mynode.size_factor[d];
    gpuMemcpyToSymbol(_d_nodefactor, s, sizeof(int) * NDIM, 0, gpuMemcpyHostToDevice);

#endif    
}

void initialize_cuda(int rank) {
    int n_devices, my_device;

    gpuGetDeviceCount(&n_devices);
    check_device_error("Could not get device count");
    // This assumes that each node has the same number of mpi ranks and GPUs
    my_device = rank % n_devices;

    output0 << "CUDA/HIP Devices accessible from node 0: " << n_devices << '\n';
    // TODO: this only for node 0?
    if (rank < 6) {
        hila::output << "Cuda: rank " << rank << " choosing device " << my_device
                     << '\n';
    }
    if (numnodes() > 6) {
        output0 << "  + " << numnodes() - 6 << " more nodes\n";
    }

    gpuSetDevice(my_device);

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
        cudaDriverGetVersion(&driverVersion);
        cudaRuntimeGetVersion(&rtVersion);
        hila::output << "CUDA driver version: " << driverVersion << ", runtime "
                     << rtVersion << '\n';
        hila::output << "CUDART_VERSION " << CUDART_VERSION << '\n';
#if defined(CUDA_MALLOC_ASYNC)
        if (CUDART_VERSION >= 11020) {
            hila::output << "Using cudaMallocAsync() to allocate memory\n";
        }
#endif

        cudaDeviceProp props;
        int my_device;
        cudaGetDevice(&my_device);
        cudaGetDeviceProperties(&props, my_device);
        hila::output << "Device on node rank 0 device " << my_device << ":\n";
        hila::output << "  " << props.name << "  capability: " << props.major << "."
                     << props.minor << '\n';
        hila::output << "  Global memory:   " << props.totalGlobalMem / mb << "MB"
                     << '\n';
        hila::output << "  Shared memory:   " << props.sharedMemPerBlock / kb << "kB"
                     << '\n';
        hila::output << "  Constant memory: " << props.totalConstMem / kb << "kB"
                     << '\n';
        hila::output << "  Block registers: " << props.regsPerBlock << '\n';

        hila::output << "  Warp size:         " << props.warpSize << '\n';
        hila::output << "  Threads per block: " << props.maxThreadsPerBlock << '\n';
        hila::output << "  Max block dimensions: [ " << props.maxThreadsDim[0] << ", "
                     << props.maxThreadsDim[1] << ", " << props.maxThreadsDim[2] << " ]"
                     << '\n';
        hila::output << "  Max grid dimensions:  [ " << props.maxGridSize[0] << ", "
                     << props.maxGridSize[1] << ", " << props.maxGridSize[2] << " ]"
                     << '\n';

        hila::output << "Threads in use: " << N_threads << '\n';

// Following should be OK in open MPI
#ifdef OPEN_MPI
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
        hila::output << "OpenMPI library supports CUDA-Aware MPI\n";
        if (MPIX_Query_cuda_support() == 1)
            hila::output << "  Runtime library supports CUDA-Aware MPI\n";
        else {
            hila::output << "  Runtime library does not support CUDA-Aware MPI!\n";
#if defined(CUDA_AWARE_MPI)
            hila::output << "CUDA_AWARE_MPI is defined -- THIS MAY CRASH IN MPI\n";
#endif
        }
#else
        hila::output << "OpenMPI library does not support CUDA-Aware MPI\n";
#if defined(CUDA_AWARE_MPI)
        hila::output << "CUDA_AWARE_MPI is defined -- THIS MAY CRASH IN MPI\n";
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
        hipDriverGetVersion(&driverVersion);
        hipRuntimeGetVersion(&rtVersion);
        hila::output << "HIP driver version: " << driverVersion << ", runtime "
                     << rtVersion << '\n';

        hipDeviceProp_t props;
        int my_device;
        hipGetDevice(&my_device);
        hipGetDeviceProperties(&props, my_device);
        hila::output << "Device on node rank 0 device " << my_device << ":\n";
        hila::output << "  " << props.name << "  capability: " << props.major << "."
                     << props.minor << '\n';
        hila::output << "  Global memory:   " << props.totalGlobalMem / mb << "MB"
                     << '\n';
        hila::output << "  Shared memory:   " << props.sharedMemPerBlock / kb << "kB"
                     << '\n';
        hila::output << "  Constant memory: " << props.totalConstMem / kb << "kB"
                     << '\n';
        hila::output << "  Block registers: " << props.regsPerBlock << '\n';

        hila::output << "  Warp size:         " << props.warpSize << '\n';
        hila::output << "  Threads per block: " << props.maxThreadsPerBlock << '\n';
        hila::output << "  Max block dimensions: [ " << props.maxThreadsDim[0] << ", "
                     << props.maxThreadsDim[1] << ", " << props.maxThreadsDim[2] << " ]"
                     << '\n';
        hila::output << "  Max grid dimensions:  [ " << props.maxGridSize[0] << ", "
                     << props.maxGridSize[1] << ", " << props.maxGridSize[2] << " ]"
                     << '\n';
        hila::output << "Threads in use: " << N_threads << '\n';
    }
}

#endif

void gpu_exit_on_error(const char *msg, const char *file, int line) {
    gpuError code = gpuGetLastError();
    if (gpuSuccess != code) {
        hila::output << "CUDA error: " << msg << " in file " << file << " line " << line
                     << '\n';
        hila::output << "CUDA error string: " << gpuGetErrorString(code) << "\n";

        hila::terminate(0);
    }
}

void gpu_exit_on_error(gpuError code, const char *msg, const char *file, int line) {
    if (gpuSuccess != code) {
        hila::output << "CUDA error in command: " << msg << " in file " << file
                     << " line " << line << '\n';
        hila::output << "CUDA error string: " << gpuGetErrorString(code) << "\n";

        hila::terminate(0);
    }
}

#endif
