
#include "plumbing/defs.h"
#include "plumbing/lattice.h"
#include "plumbing/field.h"
#include "plumbing/backend_cuda/defs.h"

// hilapp needs to transform the include files above, to make them __device__ callable...

#ifndef HILAPP

// Save "constants" lattice size and volume here
__constant__ int _d_size[NDIM];
__constant__ int64_t _d_volume;

/* Random number generator */
curandState *curandstate;
__device__ curandState *d_curandstate;
#define cuda_rand_setup_count

/* Set seed on device */
__global__ void seed_random_kernel(curandState *state, unsigned long seed,
                                   unsigned int iters_per_kernel, unsigned int stride) {
    int x = threadIdx.x + blockIdx.x * blockDim.x;
    d_curandstate = state;
    for (int i = 0; i < iters_per_kernel; i++) {
        curand_init(seed + i * stride + x, i * stride + x, 0,
                    &d_curandstate[i * stride + x]);
    }
}

/* Set seed on device and host */
void hila::seed_device_rng(unsigned long seed) {
    unsigned int iters_per_kernel = 16;
    unsigned long n_blocks =
        lattice->mynode.volume() / (N_threads * iters_per_kernel) + 1;
    unsigned long n_sites = N_threads * n_blocks * iters_per_kernel;
    unsigned long myseed = seed + hila::myrank() * n_sites;
    cudaMalloc(&curandstate, n_sites * sizeof(curandState));
    check_cuda_error("seed_random malloc");
    seed_random_kernel<<<n_blocks, N_threads>>>(curandstate, myseed, iters_per_kernel,
                                                n_blocks * N_threads);
    check_cuda_error("seed_random kernel");

}

/* Generate random numbers on device or host */
__device__ __host__ double hila::random() {
#ifdef __CUDA_ARCH__
    int x = threadIdx.x + blockIdx.x * blockDim.x;
    return curand_uniform(&d_curandstate[x]);
#else
    return mersenne();
#endif
}

// Then, define global functions loop_lattice_size() and _volume()
__device__ __host__ int loop_lattice_size(Direction dir) {
#ifdef __CUDA_ARCH__
    return _d_size[dir];
#else
    return lattice->size(dir);
#endif
}
__device__ __host__ CoordinateVector loop_lattice_size(void) {
#ifdef __CUDA_ARCH__
    CoordinateVector v;
    foralldir(d) v[d] = _d_size[d];
    return v;
#else
    return lattice->size();
#endif
}
__device__ __host__ int64_t loop_lattice_volume(void) {
#ifdef __CUDA_ARCH__
    return _d_volume;
#else
    return lattice->volume();
#endif
}

void backend_lattice_struct::setup(lattice_struct *lattice) {
    CoordinateVector *tmp;

    /* Setup neighbour fields in all directions */
    for (int d = 0; d < NDIRS; d++) {
        // For normal boundaries
        cudaMalloc((void **)&(d_neighb[d]), lattice->mynode.volume() * sizeof(unsigned));
        check_cuda_error("cudaMalloc device neighbour array");
        cudaMemcpy(d_neighb[d], lattice->neighb[d],
                   lattice->mynode.volume() * sizeof(unsigned), cudaMemcpyHostToDevice);
        check_cuda_error("cudaMemcpy device neighbour array");

        // For special boundaries
        cudaMalloc((void **)&(d_neighb_special[d]),
                   lattice->mynode.volume() * sizeof(unsigned));
        check_cuda_error("cudaMalloc device neighbour array");
        const unsigned *special_neighb =
            lattice->get_neighbour_array((Direction)d, BoundaryCondition::ANTIPERIODIC);
        cudaMemcpy(d_neighb_special[d], special_neighb,
                   lattice->mynode.volume() * sizeof(unsigned), cudaMemcpyHostToDevice);
        check_cuda_error("cudaMemcpy device neighbour array");
    }

    /* Setup the location field */
    cudaMalloc((void **)&(d_coordinates),
               lattice->mynode.volume() * sizeof(CoordinateVector));
    check_cuda_error("cudaMalloc device coordinate array");
    tmp = (CoordinateVector *)malloc(lattice->mynode.volume() * sizeof(CoordinateVector));
    for (int i = 0; i < lattice->mynode.volume(); i++)
        tmp[i] = lattice->coordinates(i);
    cudaMemcpy(d_coordinates, tmp, lattice->mynode.volume() * sizeof(CoordinateVector),
               cudaMemcpyHostToDevice);
    check_cuda_error("cudaMemcpy device coordinate array");
    free(tmp);

    // Other backend_lattice parameters
    field_alloc_size = lattice->field_alloc_size();

    int s[NDIM];
    foralldir(d) s[d] = lattice->size(d);
    cudaMemcpyToSymbol(_d_size, s, sizeof(int) * NDIM, 0, cudaMemcpyHostToDevice);
    int64_t v = lattice->volume();
    cudaMemcpyToSymbol(_d_volume, &v, sizeof(int64_t), 0, cudaMemcpyHostToDevice);
    check_cuda_error("cudaMemcpy to size or volume");
}

void initialize_cuda(int rank) {
    int n_devices, my_device;

    cudaGetDeviceCount(&n_devices);
    check_cuda_error("Could not get device count");
    // This assumes that each node has the same number of mpi ranks and GPUs
    my_device = rank % n_devices;

    output0 << "CUDA Devices accessible from node 0: " << n_devices << '\n';
    // TODO: this only for node 0?
    if (rank < 6) {
        hila::output << "Cuda: rank " << rank << " choosing device " << my_device << '\n';
    }
    if (numnodes() > 6) {
        output0 << "  + " << numnodes() - 6 << " more nodes\n";
    }

    cudaSetDevice(my_device);
}

void cuda_device_info() {
    if (hila::myrank() == 0) {
        const int kb = 1024;
        const int mb = kb * kb;

        int driverVersion, rtVersion;
        cudaDriverGetVersion(&driverVersion);
        cudaRuntimeGetVersion(&rtVersion);
        hila::output << "CUDA driver version: " << driverVersion << ", runtime "
                     << rtVersion << '\n';

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
        hila::output << "  Constant memory: " << props.totalConstMem / kb << "kB" << '\n';
        hila::output << "  Block registers: " << props.regsPerBlock << '\n';

        hila::output << "  Warp size:         " << props.warpSize << '\n';
        hila::output << "  Threads per block: " << props.maxThreadsPerBlock << '\n';
        hila::output << "  Max block dimensions: [ " << props.maxThreadsDim[0] << ", "
                     << props.maxThreadsDim[1] << ", " << props.maxThreadsDim[2] << " ]"
                     << '\n';
        hila::output << "  Max grid dimensions:  [ " << props.maxGridSize[0] << ", "
                     << props.maxGridSize[1] << ", " << props.maxGridSize[2] << " ]"
                     << '\n';
    }
}

void cuda_exit_on_error(const char *msg, const char *file, int line) {
    cudaError code = cudaGetLastError();
    if (cudaSuccess != code) {
        hila::output << "CUDA error: " << msg << " in file " << file << " line " << line
                     << '\n';
        hila::output << "CUDA error string: " << cudaGetErrorString(code) << "\n";

        hila::terminate(0);
    }
}

void cuda_exit_on_error(cudaError code, const char *msg, const char *file, int line) {
    if (cudaSuccess != code) {
        hila::output << "CUDA error: " << msg << " in file " << file << " line " << line
                     << '\n';
        hila::output << "CUDA error string: " << cudaGetErrorString(code) << "\n";

        hila::terminate(0);
    }
}

#endif
