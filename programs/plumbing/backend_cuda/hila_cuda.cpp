#include "../../plumbing/defs.h"
#include "../../plumbing/lattice.h"
#include "../../plumbing/backend_cuda/defs.h"

/* Random number generator */
curandState * curandstate;
__device__ curandState * d_curandstate;

__device__ unsigned * d_neighb[NDIRS];

/* Set seed on device */
__global__ void seed_random_kernel( curandState * state, unsigned long seed )
{
  int x = threadIdx.x + blockIdx.x * blockDim.x;
  d_curandstate = state;
  curand_init( seed, x, 0, &d_curandstate[x] );
}

/* Set seed on device and host */
void seed_random(unsigned long seed){
  unsigned long n_blocks = lattice->local_volume() / N_threads + 1;
  unsigned long n_sites = N_threads*n_blocks;
  unsigned long myseed = seed + mynode()*lattice->local_volume();
  cudaMalloc( &curandstate, n_sites*sizeof( curandState ) );
  check_cuda_error("seed_random malloc");
  seed_random_kernel<<< n_blocks, N_threads >>>( curandstate, myseed );
  check_cuda_error("seed_random kernel");
  seed_mersenne(seed+n_sites);
}

/* Generate random numbers on device or host */
#pragma transformer loop_function
double hila_random(){
  #ifdef __CUDA_ARCH__
  int x = threadIdx.x + blockIdx.x * blockDim.x;
  return curand_uniform( &d_curandstate[x] );
  #else
  return mersenne();
  #endif
}


void backend_lattice_struct::setup(lattice_struct lattice)
{
  coordinate_vector * tmp;

  /* Setup neighbour fields in all directions */
  for (int d=0; d<NDIRS; d++) {
    cudaMalloc( (void **)&(d_neighb[d]), lattice.local_volume() * sizeof(unsigned));
    check_cuda_error("cudaMalloc device neighbour array");

    cudaMemcpy( d_neighb[d], lattice.neighb[d], lattice.local_volume() * sizeof(unsigned), cudaMemcpyHostToDevice );
    check_cuda_error("cudaMemcpy device neighbour array");
  }

  /* Setup the location field */
  cudaMalloc( (void **)&(d_coordinates), lattice.local_volume() * sizeof(coordinate_vector));
  check_cuda_error("cudaMalloc device coordinate array");
  tmp = (coordinate_vector*) malloc( lattice.local_volume() * sizeof(coordinate_vector) );
  for(int i=0; i<lattice.local_volume(); i++) tmp[i] = lattice.coordinates(i);
  cudaMemcpy( d_coordinates, tmp, lattice.local_volume() * sizeof(coordinate_vector), cudaMemcpyHostToDevice );
  check_cuda_error("cudaMemcpy device coordinate array");
  free(tmp);

  // Other backend_lattice parameters
  field_alloc_size = lattice.field_alloc_size();
}



void initialize_cuda( int rank ){
  int n_devices, my_device;
  cudaGetDeviceCount(&n_devices);
  /* This assumes that each node has the same number of mpi ranks and GPUs */
  my_device = rank%n_devices;
  printf("Rank %d choosing device %d out of %d\n", rank, my_device, n_devices);
  cudaSetDevice(my_device);
}




