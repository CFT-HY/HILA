#include "../plumbing/defs.h"
#include "../plumbing/lattice.h"
#include "../plumbing/hila_cuda.h"

/* Random number generator */
curandState * curandstate;
__device__ curandState * d_curandstate;

__device__ unsigned * d_neighb[NDIRS];

/* Set seed on device */
__global__ void seed_random_kernel( curandState * state, unsigned long seed )
{
  int id = threadIdx.x;
  d_curandstate = state;
  curand_init ( seed, id, 0, &d_curandstate[id] );
}

/* Set seed on device and host */
void seed_random(unsigned long seed){
  cudaMalloc( &curandstate, N_threads*sizeof( curandState ) );
  check_cuda_error("seed_random malloc");
  seed_random_kernel<<< 1, N_threads >>>( curandstate, seed );
  check_cuda_error("seed_random kernel");
  seed_mersenne(seed+N_threads);
}

/* Generate random numbers on device or host */
#pragma transformer loop_function
double hila_random(){
  #ifdef __CUDA_ARCH__
  return curand_uniform( &d_curandstate[threadIdx.x] );
  #else
  return mersenne();
  #endif
}


void backend_lattice_struct::setup(lattice_struct lattice)
{
  location * tmp;

  /* Setup neighbour fields in all directions */
  for (int d=0; d<NDIRS; d++) {
    cudaMalloc( (void **)&(d_neighb[d]), lattice.local_volume() * sizeof(unsigned));
    check_cuda_error("cudaMalloc device neighbour array");

    cudaMemcpy( d_neighb[d], lattice.neighb[d], lattice.local_volume() * sizeof(unsigned), cudaMemcpyHostToDevice );
    check_cuda_error("cudaMemcpy device neighbour array");
  }

  /* Setup the location field */
  cudaMalloc( (void **)&(d_coordinates), lattice.local_volume() * sizeof(location));
  check_cuda_error("cudaMalloc device coordinate array");
  tmp = (location*) malloc( lattice.local_volume() * sizeof(location) );
  for(int i=0; i<lattice.local_volume(); i++) tmp[i] = lattice.coordinates(i);
  cudaMemcpy( d_coordinates, tmp, lattice.local_volume() * sizeof(location), cudaMemcpyHostToDevice );
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




