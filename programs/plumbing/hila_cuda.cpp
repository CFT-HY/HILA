#include "defs.h"
#include "hila_cuda.h"

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
loop_callable double hila_random(){
  #ifdef __CUDA_ARCH__
  return curand_uniform( &d_curandstate[threadIdx.x] );
  #else
  return mersenne();
  #endif
}



void lattice_struct::setup_lattice_device_info(){

  /* Setup neighbour fields in all directions */
  for (int d=0; d<NDIRS; d++) {
    cudaMalloc( (void **)&(device_info.d_neighb[d]), this_node.sites * sizeof(unsigned));
    check_cuda_error("cudaMalloc device neighbour array");

    cudaMemcpy( device_info.d_neighb[d], neighb[d], this_node.sites * sizeof(unsigned), cudaMemcpyHostToDevice );
    check_cuda_error("cudaMemcpy device neighbour array");
  }

  /* Setup the location field */
  cudaMalloc( (void **)&(device_info.d_coordinates), this_node.sites * sizeof(location));
  check_cuda_error("cudaMalloc device coordinate array");
  cudaMemcpy( device_info.d_coordinates, this_node.site_index_list, this_node.sites * sizeof(location), cudaMemcpyHostToDevice );
  check_cuda_error("cudaMemcpy device coordinate array");

  // Other device_info parameters
  device_info.field_alloc_size = field_alloc_size();
}






