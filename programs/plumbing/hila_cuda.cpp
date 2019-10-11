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

__global__ void set_neighbour_pointers_kernel( unsigned * neighb, int d )
{
  d_neighb[d] = neighb;
}

void set_neighbour_pointers( unsigned * neighb, int d ){
  set_neighbour_pointers_kernel<<< 1, 1 >>>( neighb, d );
  check_cuda_error("set_neighbour_pointers kernel");
}

/* Generate on device or host */
loop_callable double hila_random(){
  #ifdef __CUDA_ARCH__
  return curand_uniform( &d_curandstate[threadIdx.x] );
  #else
  return mersenne();
  #endif
}



