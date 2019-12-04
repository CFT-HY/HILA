
#include <stdlib.h>
#include "defs.h"
#include "memory.h"
#include <iostream>

// Memory alignment for "vanilla"

#ifdef VANILLA

void * allocate_field_mem(size_t size) {
  // returns nullptr if there are problems
  void * p = malloc(size);
  if (p == nullptr) {
    std::cout << "Failure in field memory allocation\n";
    exit(1);
  }
  return p;
}

void free_field_mem(void * p) {
  if (p != nullptr)
    free(p);
}


#elif AVX512

// AVX512 sweet spot?
#define FIELD_ALIGNMENT 512

void * allocate_field_mem(size_t size) {
  // guarantee size is a multiple of alignment
  if (size % FIELD_ALIGNMENT) 
    size = size - (size % FIELD_ALIGNMENT) + FIELD_ALIGNMENT;

  // returns nullptr if there are problems
  void * p = aligned_alloc( FIELD_ALIGNMENT, size);
  return p;
}

/// 
void free_field_mem(void * p) {
  if (p != nullptr)
    free(p);
}


#elif defined(CUDA) && !defined(TRANSFORMER)
#include "../plumbing/hila_cuda.h"


void * allocate_field_mem(size_t size) {
  void * p;
  cudaMalloc( &p, size );
  check_cuda_error("Allocate field memory");
  if (p == nullptr) {
    std::cout << "Failure in field memory allocation\n";
    exit(1);
  }
  return p;
}

void free_field_mem(void * p) {
  if (p != nullptr){
    cudaFree(p);
    check_cuda_error("Free field memory");
  }
}


#endif
