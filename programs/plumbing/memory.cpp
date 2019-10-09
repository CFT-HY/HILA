
#include <stdlib.h>
#include "defs.h"
#include "memory.h"


// Memory alignment for "vanilla"

#ifdef VANILLA

void * allocate_field_mem(size_t size) {
  // returns nullptr if there are problems
  return malloc(size);  
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
  return aligned_alloc( FIELD_ALIGNMENT, size);
  
}

/// 
void free_field_mem(void * p) {
  if (p != nullptr)
    free(p);
}

#endif
