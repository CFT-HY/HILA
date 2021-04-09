/// Memory allocator -- gives back aligned memory, if ALIGN defined

#include "plumbing/defs.h"

/// We'll have two prototypes, 2nd gives the file name and size for error messages
/// preprocessor substitutes memalloc() -calls with right parameters!

#if defined(VECTORIZED) || defined(AVX)
#define ALIGNED_MEMALLOC
#endif

void *memalloc(std::size_t size);
void *memalloc(std::size_t size, const char *filename, const unsigned line);

/// d_malloc allocates memory from "device": either from cpu or from gpu,
/// depending on the target.  Free with d_free()
void *d_malloc(std::size_t size);
void d_free(void * dptr);
