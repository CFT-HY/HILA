#ifndef MEMORY_H
#define MEMORY_H

// Definitions associated with field allocation

void allocate_field_mem(void ** p, size_t size);
void free_field_mem(void * p);

#endif
