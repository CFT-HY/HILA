#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>
#include <fstream>

// TODO: assertion moved somewhere where basic params
#undef NDEBUG
#include <assert.h>

#include "memory.h"


// move these somewhere - use consts?
#define NDIM 4
#define NDIRS (2*NDIM)


class lattice_struct {
public:
  // expose these directly, by far the simplest interface - who cares about c++ practices
  // use also ints instead of unsigned, just to avoid surprises in arithmetics
  // I shall assume here that int is 32 bits, and long long 64 bits.  I guess these are
  // standard for now
  // Alternative: int_32t and int_64t (or int_fast_32t  and int_fast_64t, even more generally) 
  int size[NDIM];
  long long volume;

  struct node_struct {
    int sites, evensites, oddsites;
    int field_alloc_size;          // how many sites/node in allocations 
    
    int xmin[NDIM], size[NDIM];    // node local coordinate ranges
    int nn[NDIRS];                 // nn-node of node down/up to dirs
    
  } node;

  
  void setup(int siz[NDIM]);
  
#if NDIM == 4
  void setup(int nx, int ny, int nz, int nt);
#elif NDIM == 3  
  void setup(int nx, int ny, int nz);
#elif NDIM == 2
  void setup(int nx, int ny);
#elif NDIM == 1
  void setup(int nx); 
#endif
  
};

extern lattice_struct * lattice;

#endif
