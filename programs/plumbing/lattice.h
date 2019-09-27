#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>
#include <fstream>
#include <array>
#include <vector>

// TODO: assertion moved somewhere where basic params
#undef NDEBUG
#include <assert.h>
#include "../plumbing/defs.h"
#include "../plumbing/memory.h"

using location = std::array<int,NDIM>;

#ifndef USE_MPI
static inline int mynode(){
  return 0;
}
#endif


class lattice_struct {
public:
  // expose these directly, by far the simplest interface - who cares about c++ practices
  // use also ints instead of unsigned, just to avoid surprises in arithmetics
  // I shall assume here that int is 32 bits, and long long 64 bits.  I guess these are
  // pretty much standard for now
  // Alternative: int_32t and int_64t (or int_fast_32t  and int_fast_64t, even more generally) 
  int size[NDIM];
  long long volume;

  // Information about the node stored on this process
  struct node_struct {
    unsigned index;
    unsigned sites, evensites, oddsites;
    unsigned field_alloc_size;          // how many sites/node in allocations 
    unsigned min[NDIM], size[NDIM];          // node local coordinate ranges
    unsigned nn[NDIRS];                      // nn-node of node down/up to dirs
  } mynode;

  // information about all nodes
  struct allnodes {
    unsigned number;
    unsigned num_dir[NDIM];
    // lattice division to nodes: div[d] will have num_dir[d]+1 elements, last size
    std::vector<unsigned> div[NDIM];
    unsigned * remap;                    // mapping (optional)
  } nodes;

  #ifdef USE_MPI
  int * map_node_list;
  #endif  
  
  
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
  
  bool is_on_node(const location & c);
  unsigned node_number(const location & c);
  unsigned site_index(const location & c);
  
  unsigned remap_node(const unsigned i);

  unsigned loop_begin(const parity p);
  unsigned loop_end(const parity p);

  
};

extern lattice_struct * lattice;

#endif
