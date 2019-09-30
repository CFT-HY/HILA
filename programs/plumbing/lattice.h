#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>
#include <fstream>

// TODO: assertion moved somewhere where basic params
#undef NDEBUG
#include <assert.h>
#include "../plumbing/defs.h"
#include "../plumbing/memory.h"

struct node_info {
  location min,size;
  unsigned evensites, oddsites;
};


class lattice_struct {
private:
  // expose these directly, by far the simplest interface - who cares about c++ practices
  // use also ints instead of unsigned, just to avoid surprises in arithmetics
  // I shall assume here that int is 32 bits, and long long 64 bits.  I guess these are
  // pretty much standard for now
  // Alternative: int_32t and int_64t (or int_fast_32t  and int_fast_64t, even more generally) 
  int l_size[NDIM];
  long long l_volume;

  // Information about the node stored on this process
  struct node_struct {
    unsigned index;
    unsigned sites, evensites, oddsites;
    unsigned field_alloc_size;          // how many sites/node in allocations 
    location min, size;                 // node local coordinate ranges
    unsigned nn[NDIRS];                 // nn-node of node down/up to dirs
    bool first_site_even;               // is location min even or odd?
    
    void setup(node_info & ni);
  } this_node;

  // information about all nodes
  struct allnodes {
    unsigned number;
    unsigned ndir[NDIM];  // number of node divisions to dir
    // lattice division: div[d] will have num_dir[d]+1 elements, last size
    // TODO: is this needed at all?
    std::vector<unsigned> divisors[NDIM];
    std::vector<node_info> nodelist;

    unsigned * map_array;                  // mapping (optional)
    unsigned * map_inverse;                // inv of it
    
    void create_remap();                   // create remap_node
    unsigned remap(unsigned i);            // use remap
    unsigned inverse_remap(unsigned i);    // inverse remap
    
  } nodes;

  struct comm_struct {
    
  };
  
  std::vector<comm_struct> commlist;
  
public:
  
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
  
  int size(direction d) { return l_size[d]; }
  long long volume() { return l_volume; }
  
  bool is_on_node(const location & c);
  unsigned node_number(const location & c);
  unsigned site_index(const location & c);
  unsigned site_index(const location & c, const unsigned node);
  
  const int field_alloc_size(){
    return this_node.field_alloc_size;
  }
  
  
};

extern lattice_struct * lattice;

#endif
