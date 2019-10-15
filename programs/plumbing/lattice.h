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
struct node_info {
  location min,size;
  unsigned evensites, oddsites;
};


#ifdef CUDA
struct device_lattice_info {
  unsigned * d_neighb[NDIRS];
  unsigned field_alloc_size;
  int loop_begin, loop_end;
};
#endif


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
    std::vector<location> site_index_list;

    void setup(node_info & ni, lattice_struct & lattice);
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

  struct comm_node_struct {
    unsigned index;
    unsigned sites, evensites, oddsites;
    unsigned buffer;
    std::vector<unsigned>  sitelist;
  };

  struct comminfo_struct {
    int label;    
    unsigned * index;
    std::vector<comm_node_struct> from_node;
    std::vector<comm_node_struct> to_node;
  };
  
  std::vector<comminfo_struct> comminfo;

public:

  unsigned * neighb[NDIRS];

  #ifdef CUDA
  device_lattice_info device_info;
  #endif

  void setup(int siz[NDIM]);
  void setup_layout();
  void setup_nodes();
  
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
  int size(int d) { return l_size[d]; }
  long long volume() { return l_volume; }
  int node_number() { return this_node.index; }
  
  bool is_on_node(const location & c);
  unsigned node_number(const location & c);
  unsigned site_index(const location & c);
  unsigned site_index(const location & c, const unsigned node);
  location site_location(unsigned index);
  unsigned field_alloc_size() {return this_node.field_alloc_size; }
  void create_std_gathers();

  unsigned remap_node(const unsigned i);
  
  #ifdef EVENFIRST
  const int loop_begin( parity P){
    if(P==ODD){
      return this_node.evensites;
    } else {
      return 0;
    }
  }
  const int loop_end( parity P){
    if(P==EVEN){
      return this_node.evensites;
    } else {
      return this_node.sites;
    }
  }
  #else
  const int loop_begin( parity P){
    if(P==EVEN){
      return this_node.evensites;
    } else {
      return 0;
    }
  }
  const int loop_end( parity P){
    if(P==ODD){
      return this_node.evensites;
    } else {
      return this_node.sites;
    }
  }
  #endif


  /* Communication routines. Define here in lattice? */
  template <typename T>
  void reduce_node_sum(T value, bool distribute);

  template <typename T>
  void reduce_node_product(T value, bool distribute);

  
};

/// global handle to lattice
extern lattice_struct * lattice;

#ifdef USE_MPI
#include "comm_mpi.h"
#else
#include "comm_vanilla.h"
#endif


#endif
