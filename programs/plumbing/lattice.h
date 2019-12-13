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
#include "../plumbing/inputs.h"

struct node_info {
  location min,size;
  unsigned evensites, oddsites;
};


#ifdef CUDA
struct device_lattice_info {
  unsigned * d_neighb[NDIRS];
  unsigned field_alloc_size;
  int loop_begin, loop_end;
  location * d_coordinates;

  /* Actual implementation of the coordinates function.
   * Since this is added by the transformer, the ASTVisitor
   * will not recognize it as a loop function. Since it's
   * purely for CUDA, we can just add the __device__ keyword. */
  #pragma transformer loop_function
  location coordinates( unsigned idx ){
    return d_coordinates[idx];
  }
};
#endif


/// Splits the local lattice into equal sections for vectorization
/// Contains a list of neighbour arrays for each possible vector size
struct vectorized_lattice_struct {
  std::vector<std::array<unsigned*,NDIRS>> neighbours;
  std::vector<unsigned> vector_sizes;
  std::vector<std::array<int,NDIM>> local_size;
  std::vector<std::array<int,NDIM>> splits;
  int l_size[NDIM];

  void set_size(int l[NDIM]) {
    for(int d=0; d<NDIM; d++)
      l_size[d] = l[d];
  }

  std::array<unsigned*,NDIRS> init_neighbours(int vector_size){
    std::array<int,NDIM> _local_size;
    std::array<int,NDIM> _splits;

    vector_sizes.push_back(vector_size);

    // Initialize
    for(int d=0; d<NDIM; d++){
      _local_size[d] = l_size[d];
      _splits[d] = 1;
    }
    while( vector_size > 1 ){
      // find longest direction
      int msize=1, d=0;
      for( int i=0; i<NDIM; i++ ){
        if( _local_size[i] > msize )
          d=i;
      }
      // split
      vector_size /= 2; _local_size[d] /= 2; _splits[d] *= 2;
    }
    assert(vector_size == 1 && "cannot handle vector size not divisible by 2");
    
    std::array<unsigned*,NDIRS> new_neighbours;
    for(int d=0; d<NDIM; d++){
      new_neighbours[d] = (unsigned *) malloc(sizeof(unsigned));
    }
    neighbours.push_back(new_neighbours);

    return new_neighbours;
  }

  /// Return a list of neighbours for a lattice divided into a given vector size
  std::array<unsigned*,NDIRS> neighbour_list(int vector_size){
    for( int i=0; i<neighbours.size(); i++ ){
      printf("vector size %d \n",vector_sizes[i]);
      if( vector_size == vector_sizes[i] ){
        return neighbours[i];
      }
    }
    std::array<unsigned*,NDIRS> new_neighbours = init_neighbours(vector_size);
    return(new_neighbours);
  }

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
    std::vector<location> coordinates;

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

public:

  struct comm_node_struct {
    unsigned index;
    unsigned sites, evensites, oddsites;
    unsigned buffer;
    std::vector<unsigned>  sitelist;

    // The number of sites that need to be communicated
    unsigned n_sites(parity par){
      if(par == ALL){
        return sites;
      } else if(par == EVEN){
        return evensites;
      } else {
        return oddsites;
      }
    }

    // The local index of a site that is sent to neighbour
    unsigned site_index(int site, parity par){
      if(par == ODD){
        return sitelist[evensites+site];
      } else {
        return sitelist[site];
      }
    }

    // The offset of the halo from the start of the field array
    unsigned offset(parity par){
      if(par == ODD){
        return buffer + evensites;
      } else {
        return buffer;
      }
    }
  };

  struct comminfo_struct {
    int label;    
    unsigned * index;
    std::vector<comm_node_struct> from_node;
    std::vector<comm_node_struct> to_node;
  };

  std::vector<comminfo_struct> comminfo;

  vectorized_lattice_struct vectorized_lattice;

  unsigned * neighb[NDIRS];
  unsigned char *wait_arr_;

  #ifdef CUDA
  device_lattice_info device_info;
  void setup_lattice_device_info();
  #endif

  void setup(int siz[NDIM], int &argc, char **argv);
  void setup(input & inputs);
  void setup_layout();
  void setup_nodes();
  
  #if NDIM == 4
  void setup(int nx, int ny, int nz, int nt, int &argc, char **argv);
  #elif NDIM == 3  
  void setup(int nx, int ny, int nz, int &argc, char **argv);
  #elif NDIM == 2
  void setup(int nx, int ny, int &argc, char **argv);
  #elif NDIM == 1
  void setup(int nx, int &argc, char **argv); 
  #endif


  void teardown();

  int size(direction d) { return l_size[d]; }
  int size(int d) { return l_size[d]; }
  long long volume() { return l_volume; }
  int node_number() { return this_node.index; }
  int n_nodes() { return nodes.number; }
  
  bool is_on_node(const location & c);
  unsigned node_number(const location & c);
  unsigned site_index(const location & c);
  unsigned site_index(const location & c, const unsigned node);
  location site_location(unsigned index);
  unsigned field_alloc_size() {return this_node.field_alloc_size; }
  void create_std_gathers();

  unsigned remap_node(const unsigned i);
  
  #ifdef EVENFIRST
  int loop_begin( parity P) const {
    if(P==ODD){
      return this_node.evensites;
    } else {
      return 0;
    }
  }
  int loop_end( parity P) const {
    if(P==EVEN){
      return this_node.evensites;
    } else {
      return this_node.sites;
    }
  }
  #else
  int loop_begin( parity P) const {
    if(P==EVEN){
      return this_node.evensites;
    } else {
      return 0;
    }
  }
  int loop_end( parity P) const {
    if(P==ODD){
      return this_node.evensites;
    } else {
      return this_node.sites;
    }
  }
  #endif

  location coordinates( unsigned idx ){
    return site_location(idx);
  }

  /* MPI functions and variables. Define here in lattice? */
  void initialize_wait_arrays();
  #ifdef USE_MPI
  MPI_Comm mpi_comm_lat;
  #endif

  template <typename T>
  void reduce_node_sum(T & value, bool distribute);

  template <typename T>
  void reduce_node_product(T & value, bool distribute);

  // Guarantee 64 bits for these - 32 can overflow!
  unsigned long long n_gather_done = 0, n_gather_avoided = 0;

};

/// global handle to lattice
extern lattice_struct * lattice;


// Keep track of defined lattices
extern std::vector<lattice_struct*> lattices;


#endif
