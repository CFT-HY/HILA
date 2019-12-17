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

struct vectorized_lattice_struct;


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

  std::vector<vectorized_lattice_struct*> vectorized_lattices;

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
  int local_size(int d) { return this_node.size[d]; }
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
  
  bool first_site_even() { return this_node.first_site_even; };

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

  vectorized_lattice_struct * get_vectorized_lattice(int vector_size);

};

/// global handle to lattice
extern lattice_struct * lattice;


// Keep track of defined lattices
extern std::vector<lattice_struct*> lattices;





/// Splits the local lattice into equal sections for vectorization
/// Contains a list of neighbour arrays for each possible vector size
struct vectorized_lattice_struct  {
  public:
    std::array<unsigned*,NDIRS> neighbours;
    std::array<int,NDIM> size;
    std::array<int,NDIM> split;
    int sites, evensites, oddsites, alloc_size;
    unsigned vector_size;
    std::vector<location> coordinates;
    lattice_struct * lattice;
    bool first_site_even;

    unsigned get_index(location l){
      int s = (int)first_site_even; // start at 0 for even first, 1 for odd first
      int l_index = l[NDIM-1];
      for (int d=NDIM-2; d>=0; d--){
        l_index = l_index*size[d] + (l[d] +size[d])%size[d];
        s +=  (l[d] +size[d])%size[d];
      }
      if( s%2==0 ){
        l_index /= 2;
      } else {
        l_index = evensites + l_index/2;
      }
      return l_index;
    }


    void setup(lattice_struct * _lattice, int _vector_size) {
      // Initialize
      lattice =  _lattice;
      vector_size = _vector_size;
      first_site_even = lattice->first_site_even();

      for(int d=0; d<NDIM; d++){
        size[d] = lattice->local_size(d);
        split[d] = 1;
      }

      int v = vector_size;
      while( v > 1 ){
        // find longest that will be even after split (divisible by 4)
        int msize=1, d=0;
        for( int i=0; i<NDIM; i++ ){
          if( size[i] > msize && size[i]%4==0 )
            d=i;
        }
        // split
        v /= 2; size[d] /= 2; split[d] *= 2;
      }
      assert(v == 1 && "cannot split local lattice to vector size");

      sites = 1;
      for( int i=0; i<NDIM; i++ ){
        sites *= size[i];
      }
      if( sites%2 == 0 ){
        evensites = oddsites = sites/2;
      } else {
        evensites = sites/2 + (int)first_site_even;
        oddsites  = sites/2 + (int)(!first_site_even);      
      }


      // Map index to local coordinates
      coordinates.resize(sites);
      for(unsigned i = 0; i<sites; i++){
        location l;
        unsigned l_index=i, s;
        foralldir(d){
          l[d] = l_index % size[d];
          l_index /= size[d];
        }
        coordinates[get_index(l)] = l;
      }


      // Setup neighbour array
      int halo_index = 0;
      std::vector<unsigned> orig_site(sites);
      for(int d=0; d<NDIM; d++){
        neighbours[d] = (unsigned *) malloc(sizeof(unsigned)*sites);
        for(unsigned i = 0; i<sites; i++){
          location l = coordinates[i];
          location nb = l;
          int k;
          if (is_up_dir(d)) {
            k = d;
            nb[d] = l[d] + 1;
          } else {
            k = opp_dir(d);
            nb[k] = l[k] - 1;
          }
          if( nb[k] > 0 && nb[k] < size[k] ) {
            neighbours[d][i] = get_index(nb);
          } else {
            // This is outside this split lattice (partly outside the node)
            // Fetching has already been set up
            neighbours[d][i] = sites + halo_index;
            halo_index++;
          }
        }
      }
      alloc_size = sites + halo_index;
    }


    // Return a list of neighbours for a lattice divided into a given vector size
    std::array<unsigned*,NDIRS> neighbour_list(){
      return neighbours;
    }

    int loop_begin( parity P) {
      if(P==ODD){
        return evensites;
      } else {
        return 0;
      }
    }

    int loop_end( parity P) {
      if(P==EVEN){
        return evensites;
      } else {
        return sites;
      }
    }

};




#endif
