#ifndef _BACKEND_LATTICE_H_
#define _BACKEND_LATTICE_H_

#include "../lattice.h"
#include "vector_types.h"

#ifndef SUBNODE_LAYOUT
SUBNODE_LAYOUT needs to be defined to use this
#endif

/// This stores information needed for vectorized access to lattice fields.
/// - loop_begin, loop_end
/// - neighbour arrays
/// - twist permutations when needed
/// - allocation size, including new halo sites (replaces field_alloc_size)
/// these are used by fields which are laid out with vector structure


template<int vector_size>
struct vectorized_lattice_struct  {
  public:
    lattice_struct * lattice;

    unsigned v_sites;                                // vector sites on this node
    coordinate_vector subdivisions;                  // subnode divisions to different directions
    coordinate_vector subnode_origin, subnode_size;  // origin of the 1st subnode = origin of this_node
     
    // offsets to boundary halos
    bool is_boundary_permutation[NDIM];
    int boundary_permutation[NDIRS][vector_size];
    unsigned halo_offset[NDIRS], halo_offset_odd[NDIRS], n_halo_vectors[NDIRS];
    unsigned * RESTRICT halo_index[NDIRS];

    // move data from receive buffer -- sending is fine as it is
    // takes the role of nn_comms

    unsigned * recv_list[NDIRS];
    unsigned recv_list_size[NDIRS];
    
    // coordinate offsets to nodes
    typename vector_base_type<int,vector_size>::type coordinate_offset[NDIM];
    coordinate_vector * RESTRICT coordinate_base;

    unsigned * RESTRICT neighbours[NDIRS];
    unsigned alloc_size;


    bool is_on_first_subnode( const coordinate_vector & v ) {
      foralldir(d) {
        if ( v[d] < subnode_origin[d] || 
             v[d] >= subnode_origin[d]+subnode_size[d]) 
          return false;
      }
      return true;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /// Set up the vectorized lattice:
    /// * Split the lattice and record size and splits
    /// * Map indeces into local coordinate_list
    /// * Set up neighbour vector references
    /// mapping between the sites:  here i is the original lattice site index
    ///   vector_index = i / vector_size;
    ///   index_in_vector = i % vector_size;
    /// Note that this maps site indices of 64 bit elements (double)and 32 bit (float) so that
    ///   vector_index(32) = vector_index(64)/2 and
    ///   index_in_vector(32) = index_in_vector(64) + (vector_index(64)%2)*vector_size(64)
    /// This removes the last direction where number of subnodes is divisible
    //////////////////////////////////////////////////////////////////////////////////////////////

    vectorized_lattice_struct(lattice_struct * _lattice) {
      // Initialize
      lattice =  _lattice;

      /// vector
      v_sites = lattice->local_volume() / vector_size;
      subdivisions   = lattice->this_node.subnodes.divisions;
      subnode_size   = lattice->this_node.subnodes.size;
      subnode_origin = lattice->this_node.min;

      if ( vector_size == VECTOR_SIZE/sizeof(double) ) {
        for (int i=NDIM-1; i>=0; i--) {
          if (lattice->this_node.subnodes.divisions[i] > 1) {
            subdivisions[i] /= 2;
            subnode_size[i] *= 2;
            break;
          }
        }
      }

      hila::output << "Setting up lattice struct with vector of size " << vector_size << '\n';


      // boundary permutation maps subnodes to vectors
      // fastest moving division to smallest dimension
      // a 3-dim. 4x4 node is divided into 8 subnodes as follows:  
      // here 0-3 is the index within subnode, and a-d the subnode label.
      //    0a 1a | 0b 1b    4a 5a | 4b 5b   |   0e 1e | 0f 1f    4e 5e | 4f 5f     
      //    2a 3a | 2b 3b    6a 7a | 6b 7b   |   2e 3e | 2f 3f    6e 7e | 6f 7f     
      //    -------------    -------------   |   -------------    -------------     
      //    0c 1c | 0d 1d    4c 5c | 4d 5d   |   0g 1g | 0h 1h    4g 5g | 4h 5h     
      //    2c 3c | 2d 3d    6c 7c | 6d 7d   |   2g 3g | 2h 3h    6g 7g | 6h 7h 
      //      
      //  boundary_permutation is ab ba cd dc ef fe gh hg  to 0-dir (horizontal)
      //                      and ac bd ca db eg fh ge hf  to 1-dir (vertical)
      //                      and ae bf cg dh ea fb gc hd  to 2-dir
      //
      int step = 1, sdmul = 1;
      foralldir(d) {
        is_boundary_permutation[d] = (subdivisions[d] > 1);

        sdmul *= subdivisions[d];
        for (int i=0; i<vector_size; i++) 
          boundary_permutation[d][i] = (i+step) % sdmul + (i/sdmul) * sdmul;
        step *= subdivisions[d];

        // permutation to oppsitie direction is the inverse, thus:
        for (int i=0; i<vector_size; i++) 
          boundary_permutation[-d][ boundary_permutation[d][i] ] = i;
      }



      // reserve extra storage for permutation halo sites
      // there are 2*(area) v-sites to each direction with permutation,
      //  + mpi buffers too if needed. (separately allocated)
      // to directions without permutation there is
      //  - no halo site, if no node communication
      //  - 2*(area) v-sites, filled in directly by MPI or other comm.
      // These come automatically when we tally up the neigbours below

      for (direction d=(direction)0; d<NDIRS; d++) {
        neighbours[d] = (unsigned *)memalloc(v_sites * sizeof(unsigned));
      }

      // accumulate here points off-subnode (to halo)
      int c_offset = v_sites;  
      for(direction d=(direction)0; d<NDIRS; ++d) {
        halo_offset[d] = c_offset;
        for (int i=0; i<v_sites; i++) {
          int j = vector_size*i;   // the "original lattice" index for the 1st site of vector
          coordinate_vector here = lattice->coordinates(j);
          // std::cout << here << '\n';

          if (is_on_first_subnode(here+d)) {
            assert(lattice->neighb[d][j] % vector_size == 0);   // REMOVE THIS LATER, consistency check
            neighbours[d][i] = lattice->neighb[d][j]/vector_size;
          } else {
            neighbours[d][i] = c_offset++;  // now points beyond the lattice
          }
        }
        n_halo_vectors[d] = c_offset - halo_offset[d];
        halo_offset_odd[d] = halo_offset[d] + n_halo_vectors[d]/2;
        assert(n_halo_vectors[d] % 2 == 0);

        // set also the index array, if needed
        // halo_index[d] points to the subnode modded neighbour site to dir d, if 
        // there is boundary twist (meaning there are on-node subnodes to this dir)
        // we'll use the standard neighb array to do this.
        if (n_halo_vectors[d] > 0 && is_boundary_permutation[abs(d)]) {
          halo_index[d] = (unsigned *)memalloc(n_halo_vectors[d] * sizeof(unsigned));
          int j=0;
          for (int i=0; i<v_sites; i++) {
            if (neighbours[d][i] >= v_sites) {
              halo_index[d][j++] = lattice->neighb[d][i*vector_size]/vector_size;
            }
          }
          assert( j == n_halo_vectors[d] );
        }
      }

      /// and then possible neighbour receive indices
      for (direction d=(direction)0; d<NDIRS; d++) {
        if (is_boundary_permutation[abs(d)] && lattice->nodes.n_divisions[abs(d)] > 1) {

          // now need to receive and copy - note: now this is in terms of 
          // non-vector sites.   Set the recv_list to point to where to move the stuff
          // Note: now the stuff has to be moved to boundary_halo, not to lattice n!
    
          recv_list_size[d] = lattice->this_node.sites / lattice->this_node.size[abs(d)];
          recv_list[d] = (unsigned *)memalloc( recv_list_size[d] * sizeof(unsigned) );

          int j=0;
          for (int i=0; i<lattice->this_node.sites; i++) {
            if (lattice->neighb[d][i] >= lattice->this_node.sites) {
          
              // i/vector_size is the "vector index" of site, and 
              // i % vector_size the index within the vector.  
              // vector neighbour is neighbours[d][i/vector_size]

              recv_list[d][j++] = neighbours[d][i/vector_size]*vector_size + i % vector_size;
            }
          }
          assert(j == recv_list_size[d]);

        } else {
          // now use halo_offset directly for buffer 
          recv_list[d] = nullptr;
          recv_list_size[d] = 0;
        }

      }

      /// how big the field allocation should be - IN SITES, not vectors
      alloc_size = c_offset * vector_size;

      /// and set the coordinates
      /// first vector_size elements should give the coordinates of vector offsets
      coordinate_vector base = lattice->coordinates(0);
      for (int i=0; i<vector_size; i++) {
        coordinate_vector diff = lattice->coordinates(i) - base;
        foralldir(d) 
          coordinate_offset[d].insert(i, diff[d]);
      }

      // and then set the coordinate_base with the original coords
      coordinate_base = (coordinate_vector *)memalloc(v_sites * sizeof(coordinate_vector));
      for (int i=0; i<v_sites; i++) {
        coordinate_base[i] = lattice->coordinates(vector_size*i);
      }


    }  // end of initialization

    /// Return the communication info
    lattice_struct::nn_comminfo_struct get_comminfo(int d){
      return lattice->get_comminfo(d);
    }

    /// get neighbours for this, with 2 different methods:
    unsigned vector_neighbour(direction d, int idx) {
      return neighbours[d][idx];
    }

    /// this gives the neighbour when the lattice is traversed
    /// site-by-site.  Now idx is the "true" site index, not vector index
    unsigned site_neighbour(direction d, int idx) {
      return vector_size * neighbours[d][idx/vector_size] + idx % vector_size;
    }

    /// Return the number of sites that need to be allocated 
    /// returns sites, not vectors!
    unsigned field_alloc_size() const {
      return alloc_size;
    }


    /// Return the coordinates of each vector nested as
    /// coordinate[direction][vector_index]
    auto coordinates(int idx) {
      std::array<typename vector_base_type<int,vector_size>::type ,NDIM> r;
      foralldir(d) r[d] = coordinate_offset[d] + coordinate_base[idx][d];
      return r;
    }

    auto coordinate(direction d, int idx) {
      return coordinate_offset[d] + coordinate_base[idx][d];
    }

    // parity is the same for all elements in vector, return scalar
    ::parity parity(int idx) {
      return coordinate_base[idx].parity();
    }

    /// First index in a lattice loop
    int loop_begin( ::parity P) {
      if(P==ODD){
        return v_sites/2;
      } else {
        return 0;
      }
    }

    // Last index in a lattice loop
    int loop_end( ::parity P) {
      if(P==EVEN){
        return v_sites/2;
      } else {
        return v_sites;
      }
    }

    template <typename T>
    void reduce_node_sum(T * value, int N, bool distribute){
      lattice->reduce_node_sum(value, N, distribute);
    };

    template <typename T>
    void reduce_node_product(T * value, int N, bool distribute){
      lattice->reduce_node_product(value, N, distribute);
    };

};




struct backend_lattice_struct {
  lattice_struct * latticep;

  void setup(lattice_struct & _lattice){
    latticep = &_lattice;
  }

  template< int vector_size >
  vectorized_lattice_struct<vector_size> * get_vectorized_lattice() {
    // Create one if not already created
    static bool init = true;
    static vectorized_lattice_struct<vector_size> * vlat; 
    if(init){
      vlat = new vectorized_lattice_struct<vector_size>(latticep);
      init = false;
    }

    return vlat;
  }
};




#endif