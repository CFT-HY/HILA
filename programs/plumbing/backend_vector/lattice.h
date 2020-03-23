#ifndef _BACKEND_LATTICE_H_
#define _BACKEND_LATTICE_H_


/// Splits the local lattice into equal sections for vectorization
template<int vector_size>
struct vectorized_lattice_struct  {
  public:
    std::array<int,NDIM> size;
    std::array<int,NDIM> split;
    int sites, evensites, oddsites, alloc_size;
    std::vector<coordinate_vector> coordinate_list;
    coordinate_vector min; // Coordinate of first site
    std::vector<int**> coordinate;
    lattice_struct * lattice;
    bool first_site_even;

    // Map full lattice index into a vectorized lattice index
    // and a vector index
    std::vector<int> lattice_index;
    std::vector<int> vector_index;

    std::array<unsigned*,NDIRS> neighbours;
    std::array<int*, NDIRS> boundary_permutation;

    struct halo_sites_struct {
      unsigned first_index;
      std::vector<unsigned> nb_index;
    };
    halo_sites_struct halo_sites[2][NDIRS];


    /// Set up the vectorized lattice:
    /// * Split the lattice and record size and splits
    /// * Map indeces into local coordinate_list
    /// * Set up neighbour vector references
    vectorized_lattice_struct(lattice_struct * _lattice) {
      // Initialize
      lattice =  _lattice;
      first_site_even = lattice->first_site_even();
      min = lattice->coordinates(0);

      for(int d=0; d<NDIM; d++){
        size[d] = lattice->local_size(d);
        split[d] = 1;
      }

      int v = vector_size;
      while( v > 1 ){
        // find longest that will be even after split (divisible by 4)
        int msize=1, d=-1;
        for( int i=0; i<NDIM; i++ ){
          if( size[i] > msize && size[i]%4==0 ) {
            d=i;
            msize = size[i];
          }
        }
        assert(d>=0 && "cannot split local lattice to vector size");
        // split
        v /= 2; size[d] /= 2; split[d] *= 2;
      }

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

      // Map index to local coordinate_list
      coordinate_list.resize(sites);
      coordinate.resize(sites);
      for(unsigned i = 0; i<sites; i++){
        coordinate_vector l;
        unsigned l_index=i;
        for(int d=0; d<NDIM; d++){
          l[d] = l_index % size[d];
          l_index /= size[d];
        }
        int index = get_index(l);
        coordinate_list[index] = l;
      }

      // Setup neighbour array and halo copy
      int halo_index = 0;
      for(int d=0; d<NDIRS; d++) {
        neighbours[d] = (unsigned *) malloc(sizeof(unsigned)*sites);
        int updir;
        if (is_up_dir(d)) {
          updir = d;
        } else {
          updir = opp_dir(d);
        }

        // If the lattice is split by MPI or by vectorization we need
        // to create a halo copy of neighbours that cross the boundary
        // Check here if that is the case for this direction
        bool need_halo = false;
        if( lattice->comminfo[d].from_node.size() > 0 
            || split[updir] > 1 ){
          need_halo = true;
        }

        // Loop over parities
        for( int par_int = 0; par_int < 2; par_int++ ){ 
          int first_index = halo_index;
          halo_sites[par_int][d].first_index = first_index;
          std::vector<unsigned> halo_temp;
          std::vector<unsigned> & nb_index = halo_sites[par_int][d].nb_index;
          
          // Check each site in the parity (for each direction separately)
          for(unsigned i = par_int*evensites; i<evensites+par_int*oddsites; i++){
            coordinate_vector l = coordinate_list[i];
            coordinate_vector nb = l;
            if (is_up_dir(d)) {
              nb[d] = l[d] + 1;
            } else {
              nb[updir] = l[updir] - 1;
            }
            if( nb[updir] >= 0 && nb[updir] < size[updir] ) {
              neighbours[d][i] = get_index(nb);
            } else {
              // This is outside this split lattice (maybe partly outside the node)              
              if( need_halo ) {
                // Add to halo site list
                nb_index.push_back(get_index(nb));
                halo_temp.push_back(i);
                halo_index++;
              } else {
                // It's accross the boundary but we can treat it like
                // a normal neighbour
                neighbours[d][i] = get_index(nb);
              }
            }
          }

          // Bubble sort the halo array according to the
          // neighbor index
          if(nb_index.size() > 0)
          for (int i=0; i < nb_index.size(); i++) {
            for (int k=0; k < nb_index.size()-i-1; k++) {
              if( nb_index[k] > nb_index[k+1] ){
                int t;
                t=nb_index[k];
                nb_index[k] = nb_index[k+1];
                nb_index[k+1] = t;
                t=halo_temp[k];
                halo_temp[k] = halo_temp[k+1];
                halo_temp[k+1] = t;
              }
            }
          }

          for( int i=0; i<halo_temp.size(); i++ ) {
            neighbours[d][halo_temp[i]] = sites + first_index + i;
          }
          
        }
      }
      alloc_size = sites + halo_index;


      // Setup permutation vectors for setting the halo
      for(int d=0; d<NDIM; d++){
        boundary_permutation[d] = (int *) malloc(sizeof(int)*vector_size);
        boundary_permutation[opp_dir(d)] = (int *) malloc(sizeof(int)*vector_size);
      }
      // This needs to be done only once, so use the most straightforward method
      // Check each vector index in each direction separately. Find the coordinate_vector
      // of the vectorized sublattice and figure out each neighbour from there.
      for(int v=0; v<vector_size; v++){
        coordinate_vector vl, nb;
        int v_ind = v;
        // Find coordinate_vector of this vector index
        for(int d=0; d<NDIM; d++){
          vl[d] = v_ind % split[d];
          v_ind /= split[d];
        }
        for(int d=0; d<NDIM; d++){
          // Find neighbour in positive direction
          for(int d2=0; d2<NDIM; d2++)
            nb[d2] = vl[d2];
          nb[d] = (vl[d] + 1) % split[d];
          int nb_index = nb[0];
          int step = split[0];
          for(int d2=1; d2<NDIM; d2++){
            nb_index += step*nb[d2];
            step *= split[d2];
          }
          boundary_permutation[d][v] = nb_index;

          // And negative direction
          for(int d2=0; d2<NDIM; d2++)
            nb[d2] = vl[d2];
          nb[d] = (vl[d] + split[d] - 1) % split[d];
          nb_index = nb[0];
          step = split[0];
          for(int d2=1; d2<NDIM; d2++){
            nb_index += step*nb[d2];
            step *= split[d2];
          }
          boundary_permutation[opp_dir(d)][v] = nb_index;
        }
      }
      

      // Map full lattice index to local index
      lattice_index.resize(lattice->field_alloc_size());
      vector_index.resize(lattice->field_alloc_size());
      for(unsigned i = lattice->loop_begin(ALL); i<lattice->loop_end(ALL); i++){
        coordinate_vector fl = lattice->local_coordinates(i);
        coordinate_vector vl;
        int step=1, v_index=0;
        for( int d=0; d<NDIM; d++ ){
          vl[d] = fl[d] % size[d];
          v_index += step * (fl[d] / size[d]);
          step *= split[d];
        }
        lattice_index[i] = get_index(vl);
        vector_index[i] = v_index;

        // Check for neighbours in the halo
        for(int d=0; d<NDIRS; d++){
          int fl_index = lattice->neighb[d][i];
          if( fl_index >= lattice->local_volume() ){
            // There is an off-node neighbour. Map this to a halo-index
            int l_index = neighbours[d][get_index(vl)];
            lattice_index[fl_index] = l_index;
            vector_index[fl_index] = v_index; 
          }
        }
      }
    }


    /// Return the communication info
    lattice_struct::comminfo_struct get_comminfo(int d){
      return lattice->get_comminfo(d);
    }


    /// Return the number of sites that need to be allocated
    /// (1 vector for each site)
    unsigned field_alloc_size() {
      return alloc_size;
    }


    /// Return the coordinates of each vector nested as
    /// coordinate[direction][vector_index]
    auto coordinates(int idx);


    /// Translate a local coordinate_vector into an index 
    unsigned get_index(coordinate_vector l){
      int s = 1-(int)first_site_even; // start at 0 for even first, 1 for odd first
      int l_index = (l[NDIM-1] +size[NDIM-1])%size[NDIM-1];
      s += (l[NDIM-1] +size[NDIM-1])%size[NDIM-1];
      for (int d=NDIM-2; d>=0; d--){
        l_index = l_index*size[d] + (l[d] +size[d])%size[d];
        s += (l[d] +size[d])%size[d];
      }
      if( s%2==0 ){
        l_index /= 2;
      } else {
        l_index = evensites + l_index/2;
      }
      return l_index;
    }


    /// Return a list of neighbours for a lattice divided into a given vector size
    std::array<unsigned*,NDIRS> neighbour_list(){
      return neighbours;
    }

    /// First index in a lattice loop
    int loop_begin( parity P) {
      if(P==ODD){
        return evensites;
      } else {
        return 0;
      }
    }

    // Last index in a lattice loop
    int loop_end( parity P) {
      if(P==EVEN){
        return evensites;
      } else {
        return sites;
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


template<>
inline auto vectorized_lattice_struct<4>::coordinates(int idx){
  std::array<Vec4i,NDIM> r;
  int step=1;
  for(int d=0; d<NDIM; d++){
    int first = coordinate_list[idx][d] + min[d];
    for(int v=0; v<4; v++){
      r[d].insert(v, first + size[d]*((v/step)%split[d]));
    }
    step *= split[d];
  }
  return r;
}

template<>
inline auto vectorized_lattice_struct<8>::coordinates(int idx){
  std::array<Vec8i,NDIM> r;
  int step=1;
  for(int d=0; d<NDIM; d++){
    int first = coordinate_list[idx][d] + min[d];
    for(int v=0; v<8; v++){
      r[d].insert(v, first + size[d]*((v/step)%split[d]));
    }
    step *= split[d];
  }
  return r;
}

template<>
inline auto vectorized_lattice_struct<16>::coordinates(int idx){
  std::array<Vec16i,NDIM> r;
  int step=1;
  for(int d=0; d<NDIM; d++){
    int first = coordinate_list[idx][d] + min[d];
    for(int v=0; v<16; v++){
      r[d].insert(v, first + size[d]*((v/step)%split[d]));
    }
    step *= split[d];
  }
  return r;
}







struct backend_lattice_struct {
  lattice_struct lattice;

  void setup(lattice_struct _lattice){
    lattice = _lattice;
  }

  template< int vector_size >
  vectorized_lattice_struct<vector_size> * get_vectorized_lattice() {
    // Create one if not already created
    static bool init = true;
    static vectorized_lattice_struct<vector_size> * vlat; 
    if(init){
      vlat = new vectorized_lattice_struct<vector_size>(&lattice);
      init = false;
    }

    return vlat;
  }
};




#endif