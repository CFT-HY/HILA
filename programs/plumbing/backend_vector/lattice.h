#ifndef _BACKEND_LATTICE_H_
#define _BACKEND_LATTICE_H_


/// Splits the local lattice into equal sections for vectorization
struct vectorized_lattice_struct  {
  public:
    std::array<int,NDIM> size;
    std::array<int,NDIM> split;
    int sites, evensites, oddsites, alloc_size;
    unsigned vector_size;
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

    struct halo_site {
      unsigned nb_index;
      unsigned halo_index;
      parity par;
      direction dir;
    };
    std::vector<halo_site> halo_sites;


    /// Set up the vectorized lattice:
    /// * Split the lattice and record size and splits
    /// * Map indeces into local coordinate_list
    /// * Set up neighbour vector references
    vectorized_lattice_struct(lattice_struct * _lattice, int _vector_size) {
      // Initialize
      lattice =  _lattice;
      vector_size = _vector_size;
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

      // Setup neighbour array
      int halo_index = 0;
      for(int d=0; d<NDIRS; d++){
        neighbours[d] = (unsigned *) malloc(sizeof(unsigned)*sites);
        for(unsigned i = 0; i<sites; i++){
          coordinate_vector l = coordinate_list[i];
          coordinate_vector nb = l;
          int k;
          if (is_up_dir(d)) {
            k = d;
            nb[d] = l[d] + 1;
          } else {
            k = opp_dir(d);
            nb[k] = l[k] - 1;
          }
          if( nb[k] >= 0 && nb[k] < size[k] ) {
            neighbours[d][i] = get_index(nb);
          } else {
            // This is outside this split lattice (maybe partly outside the node)
            neighbours[d][i] = sites + halo_index;
            // Find the corresponding site on the other side of the lattice
            halo_site hs;
            hs.nb_index = get_index(nb);
            hs.halo_index = halo_index;
            hs.dir = (direction)d;
            if( i < evensites ){
              hs.par = EVEN;
            } else {
              hs.par = ODD;
            }
            halo_sites.push_back(hs);
            halo_index++;
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
      
      //printf(" Vectorized lattice size: (%d %d %d %d)\n",
      //  size[0], size[1], size[2], size[3]);
      //printf(" Vectorized lattice split: (%d %d %d %d)\n",
      //  split[0], split[1], split[2], split[3]);
      //for(int d=0; d<NDIRS; d++){
      //  printf(" permutation %d: (",(int)d);
      //  for(int v=0;v<vector_size; v++){
      //    printf("%d, ", boundary_permutation[d][v]);
      //  }
      //  printf(")\n");
      //}

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
    std::array<Vec16i,NDIM> coordinates_Vec16i(int idx){
      assert(vector_size == 16);
      std::array<Vec16i,NDIM> r;
      int step=1;
      for(int d=0; d<NDIM; d++){
        int first = coordinate_list[idx][d] + min[d];
        for(int v=0; v<vector_size; v++){
          r[d].insert(v, first + size[d]*((v/step)%split[d]));
        }
        step *= split[d];
      }
      return r;
    }
    std::array<Vec8i,NDIM> coordinates_Vec16f(int idx){
      return coordinates_Vec8i(idx);
    }

    std::array<Vec8i,NDIM> coordinates_Vec8i(int idx){
      assert(vector_size == 8);
      std::array<Vec8i,NDIM> r;
      int step=1;
      for(int d=0; d<NDIM; d++){
        int first = coordinate_list[idx][d] + min[d];
        for(int v=0; v<vector_size; v++){
          r[d].insert(v, first + size[d]*((v/step)%split[d]));
        }
        step *= split[d];
      }
      return r;
    }
    std::array<Vec8i,NDIM> coordinates_Vec8f(int idx){
      return coordinates_Vec8i(idx);
    }
    std::array<Vec8i,NDIM> coordinates_Vec8d(int idx){
      return coordinates_Vec8i(idx);
    }

    std::array<Vec4i,NDIM> coordinates_Vec4d(int idx){
      assert(vector_size == 4);
      std::array<Vec4i,NDIM> r;
      int step=1;
      for(int d=0; d<NDIM; d++){
        int first = coordinate_list[idx][d] + min[d];
        for(int v=0; v<vector_size; v++){
          r[d].insert(v, first + size[d]*((v/step)%split[d]));
        }
        step *= split[d];
      }
      return r;
    }


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
    void reduce_node_sum(T & value, bool distribute){
      lattice->reduce_node_sum(value, distribute);
    };

    template <typename T>
    void reduce_node_product(T & value, bool distribute){
      lattice->reduce_node_product(value, distribute);
    };

};



struct backend_lattice_struct {
  std::vector<vectorized_lattice_struct*> vectorized_lattices;
  lattice_struct lattice;

  void setup(lattice_struct _lattice){
    lattice = _lattice;
  }

  vectorized_lattice_struct * get_vectorized_lattice(int vector_size) {
    // Check if the vectorized lattice has been created
    for( vectorized_lattice_struct * vl : vectorized_lattices ) {
      if( vl->vector_size == vector_size )
        return vl;
    }

    // Not found, setup here
    vectorized_lattice_struct * vectorized_lattice = new vectorized_lattice_struct(&lattice, vector_size);
    vectorized_lattices.push_back(vectorized_lattice);
    return vectorized_lattice;
  }
};




#endif