#ifndef _BACKEND_LATTICE_H_
#define _BACKEND_LATTICE_H_

#include "../lattice.h"
#include "vector_types.h"

#ifndef SUBNODE_LAYOUT
static_assert(0, "SUBNODE_LAYOUT needs to be defined to use vectorized lattice");
#endif

/////////////////////////////////////////////////////////////////////////////////////////
/// Vectorized lattice struct stores information needed for vectorized access to lattice fields.
/// - loop_begin, loop_end
/// - neighbour arrays
/// - twist permutations when needed
/// - allocation size, including new halo sites (replaces field_alloc_size)
/// these are used by fields which are laid out with vector structure
template <int vector_size>
struct vectorized_lattice_struct {
    static_assert(vector_size > 0, "Vector size in vectorized_lattice_struct");

  public:
    /// pointer to the original lattice
    // lattice_struct *lattice;

    /// vector sites on this node
    size_t v_sites;
    /// subnode divisions to different directions
    CoordinateVector subdivisions;
    /// origin of the 1st subnode = origin of mynode
    CoordinateVector subnode_origin, subnode_size;

    /// True if boundary needs a permutation
    bool is_boundary_permutation[NDIM];
    /// True if the boundary elements are local
    bool only_local_boundary_copy[NDIRS];
    /// permutation vectors
    int boundary_permutation[NDIRS][vector_size];
    /// offsets to boundary halos
    unsigned halo_offset[NDIRS], halo_offset_odd[NDIRS], n_halo_vectors[NDIRS];
    /// storage for indexes to halo sites
    unsigned *RESTRICT halo_index[NDIRS];

    /// move data from receive buffer -- sending is fine as it is
    /// takes the role of nn_comms
    unsigned *recv_list[NDIRS];
    /// The size of the receive list in each direction
    unsigned recv_list_size[NDIRS];

    // coordinate offsets to nodes
    using int_vector_t = typename hila::vector_base_type<int, vector_size>::type;
    int_vector_t coordinate_offset[NDIM];
    // using coordinate_compound_vec_type = typename Vector<NDIM,
    // typename vector_base_type<int,vector_size>::type>; coordinate_compound_vec_type
    // coordinate_offset;
    CoordinateVector *RESTRICT coordinate_base;

    /// Storage for neighbour indexes on each site
    unsigned *RESTRICT neighbours[NDIRS];

    /// The storage size of a field
    size_t alloc_size;
    /// A wait array for the vectorized field
    unsigned char *RESTRICT vec_wait_arr_;

    /// Check if this is the first subnode
    bool is_on_first_subnode(CoordinateVector v) {
        v = v.mod(lattice.size());
        foralldir (d) {
            if (v[d] < subnode_origin[d] || v[d] >= subnode_origin[d] + subnode_size[d])
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
    /// Note that this maps site indices of 64 bit elements (double)and 32 bit (float) so
    /// that
    ///   vector_index(32) = vector_index(64)/2 and
    ///   index_in_vector(32) = index_in_vector(64) + (vector_index(64)%2)*vector_size(64)
    /// This removes the last direction where number of subnodes is divisible
    /////////////////////////////////////////////////////////////////////////////////////////////
    vectorized_lattice_struct() {
        /// Initialize

        /// sites on vector
        v_sites = lattice.mynode.volume() / vector_size;
        subdivisions = lattice.mynode.subnodes.divisions;
        subnode_size = lattice.mynode.subnodes.size;
        subnode_origin = lattice.mynode.min;

        // hila::out0 << "Subdivisions " << subdivisions << '\n';

        // the basic division is done using "float" vectors -
        // for "double" vectors the vector_size and number of subnodes
        // is halved to direction lattice.mynode.subnodes.lastype_divided.dir

        if (vector_size == VECTOR_SIZE / sizeof(double)) {
            subdivisions[lattice.mynode.subnodes.merged_subnodes_dir] /= 2;
            subnode_size[lattice.mynode.subnodes.merged_subnodes_dir] *= 2;
        }

        hila::out0 << "Setting up lattice struct with vector of size " << vector_size
                   << " elements\n";

        get_boundary_permutations();

        for (Direction d = (Direction)0; d < NDIRS; d++) {
            neighbours[d] = (unsigned *)memalloc(v_sites * sizeof(unsigned));
        }

        get_neighbours_and_local_halo();

        get_receive_lists();
        build_wait_arrays();

        set_coordinates();

    } // end of initialization

    /////////////////////////////////////////////////////////////////////
    /// Find the boundary permutations to different directions
    /////////////////////////////////////////////////////////////////////
    void get_boundary_permutations() {

        // boundary permutation is done in a "layout-agnostic" way:
        // go to the edge of the subnode, check the neighbour, and
        // identify the permutation as index % vector_length
        // Uses

        int step = 1, sdmul = 1;
        foralldir (d) {
            is_boundary_permutation[d] = (subdivisions[d] > 1);

            // upper corner of 1st subnode
            CoordinateVector here;
            here = subnode_origin + subnode_size;
            here.asArray() -= 1;

            unsigned idx = lattice.site_index(here);
            assert(idx % vector_size == 0); // is it really on 1st subnode

            // loop over subnodes
            for (int i = 0; i < vector_size; i++) {

                // get the site index of the neighbouring site
                CoordinateVector h = lattice.coordinates(idx + i) + d;

                // remember to mod the coordinate on lattice
                h = h.mod(lattice.size());
                int rank = lattice.node_rank(h);
                unsigned nn = lattice.site_index(h, rank);
                boundary_permutation[d][i] = nn % vector_size;
            }

            // permutation to opposite direction is the inverse, thus:
            for (int i = 0; i < vector_size; i++)
                boundary_permutation[-d][boundary_permutation[d][i]] = i;
        }
    }

    /////////////////////////////////////////////////////////////////////////
    /// Set the neighbour array, and also local halo mapping
    /// reserve extra storage for permutation halo sites
    /// there are 2*(area) v-sites to each direction with permutation,
    ///  + mpi buffers too if needed. (separately allocated)
    /// to directions without permutation there is also
    ///  - 2*(area) v-sites, filled in directly by MPI or copying from local node
    ///  if no mpi comm.  The copying is done to ease implementing different boundary
    ///  conditions
    /// These come automatically when we tally up the neigbours below
    /////////////////////////////////////////////////////////////////////////
    void get_neighbours_and_local_halo() {

        // check special case: 1st subnode is across the whole lattice to Direction d and
        // no boundary permutation
        // we do the copy also in this case, in order to implement other boundary
        // conditions Slows down a bit the periodic case, but with MPI comms this should
        // make no difference

        foralldir (d) {
            if (lattice.nodes.n_divisions[d] == 1 && !is_boundary_permutation[d]) {
                only_local_boundary_copy[d] = only_local_boundary_copy[-d] = true;
            } else {
                only_local_boundary_copy[d] = only_local_boundary_copy[-d] = false;
            }
        }

        // accumulate here points off-subnode (to halo)
        size_t c_offset = v_sites;
        for (Direction d = (Direction)0; d < NDIRS; ++d) {

            halo_offset[d] = c_offset;
            for (int i = 0; i < v_sites; i++) {
                int j = vector_size * i; // the "original lattice" index for the 1st site of vector
                CoordinateVector here = lattice.coordinates(j);
                // std::cout << here << '\n';

                if (is_on_first_subnode(here + d)) {

                    assert(lattice.neighb[d][j] % vector_size == 0); // consistency check
                    Direction ad = abs(d);

                    if (only_local_boundary_copy[d] &&
                        ((is_up_dir(d) && here[ad] == lattice.size(ad) - 1) ||
                         (is_up_dir(-d) && here[ad] == 0))) {
                        neighbours[d][i] = c_offset++;
                    } else {
                        // standard branch, within the subnode
                        neighbours[d][i] = lattice.neighb[d][j] / vector_size;
                    }
                } else {
                    neighbours[d][i] = c_offset++; // now points beyond the lattice
                }
            }
            n_halo_vectors[d] = c_offset - halo_offset[d];
            halo_offset_odd[d] = halo_offset[d] + n_halo_vectors[d] / 2;
            assert(n_halo_vectors[d] % 2 == 0);

            /// set also the index array, if needed
            /// halo_index[d] points to the subnode modded neighbour site to dir d, if
            /// there is boundary twist (meaning there are on-node subnodes to this dir)
            /// we'll use the standard neighb array to do this.
            if (n_halo_vectors[d] > 0 && is_boundary_permutation[abs(d)]) {
                halo_index[d] = (unsigned *)memalloc(n_halo_vectors[d] * sizeof(unsigned));
                int j = 0;
                for (int i = 0; i < v_sites; i++) {
                    if (neighbours[d][i] >= v_sites) {
                        // scan neighbour sites within the vector i -- some must be inside
                        // the node, which we want in this case
                        int k, n = -1;
                        bool found = false;
                        for (k = 0; k < vector_size; k++) {
                            if (lattice.neighb[d][i * vector_size + k] < lattice.mynode.sites) {
                                if (!found) {
                                    n = lattice.neighb[d][i * vector_size + k] / vector_size;
                                    found = true;
                                } else
                                    assert(n ==
                                           lattice.neighb[d][i * vector_size + k] / vector_size);
                            }
                        }
                        assert(n >= 0);
                        halo_index[d][j++] = n;
                    }
                }
                assert(j == n_halo_vectors[d]);

            } else if (only_local_boundary_copy[d]) {
                // set the local untwisted array copy here

                halo_index[d] = (unsigned *)memalloc(n_halo_vectors[d] * sizeof(unsigned));
                int j = 0;
                for (int i = 0; i < v_sites; i++) {
                    if (neighbours[d][i] >= v_sites) {
                        halo_index[d][j++] = lattice.neighb[d][i * vector_size] / vector_size;
                        assert(lattice.neighb[d][i * vector_size] % vector_size == 0);
                    }
                }

            } else
                halo_index[d] = nullptr; // no special copy here - mpi fills
        }

        /// Finally, how big the field allocation should be - IN SITES, not vectors
        alloc_size = c_offset * vector_size;

        if (alloc_size >= (1ULL << 32)) {
            report_too_large_node();
        }
    }

    //////////////////////////////////////////////////////////////////////////////
    /// Get neighbour receive indices for MPI
    //////////////////////////////////////////////////////////////////////////////
    void get_receive_lists() {

        for (Direction d = (Direction)0; d < NDIRS; d++) {
            if (is_boundary_permutation[abs(d)] && lattice.nodes.n_divisions[abs(d)] > 1) {

                // now need to receive and copy - note: now this is in terms of
                // non-vector sites.   Set the recv_list to point to where to move the
                // stuff Note: now the stuff has to be moved to boundary_halo, not to
                // lattice n!

                recv_list_size[d] = lattice.mynode.sites / lattice.mynode.size[abs(d)];
                recv_list[d] = (unsigned *)memalloc(recv_list_size[d] * sizeof(unsigned));

                int j = 0;
                for (int i = 0; i < lattice.mynode.sites; i++) {
                    if (lattice.neighb[d][i] >= lattice.mynode.sites) {

                        // i/vector_size is the "vector index" of site, and
                        // i % vector_size the index within the vector.
                        // vector neighbour is neighbours[d][i/vector_size]
                        // remember to do the permutation too

                        recv_list[d][j++] =
                            neighbours[d][i / vector_size] * vector_size + (i % vector_size);
                    }
                }
                assert(j == recv_list_size[d]);

            } else {
                // now use halo_offset directly for buffer
                recv_list[d] = nullptr;
                recv_list_size[d] = 0;
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////
    /// Build the structs for coordinates
    /////////////////////////////////////////////////////////////////////////
    void set_coordinates() {

        /// first vector_size elements should give the coordinates of vector offsets
        CoordinateVector base = lattice.coordinates(0);
        for (int i = 0; i < vector_size; i++) {
            CoordinateVector diff = lattice.coordinates(i) - base;
            foralldir (d)
                coordinate_offset[d].insert(i, diff[d]);
        }

        // and then set the coordinate_base with the original coords
        coordinate_base = (CoordinateVector *)memalloc(v_sites * sizeof(CoordinateVector));
        for (int i = 0; i < v_sites; i++) {
            coordinate_base[i] = lattice.coordinates(vector_size * i);
        }
    }

////////////////////////////////////////////////////////////////////////////
/// Finally, initialize wait arrays
/// it is a bit mask array containing a bit at location dir if the neighbour
/// at that dir is out of the local volume
    void build_wait_arrays() {
        vec_wait_arr_ = (dir_mask_t *)memalloc(v_sites * sizeof(dir_mask_t));

        for (int i = 0; i < v_sites; i++) {
            vec_wait_arr_[i] = 0; /* basic, no wait */
            foralldir (dir) {
                Direction odir = -dir;
                if (lattice.nodes.n_divisions[dir] > 1) {
                    if (neighbours[dir][i] >= v_sites)
                        vec_wait_arr_[i] = vec_wait_arr_[i] | (1 << dir);
                    if (neighbours[odir][i] >= v_sites)
                        vec_wait_arr_[i] = vec_wait_arr_[i] | (1 << odir);
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////
    /// Return the communication info
    lattice_struct::nn_comminfo_struct get_comminfo(int d) {
        return lattice.get_comminfo(d);
    }

    /////////////////////////////////////////////////////////////////////////
    /// get neighbours for this, with 2 different methods:
    /// First vector neighbour.  Now idx is the vector index
    unsigned vector_neighbour(Direction d, int idx) const {
        return neighbours[d][idx];
    }

    /// this gives the neighbour when the lattice is traversed
    /// site-by-site.  Now idx is the "true" site index, not vector index
    unsigned site_neighbour(Direction d, int idx) const {
        return vector_size * neighbours[d][idx / vector_size] + idx % vector_size;
    }

    //////////////////////////////////////////////////////////////////////////
    /// Return the number of sites that need to be allocated
    /// returns sites, not vectors!
    unsigned field_alloc_size() const {
        return alloc_size;
    }

    //////////////////////////////////////////////////////////////////////////
    /// Return the coordinates of each vector nested as
    /// coordinate[Direction][vector_index]
    auto coordinates(int idx) const {
        // std::array<typename vector_base_type<int,vector_size>::type ,NDIM> r;
        CoordinateVector_t<int_vector_t> r;
        // Vector<NDIM,int_vector_t> r;
        foralldir (d)
            r.e(d) = coordinate_offset[d] + coordinate_base[idx][d];
        return r;
    }

    auto coordinate(unsigned idx, Direction d) const {
        return coordinate_offset[d] + coordinate_base[idx][d];
    }

    // parity is the same for all elements in vector, return scalar
    ::Parity site_parity(int idx) {
        return coordinate_base[idx].parity();
    }

    /// First index in a lattice loop
    unsigned loop_begin(::Parity P) const {
        if (P == ODD) {
            return v_sites / 2;
        } else {
            return 0;
        }
    }

    /// Last index in a lattice loop
    unsigned loop_end(::Parity P) const {
        if (P == EVEN) {
            return v_sites / 2;
        } else {
            return v_sites;
        }
    }
};

///////////////////////////////////////////////////////////////////////////////
/// Helper class for loading the vectorized lattice
///////////////////////////////////////////////////////////////////////////////

struct backend_lattice_struct {

    void setup(const lattice_struct &lat) {}

    /// Returns a vectorized lattice with given vector size
    template <int vector_size>
    vectorized_lattice_struct<vector_size> *get_vectorized_lattice() {
        // Create one if not already created

        assert(lattice.id() == 0 &&
               "Vectorized lattice layout only possible with main (original) lattice");


        static vectorized_lattice_struct<vector_size> *vlat = nullptr;
        if (vlat == nullptr) {
            vlat = new vectorized_lattice_struct<vector_size>();
        }

        return vlat;
    }
};

#endif