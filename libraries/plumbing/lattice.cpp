
#include "plumbing/defs.h"
#include "plumbing/lattice.h"
#include "plumbing/field.h"

// Reporting on possibly too large node: stop if
// node size (with buffers) is larger than 2^31 - too close for comfort!

void report_too_large_node() {
    if (hila::myrank() == 0) {
        hila::out << "Node size too large: size = " << lattice->mynode.size[0];
        for (int d = 1; d < NDIM; d++)
            hila::out << " x " << lattice->mynode.size[d];
        hila::out << " + communication buffers = " << lattice->mynode.field_alloc_size;
        hila::out << "\nConsider using more nodes (smaller node size).\n";
        hila::out << "[TODO: allow 64bit index?]\n";
    }
    hila::finishrun();
}


// Define the global lattice handle
Lattice lattice;
/// A list of all defined lattices (for the future expansion)

// Keep track of defined lattices
std::vector<lattice_struct *> defined_lattices;

// global bookkeeping
namespace hila {
int64_t n_gather_avoided = 0, n_gather_done = 0;
}

/// General lattice setup
void lattice_struct::setup_base_lattice(const CoordinateVector &siz) {

    assert(defined_lattices.size() == 0 && "lattice.setup() can be called only once");

    l_label = 0;
    defined_lattices.push_back(this);
    parent = nullptr; // this has no parent lattice

    l_volume = 1;
    foralldir (i) {
        l_size[i] = siz[i];
        l_volume *= siz[i];
    }

    setup_layout();

    setup_nodes();

    // set up the comm arrays
    create_std_gathers();

    // Initialize wait_array structures - has to be after std gathers()
    initialize_wait_arrays();

#ifdef SPECIAL_BOUNDARY_CONDITIONS
    // do this after std. boundary is done
    init_special_boundaries();
#endif

    // Alignment: set field_alloc_size to be divisible by 256
    if (mynode.field_alloc_size % 256 > 0)
        mynode.field_alloc_size += 256 - mynode.field_alloc_size % 256;

#ifndef VANILLA
    /* Setup backend-specific lattice info if necessary */
    backend_lattice = new backend_lattice_struct;
    backend_lattice->setup(*this);
#endif

    if (hila::check_input) {
        hila::out << "***** Input check done *****\n";
        hila::finishrun();
    }

    test_std_gathers();
}

///////////////////////////////////////////////////////////////////////
/// The routines is_on_mynode(), node_rank(), site_index()
/// implement the "layout" of the nodes and sites of the lattice.
/// To be changed in different implementations!
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
/// Get the node rank for coordinates
/// This is the fundamental routine which defines how the nodes
/// are mapped.  map_node_layout MUST BE compatible with this
/// algorithm!  So if you change this, change that too
/// Divisor moved now below
///
/// Here the node number along one Direction is calculated with
///    loc * n_divisions/l_size  (with integer division)
/// example: with l_size=14 and n_divisions=3,
/// the dividers are at 0, 5, 10, 14
///
///////////////////////////////////////////////////////////////////////

int lattice_struct::node_rank(const CoordinateVector &loc) const {

    int i = (loc[NDIM - 1] * nodes.n_divisions[NDIM - 1]) / l_size[NDIM - 1];

    for (int dir = NDIM - 2; dir >= 0; dir--) {
        i = i * nodes.n_divisions[dir] + ((loc[dir] * nodes.n_divisions[dir]) / l_size[dir]);
    }

    /* do we want to remap this?  YES PLEASE */
    if (i < 0 || i >= hila::number_of_nodes())
        hila::out << "NODE NUMBER " << i << " CV " << loc << " MYNODE " << hila::myrank() << '\n';
    i = nodes.remap(i);

    return (i);
}

////////////////////////////////////////////////////////////////////////
/// In this routine we set up the node divisors
/// this must be compatible with the node_rank above
////////////////////////////////////////////////////////////////////////

void lattice_struct::setup_node_divisors() {

    foralldir (dir) {
        nodes.divisors[dir].resize(nodes.n_divisions[dir] + 1);
        // to be sure, we use naively the same method than in node_rank
        // last element will be size(dir), for convenience
        int n = -1;
        for (int i = 0; i <= l_size[dir]; i++) {
            while (n < (i * nodes.n_divisions[dir]) / l_size[dir]) {
                ++n;
                nodes.divisors[dir][n] = i;
            }
        }
        // hila::out0 << "Divisors ";
        // for (int i=0;i<nodes.n_divisions[dir]; i++) hila::out0 << nodes.divisors[dir][i]
        // << " "; hila::out0 << '\n';
    }
}


///////////////////////////////////////////////////////////////////////
/// Is the coordinate on THIS node ?
///////////////////////////////////////////////////////////////////////

bool lattice_struct::is_on_mynode(const CoordinateVector &loc) const {
    int d;

    for (int dir = 0; dir < NDIM; dir++) {
        d = loc[dir] - mynode.min[dir];
        if (d < 0 || d >= mynode.size[dir])
            return false;
    }
    return true;
}

///////////////////////////////////////////////////////////////////////
/// give site index for ON NODE sites
/// Note: loc really has to be on this node
/// Separate versions for standard layout, boundary layer layout, and
/// for vector layout.
///////////////////////////////////////////////////////////////////////

#ifndef SUBNODE_LAYOUT

#ifndef GPU_OVERLAP_COMM

unsigned lattice_struct::site_index(const CoordinateVector &loc) const {


    unsigned i = mynode.get_logical_index(loc);

    // now i contains the `running index' for site
#if defined(EVEN_SITES_FIRST)
    if (loc.parity() == EVEN)
        return (i / 2); /* even site index */
    else
        return (i / 2 + mynode.evensites); /* odd site */
#else
    return (i);
#endif
}

#else // Now GPU_VERLAP_COMM layout version

unsigned lattice_struct::site_index(const CoordinateVector &loc) const {

    return mynode.map_site_index[mynode.get_logical_index(loc)];
}

#endif


///////////////////////////////////////////////////////////////////////
/// give site index for nodeid sites
/// compare to above
///////////////////////////////////////////////////////////////////////

// unsigned lattice_struct::site_index(const CoordinateVector &loc, const unsigned nodeid) const {
//     const node_info ni = nodes.nodeinfo(nodeid);

//     unsigned i = loc[NDIM - 1] - ni.min[NDIM - 1];
//     int s = loc[NDIM - 1];
//     for (int dir = NDIM - 2; dir >= 0; dir--) {
//         int l = loc[dir] - ni.min[dir];
//         i = i * ni.size[dir] + l;
//         s += loc[dir];
//     }

//     // now i contains the `running index' for site
// #if defined(EVEN_SITES_FIRST)
//     if (s % 2 == 0)
//         return (i / 2); /* even site index */
//     else
//         return (i / 2 + ni.evensites); /* odd site */
// #else
//     return (i);
// #endif
// }


#else // now SUBNODE_LAYOUT

///////////////////////////////////////////////////////////////////////
/// The (AVX) vectorized version of the site_index function.
/// Now there are two stages: an "inner" vect index, which goes over
/// over "virtual nodes", and the "outer" index inside the virtual node.
/// This is to enable the (float/32bit) and the (double/64bit) vectorized
/// structures, where the latter is achieve by merging two of the 32bit subnodes,
/// along the Direction merged_subnodes_dir
/// E.g.
/// a 2-dim. 4x8 node is divided into 8 / 4 subnodes as follows:
/// here 0-3 / 0-7 is the index within subnode, and a-h / a-d the subnode label.
///
///     32bit(float)   64bit(double)
///                                      32bit storage: 64bit storage:
///    0a 1a | 0b 1b   0a 2a | 0b 2b     0(abcdefgh)    0(abcd) 1(abcd)
///    2a 3a | 2b 3b   4a 6a | 4b 6b     1(abcdefgh)    2(abcd) 3(abcd)
///    -------------                     2(abcdefgh)    4(abcd) 5(abcd)
///    0e 1e | 0f 1f   1a 3a | 1b 3b     3(abcdefgh)    6(abcd) 7(abcd)
///    2e 3e | 2f 3f   5a 7a | 5b 7b
///    -------------   -------------    32bit vec 1st half <-> even 64bit
///    0c 1c | 0d 1d   0c 2c | 0d 2d    32bit 2nd half <-> odd 64bit
///    2c 3c | 2d 3d   4c 6c | 4d 6d
///    -------------
///    0g 1g | 0h 1h   1c 3c | 1d 3d
///    2g 3g | 2h 3h   5c 7c | 5d 7d
///
/// The "storage" order above is the site-by-site order, enabling site
/// traversal with maximum locality. It also enables mixed 32bit/64bit vector
/// operations with half vectors.
///
/// Direction where the "doubling" is done is the last Direction where subnodes
/// are divided. In layout, this will become the "slowest" Direction
///
///////////////////////////////////////////////////////////////////////

unsigned lattice_struct::site_index(const CoordinateVector &loc) const {
    return site_index(loc, hila::myrank());
}

///////////////////////////////////////////////////////////////////////
/// give site index for nodeid sites
/// compare to above
///////////////////////////////////////////////////////////////////////

unsigned lattice_struct::site_index(const CoordinateVector &loc, const unsigned nodeid) const {
    int dir, l, s, subl;
    unsigned i;

    assert(nodeid < nodes.number);
    const node_info ni = nodes.nodeinfo(nodeid);

    // let's mod the coordinate to partition
    // get subnode size - divisons are the same in all nodes

    // foralldir(d) assert( ni.size[d] % mynode.subnodes.divisions[d] == 0);

    CoordinateVector subsize;
    subsize.asArray() = ni.size.asArray() / mynode.subnodes.divisions.asArray();

    dir = mynode.subnodes.merged_subnodes_dir;
    l = loc[dir] - ni.min[dir];
    subl = l / subsize[dir];
    // flip to interleaved ordering, see above  --  this should work
    // automatically w. mapping between 32 and 64 bit vector elems
    subl = subl / 2 + (subl % 2) * (mynode.subnodes.divisions[dir] / 2);

    i = 0;
    s = 0;

    for (dir = NDIM - 1; dir >= 0; dir--) {
        l = loc[dir] - ni.min[dir];
        if (dir != mynode.subnodes.merged_subnodes_dir) {
            subl = subl * mynode.subnodes.divisions[dir] + l / subsize[dir];
        }
        i = i * subsize[dir] + l % subsize[dir];
        s += loc[dir];
    }

#if defined(EVEN_SITES_FIRST)
    if (s % 2 == 0)
        i = i / 2; /* even site index */
    else
        i = i / 2 + ni.evensites / hila::number_of_subnodes; /* odd site */
#endif

    return (subl + hila::number_of_subnodes * i);
}

#endif // SUBNODE_LAYOUT


/////////////////////////////////////////////////////////////////////
/// Define nodeinfo(nodeid)
/// Gives the coordinate range managed by the node number nodeid
/////////////////////////////////////////////////////////////////////

node_info lattice_struct::allnodes::nodeinfo(int nodeid) const {
    // get the logical node id
    assert(nodeid < number);
    nodeid = inverse_remap(nodeid);

    int64_t nodevol = 1;
    node_info ninfo;

    foralldir (d) {
        int nodecoord = nodeid % n_divisions[d];
        nodeid = nodeid / n_divisions[d];

        ninfo.min[d] = divisors[d].at(nodecoord);
        ninfo.size[d] = divisors[d].at(nodecoord + 1) - ninfo.min[d];
        nodevol *= ninfo.size[d];
    }

    if (nodevol >= (1ULL << 32)) {
        // node size is larger than 2^32-1 - does not fit to largest unsigned int!
        report_too_large_node();
    }

    if (nodevol % 2 == 0) {
        ninfo.evensites = ninfo.oddsites = nodevol / 2;
    } else {
        if (ninfo.min.parity() == EVEN) {
            ninfo.evensites = nodevol / 2 + 1;
            ninfo.oddsites = nodevol / 2;
        } else {
            ninfo.evensites = nodevol / 2;
            ninfo.oddsites = nodevol / 2 + 1;
        }
    }

    return ninfo;
}


/////////////////////////////////////////////////////////////////////
/// only thing setup_nodes now does is to call mynode.setup and
/// set the max size field of nodes-struct
/////////////////////////////////////////////////////////////////////

void lattice_struct::setup_nodes() {

    nodes.number = hila::number_of_nodes();
    // Loop over all node "origins"

    // n keeps track of the node "root coordinates"
    CoordinateVector n(0);

    nodes.max_size = 0;

    // use nodes.divisors - vectors to fill in stuff
    for (int i = 0; i < nodes.number; i++) {
        auto ni = nodes.nodeinfo(i);
        foralldir (d) {
            if (ni.size[d] > nodes.max_size[d])
                nodes.max_size[d] = ni.size[d];
        }

        // use the opportunity to set up mynode when it is met
        if (i == hila::myrank())
            mynode.setup(ni, *this);
    }
}

////////////////////////////////////////////////////////////////////////
/// Obtain the "logical index" on this node from CoordinateVector c
/// returns c_x + c_y * nx + c_z * nx * ny + ...
////////////////////////////////////////////////////////////////////////

unsigned lattice_struct::node_struct::get_logical_index(const CoordinateVector &l) const {

    unsigned res = 0;
    foralldir (d) {
        res += static_cast<unsigned>(l[d] - min[d]) * size_factor[d];
    }
    return res;
}

////////////////////////////////////////////////////////////////////////
/// Helper function to scan trough the coordinates of the local node
/// advances fastest x, then y etc. coordinates
/// To be called inside loop up to node_struct.size times
////////////////////////////////////////////////////////////////////////

void lattice_struct::node_struct::advance_local_coordinate(CoordinateVector &l) const {
    foralldir (d) {
        if (++l[d] < (min[d] + size[d]))
            return;
        l[d] = min[d];
    }
}

////////////////////////////////////////////////////////////////////////

#ifdef GPU_OVERLAP_COMM

int lattice_struct::loop_ranges(Parity par, bool gather_on, hila::iter_range_t &ranges) const {

    // can we do it in one sweep?
    if (par == Parity::all && !gather_on) {
        ranges.min[0] = 0;
        ranges.max[0] = mynode.volume;
        return 1;
    }

    if (par == Parity::even) {
        ranges.min[0] = 0;
        ranges.max[0] = mynode.inner_start_odd;
        ranges.min[1] = mynode.boundary_start_even;
        ranges.max[1] = mynode.boundary_start_odd;
    } else if (par == Parity::odd) {
        ranges.min[0] = mynode.inner_start_odd;
        ranges.max[0] = mynode.boundary_start_even;
        ranges.min[1] = mynode.boundary_start_odd;
        ranges.max[1] = mynode.volume;
    } else {
        ranges.min[0] = 0;
        ranges.max[0] = mynode.boundary_start_even;
        ranges.min[1] = mynode.boundary_start_even;
        ranges.max[1] = mynode.volume;
    }
    return 2;
}

static bool is_boundary_site(const CoordinateVector &cv, const CoordinateVector &min,
                             const CoordinateVector &size) {
    bool boundary_site = false;
    foralldir (d) {
        boundary_site = (cv[d] == min[d] || cv[d] == (min[d] + size[d] - 1));
        if (boundary_site)
            break;
    }
    return boundary_site;
}

////////////////////////////////////////////////////////////////////////
/// @internal construct the mapping from logical to real index,
/// separating boundary sites from "bulk"
///
////////////////////////////////////////////////////////////////////////

void lattice_struct::node_struct::construct_index_map() {

    // alloc resize map
    map_site_index.resize(volume);

    // scan the sites, accumulate as as classified by
    // inner and boundary sites. For EVEN_SITES_FIRST
    // also by parity

    CoordinateVector l = min;

    unsigned inner_even = 0, inner_odd = 0, boundary_even = 0, boundary_odd = 0;
    // scan through the sites first just to count
    for (size_t i = 0; i < volume; ++i) {
        bool boundary_site = is_boundary_site(l, min, size);

        if (l.parity() == EVEN) {
            if (boundary_site)
                boundary_even++;
            else
                inner_even++;
        } else {
            if (boundary_site)
                boundary_odd++;
            else
                inner_odd++;
        }

        advance_local_coordinate(l);
    }

    inner_start_odd = inner_even;
    boundary_start_even = inner_start_odd + inner_odd;
    boundary_start_odd = boundary_start_even + boundary_even;
    inner_vol = inner_even + inner_odd;
    boundary_vol = boundary_even + boundary_odd;

    inner_even = 0;
    inner_odd = inner_start_odd;
    boundary_even = boundary_start_even;
    boundary_odd = boundary_start_odd;

    l = min;
    for (size_t i = 0; i < volume; ++i) {
        bool boundary_site = is_boundary_site(l, min, size);

        // auto index = get_logical_index(l);
        // i is the logical index

#if defined(EVEN_SITES_FIRST)
        if (l.parity() == EVEN) {
            if (boundary_site)
                map_site_index[i] = boundary_even++;
            else
                map_site_index[i] = inner_even++;
        } else {
            if (boundary_site)
                map_site_index[i] = boundary_odd++;
            else
                map_site_index[i] = inner_odd++;
        }
#else
        if (boundary_site)
            map_site_index[i] = boundary_even++;
        else
            map_site_index[i] = inner_even++;
#endif

        advance_local_coordinate(l);
    }

    assert(boundary_odd == volume);
}

#endif

////////////////////////////////////////////////////////////////////////
/// Fill in mynode fields -- node_rank() must be set up OK
////////////////////////////////////////////////////////////////////////
void lattice_struct::node_struct::setup(node_info &ni, lattice_struct &lattice) {

    rank = hila::myrank();

    min = ni.min;
    size = ni.size;

    evensites = ni.evensites;
    oddsites = ni.oddsites;
    volume = ni.evensites + ni.oddsites;

    // first_site_even = (min.parity() == EVEN);
    parent = &lattice;

    // set up the auxiliary size_factor array
    // Now l.dot(size_factor) gives logical running index
    unsigned v = 1;
    foralldir (d) {
        size_factor[d] = v; // = size[d-1] * size[d-2] * ..
        v *= size[d];
    }

#ifdef GPU_OVERLAP_COMM
    construct_index_map();
#endif

    // map site indexes to locations -- coordinates array
    // after the above site_index should work

#ifdef EVEN_SITES_FIRST
    coordinates.resize(volume);
    CoordinateVector l = min;
    for (unsigned i = 0; i < volume; i++) {
        coordinates[lattice.site_index(l)] = l;
        advance_local_coordinate(l);
    }

#endif

#ifdef SUBNODE_LAYOUT

    // set up the subnodes
    subnodes.setup(*this);
#endif
}

#ifdef SUBNODE_LAYOUT

////////////////////////////////////////////////////////////////////////
/// Fill in subnodes -struct
////////////////////////////////////////////////////////////////////////
void lattice_struct::node_struct::subnode_struct::setup(const node_struct &tn) {
    size.asArray() = tn.size.asArray() / divisions.asArray();
    evensites = tn.evensites / hila::number_of_subnodes;
    oddsites = tn.oddsites / hila::number_of_subnodes;
    sites = evensites + oddsites;
}

#endif

bool lattice_struct::is_this_odd_boundary(Direction d) const {
    // return false unless lattice size is odd to direction d
    // and mynode is on the boundary of the lattice

    return (lattice.size(d) % 2 > 0 && lattice->mynode.is_on_edge(d));
}

/////////////////////////////////////////////////////////////////////
/// Create the neighbour index arrays
/// This is for the index array neighbours
/// TODO: implement some other neighbour schemas!
/////////////////////////////////////////////////////////////////////

// unsigned lattice_struct::count_off_node_sites(const CoordinateVector &offset) {

//     CoordinateVector ln, l;
//     unsigned count;

//     for (int i = 0; i < mynode.volume; i++) {
//         l = coordinates(i);
//         // set ln to be the neighbour of the site
//         // TODO: FIXED BOUNDARY CONDITIONS DO NOT WRAP
//         ln = (l + offset).mod(size());

//         if (!is_on_mynode(ln))
//             count++;
//     }
//     return count;
// }


void lattice_struct::create_std_gathers() {

    // allocate neighbour arrays - TODO: these should
    // be allocated on "device" memory too!

    for (int d = 0; d < NDIRS; d++) {
        neighb[d] = (unsigned *)memalloc(((size_t)mynode.volume) * sizeof(unsigned));
    }

    size_t c_offset = mynode.volume; // current offset in field-arrays

    // We set the communication and the neigbour-array here
    int too_large_node = 0;

    for (Direction d = e_x; d < NDIRS; ++d) {

        nn_comminfo[d].index = neighb[d]; // this is not really used for nn gathers

        comm_node_struct &from_node = nn_comminfo[d].from_node;
        // we can do the opposite send during another pass of the sites.
        // This is just the gather inverted
        // NOTE: this is not the send to Direction d, but to -d!
        comm_node_struct &to_node = nn_comminfo[-d].to_node;

        from_node.init();
        to_node.init();

        // if there are no communications the rank is left as is

        // first pass through the sites
        // - set the neighb array, anything to communicate?

        for (size_t i = 0; i < mynode.volume; i++) {
            CoordinateVector ln, l;
            l = coordinates(i);
            // set ln to be the neighbour of the site
            ln = (l + d).mod(l_size);

            if (is_on_mynode(ln)) {
                neighb[d][i] = site_index(ln);
            } else {
                // short-circuit neighbour array, for later handling
                neighb[d][i] = mynode.volume;

                // Now site is off-node, this leads to gathering
                if (from_node.rank == mynode.rank) {
                    from_node.rank = to_node.rank = node_rank(ln);
                }

                if (l.parity() == EVEN)
                    from_node.evensites++;
                else
                    from_node.oddsites++;

                // evensites / oddsites classified by the parity of receiving site
                if (ln.parity() == EVEN)
                    to_node.evensites++;
                else
                    to_node.oddsites++;
            }
        }

        // store here the buffer index
        from_node.buffer = c_offset;

        from_node.sites = from_node.evensites + from_node.oddsites;
        to_node.sites = to_node.evensites + to_node.oddsites;

        assert(from_node.sites == to_node.sites);

        if (to_node.sites > 0) {

            // sitelist tells us which sites to send
            to_node.sitelist = (unsigned *)memalloc(to_node.sites * sizeof(unsigned));

            // set the remaining neighbour array indices and sitelists in another go
            // over sites. temp counters NOTE: ordering is automatically right: with a
            // given parity, neighbour node indices come in ascending order of host node
            // index - no sorting needed
            size_t to_even = 0, to_odd = 0, from_even = 0, from_odd = 0;

            for (size_t i = 0; i < mynode.volume; i++) {
                if (neighb[d][i] == mynode.volume) {
                    CoordinateVector l, ln;
                    l = coordinates(i);

                    if (l.parity() == EVEN) {
                        neighb[d][i] = c_offset + from_even;
                        ++from_even;

                    } else {
                        neighb[d][i] = c_offset + from_node.evensites + from_odd;
                        ++from_odd;
                    }

                    ln = (l + d).mod(l_size);
                    if (ln.parity() == EVEN) {
                        to_node.sitelist[to_even++] = i;
                    } else {
                        to_node.sitelist[(to_odd++) + to_node.evensites] = i;
                    }
                }
            }

#if defined(CUDA) || defined(HIP)
            auto p = copy_array_to_gpu(to_node.sitelist, to_node.sites);
            free(to_node.sitelist);
            to_node.sitelist = p;
#endif
        }

        // and advance offset for next direction
        c_offset += from_node.sites;

        if (c_offset >= (1ULL << 32))
            too_large_node = 1;

    } /* directions */

    /* Finally, set the site to the final offset (better be right!) */
    mynode.field_alloc_size = c_offset;

    if (hila::reduce_node_sum(too_large_node) > 0) {
        report_too_large_node();
    }
}


/************************************************************************/

/* this formats the wait_array, used by forallsites_waitA()site_neighbour
 * wait_array[i] contains a bit at position 1<<dir if nb(dir,i) is out
 * of lattice.
 * Site OK if ((wait_arr ^ xor_mask ) & and_mask) == 0
 */

static_assert(NDIM <= 4 && "Dimensions at most 4 in dir_mask_t = unsigned char!  Use "
                           "larger type to circumvent");

void lattice_struct::initialize_wait_arrays() {

#if !(defined(CUDA) || defined(HIP))

    /* Allocate here the mask array needed for forallsites_wait
     * This will contain a bit at location dir if the neighbour
     * at that dir is out of the local volume
     */

    wait_arr_ = (dir_mask_t *)memalloc(mynode.volume * sizeof(unsigned char));

    for (size_t i = 0; i < mynode.volume; i++) {
        wait_arr_[i] = 0; /* basic, no wait */
        foralldir (dir) {
            Direction odir = -dir;
            if (neighb[dir][i] >= mynode.volume)
                wait_arr_[i] = wait_arr_[i] | (1 << dir);
            if (neighb[odir][i] >= mynode.volume)
                wait_arr_[i] = wait_arr_[i] | (1 << odir);
        }
    }
#endif
}


#ifdef SPECIAL_BOUNDARY_CONDITIONS

/////////////////////////////////////////////////////////////////////
/// set up special boundary
/// sets up the bool array which tells if special neighbour indices
/// are needed.  Note that this is not uniform across the nodes,
/// not all nodes have them.
/////////////////////////////////////////////////////////////////////

void lattice_struct::init_special_boundaries() {
    for (Direction d = (Direction)0; d < NDIRS; ++d) {

        // default values, nothing interesting happens
        special_boundaries[d].n_even = special_boundaries[d].n_odd = special_boundaries[d].n_total =
            0;
        special_boundaries[d].is_needed = false;

        Direction od = -d;
        int coord = -1;
        // do we get up/down boundary?
        if (is_up_dir(d) && mynode.min[d] + mynode.size[d] == l_size[d])
            coord = l_size[d] - 1;
        if (is_up_dir(od) && mynode.min[od] == 0)
            coord = 0;

        if (coord >= 0) {
            // now we got it

            if (nodes.n_divisions[abs(d)] == 1) {
                special_boundaries[d].is_needed = true;
                special_boundaries[d].offset = mynode.field_alloc_size;

                for (unsigned i = 0; i < mynode.volume; i++)
                    if (coordinate(i, abs(d)) == coord) {
                        // set buffer indices
                        special_boundaries[d].n_total++;
                        if (site_parity(i) == EVEN)
                            special_boundaries[d].n_even++;
                        else
                            special_boundaries[d].n_odd++;
                    }
                mynode.field_alloc_size += special_boundaries[d].n_total;
            }
        }

        // hila::out << "Node " << hila::myrank() << " dir " << d << " min " <<
        // mynode.min << " is_on_edge "
        //   << lattice->mynode.is_on_edge(d) << '\n';

        // allocate neighbours only on 1st use, otherwise unneeded
        special_boundaries[d].neighbours = nullptr;
    }

    int toolarge = 0;
    if (mynode.field_alloc_size >= (1ULL << 32))
        toolarge = 1;
    if (hila::reduce_node_sum(toolarge) > 0) {
        report_too_large_node();
    }
}

/////////////////////////////////////////////////////////////////////
/// give the neighbour array pointer.  Allocate if needed

const unsigned *lattice_struct::get_neighbour_array(Direction d, hila::bc bc) {

#ifndef SPECIAL_BOUNDARY_CONDITIONS
    assert(bc == hila::bc::PERIODIC &&
           "non-periodic BC only if SPECIAL_BOUNDARY_CONDITIONS defined");
    return neighb[d];
#else

    // regular bc exit, should happen almost always
    if (special_boundaries[d].is_needed == false || bc == hila::bc::PERIODIC)
        return neighb[d];

    if (special_boundaries[d].neighbours == nullptr) {
        setup_special_boundary_array(d);
    }
    return special_boundaries[d].neighbours;

#endif
}

//////////////////////////////////////////////////////////////////////
/// and set up the boundary to one Direction
//////////////////////////////////////////////////////////////////////

void lattice_struct::setup_special_boundary_array(Direction d) {
    // if it is not needed or already done...
    if (special_boundaries[d].is_needed == false || special_boundaries[d].neighbours != nullptr)
        return;

    // now allocate neighbour array and the gathering array
    special_boundaries[d].neighbours = (unsigned *)memalloc(sizeof(unsigned) * mynode.volume);
    special_boundaries[d].move_index =
        (unsigned *)memalloc(sizeof(unsigned) * special_boundaries[d].n_total);

    int coord;
    int offs = special_boundaries[d].offset;
    if (is_up_dir(d))
        coord = l_size[d] - 1;
    else
        coord = 0;

    int k = 0;
    for (int i = 0; i < mynode.volume; i++) {
        if (coordinate(i, abs(d)) != coord) {
            special_boundaries[d].neighbours[i] = neighb[d][i];
        } else {
            special_boundaries[d].neighbours[i] = offs++;
            special_boundaries[d].move_index[k++] = neighb[d][i];
        }
    }

    assert(k == special_boundaries[d].n_total);
}

#endif // SPECIAL_BOUNDARY_CONDITIONS

/**
 * @internal set lattice global variables, useful for GPUs
 */

void lattice_struct::set_lattice_globals() const {
#if defined(CUDA) || defined(HIP)
    this->backend_lattice->set_device_globals(*this);
#endif
}


/**
 * @internal implementation of lattice.can_block()
 */

bool lattice_struct::can_block_by_factor(const CoordinateVector &blocking_factor) const {

#ifdef SUBNODE_LAYOUT
    return false; // blocking not implemented for vector layout!
#else

    bool ok = true;
    foralldir (d) {
        if (blocking_factor[d] <= 0) {
            hila::out0 << "ERROR: invalid blocking factor " << blocking_factor[d] << '\n';
            hila::terminate(0);
        }
        ok = ok && (l_size[d] % blocking_factor[d] == 0);
    }
    return ok;

#endif
}

/**
 * @internal implementation of lattice.block()
 */

lattice_struct *lattice_struct::block_by_factor(const CoordinateVector &blocking_factor) {
    if (!can_block_by_factor(blocking_factor)) {
        hila::out0 << "Cannot block lattice with factor " << blocking_factor << '\n';
        hila::terminate(0);
    }

    CoordinateVector blockvol;
    foralldir (d)
        blockvol[d] = l_size[d] / blocking_factor[d];

    hila::synchronize();

    for (auto *l : defined_lattices) {
        if (l->l_size == blockvol) {
            lattice.switch_to(l);
            return l;
        }
    }

    // Now did not find a lattice, make one

    int i = defined_lattices.size(); // label for new lattice

    auto *lp = new lattice_struct;

    defined_lattices.push_back(lp);
    lp->setup_blocked_lattice(blockvol, i, *this);

    return lp;
}

/**
 * @internal create blocked lattice
 */

void lattice_struct::setup_blocked_lattice(const CoordinateVector &siz, int label,
                                           lattice_struct &orig) {

    lattice.set_lattice_pointer(this);

    l_label = label;

    l_volume = 1;
    foralldir (i) {
        l_size[i] = siz[i];
        l_volume *= siz[i];
    }

    parent = &orig;

    mpi_comm_lat = orig.mpi_comm_lat;
#ifdef GPU_CCL
    gccl_comm_lat = orig.gccl_comm_lat;
#endif

    // set the layout by hand from orig lattice
    nodes.n_divisions = orig.nodes.n_divisions;
    nodes.number = orig.nodes.number;
    setup_node_divisors();

    nodes.map_array = orig.nodes.map_array;
    nodes.map_inverse = orig.nodes.map_inverse;

    setup_nodes();
    create_std_gathers();

    initialize_wait_arrays();

#ifdef SPECIAL_BOUNDARY_CONDITIONS
    // do this after std. boundary is done
    init_special_boundaries();
#endif

    // Alignment: set field_alloc_size to be divisible by 256
    if (mynode.field_alloc_size % 256 > 0)
        mynode.field_alloc_size += 256 - mynode.field_alloc_size % 256;

#ifndef VANILLA
    /* Setup backend-specific lattice info if necessary */
    backend_lattice = new backend_lattice_struct;
    backend_lattice->setup(*this);
#endif

    test_std_gathers();
}


/////////////////////////////////////////////////////////////////////
/// Create the neighbour index arrays
/// This is for the index array neighbours
/// TODO: implement some other neighbour schemas!
/////////////////////////////////////////////////////////////////////

#if 0

/// This is a helper routine, returning a vector of comm_node_structs for all nodes
/// involved with communication.
/// If receive == true, this is "receive" end and index will be filled.
/// For receive == false the is "send" half is done.

std::vector<lattice_struct::comm_node_struct>
lattice_struct::create_comm_node_vector(CoordinateVector offset, unsigned *index,
                                        bool receive) {

    // for send flip the offset
    if (!receive)
        offset = -offset;

    // temp work array: np = node points
    std::vector<unsigned> np_even(nodes.number); // std::vector initializes to zero
    std::vector<unsigned> np_odd(nodes.number);

    // we'll go through the sites twice, in order to first resolve the size of the
    // needed buffers, then to fill them.  Trying to avoid memory fragmentation a bit

    // pass over sites
    int num = 0; // number of sites off node
    for (unsigned i = 0; i < mynode.volume; i++) {
        CoordinateVector ln, l;
        l = coordinates(i);
        ln = (l + offset).mod(size());

        if (is_on_mynode(ln)) {
            if (receive)
                index[i] = site_index(ln);
        } else {
            // Now site is off-node, this will leads to gathering
            // use ci.index to store the node rank
            unsigned r = node_rank(ln);

            if (receive) {
                index[i] = mynode.volume + r;

                // using parity of THIS
                if (l.parity() == EVEN)
                    np_even[r]++;
                else
                    np_odd[r]++;

            } else {
                // on the sending side - we use parity of target
                if (ln.parity() == EVEN)
                    np_even[r]++;
                else
                    np_odd[r]++;
            }

            num++;
        }
    }

    // count the number of nodes taking part
    unsigned nnodes = 0;
    for (int r = 0; r < nodes.number; r++) {
        if (np_even[r] > 0 || np_odd[r] > 0)
            nnodes++;
    }

    // allocate the vector
    std::vector<comm_node_struct> node_v(nnodes);

    int n = 0;
    int c_buffer = 0;
    for (int r = 0; r < nodes.number; r++) {
        if (np_even[r] > 0 || np_odd[r] > 0) {
            // add the rank
            node_v[n].rank = r;
            node_v[n].evensites = np_even[r];
            node_v[n].oddsites = np_odd[r];
            node_v[n].sites = np_even[r] + np_odd[r];

            // pre-allocate the sitelist for sufficient size
            if (!receive)
                node_v[n].sitelist =
                    (unsigned *)memalloc(node_v[n].sites * sizeof(unsigned));

            node_v[n].buffer =
                c_buffer; // running idx to comm buffer - used from receive
            c_buffer += node_v[n].sites;
            n++;
        }
    }

    // ci.receive_buf_size = c_buffer;  // total buf size

    // we'll reuse np_even and np_odd as counting arrays below
    for (int i = 0; i < nnodes; i++)
        np_even[i] = np_odd[i] = 0;

    if (!receive) {
        // sending end -- create sitelists

        for (unsigned i = 0; i < mynode.volume; i++) {
            CoordinateVector ln, l;
            l = coordinates(i);
            ln = (l + offset).mod(size());

            if (!is_on_mynode(ln)) {
                unsigned r = node_rank(ln);
                int n = 0;
                // find the node from the list
                while (node_v[n].rank != r)
                    n++;
                // we'll fill the buffers according to the parity of receieving node
                // first even, then odd sites in the buffer
                unsigned k;
                if (ln.parity() == EVEN)
                    k = np_even[n]++;
                else
                    k = node_v[n].evensites + np_odd[n]++;

                // and set the ptr to the site to be communicated
                node_v[n].sitelist[k] = i;
            }
        }

    } else {
        // receive end
        // fill in the index pointers

        for (unsigned i = 0; i < mynode.volume; i++) {
            if (index[i] >= mynode.volume) {
                int r = index[i] - mynode.volume;
                int n = 0;
                // find the node which sends this
                while (node_v[n].rank != r)
                    n++;

                CoordinateVector l = coordinates(i);
                if (l.parity() == EVEN)
                    index[i] = node_v[n].buffer + (np_even[n]++);
                else
                    index[i] = node_v[n].buffer + node_v[n].evensites + (np_odd[n]++);
            }
        }
    }

    return node_v;
}

lattice_struct::gen_comminfo_struct
lattice_struct::create_general_gather(const CoordinateVector &offset) {
    // allocate neighbour arrays - TODO: these should
    // be allocated on "device" memory too!

    gen_comminfo_struct ci;

    // communication buffer
    ci.index = (unsigned *)memalloc(mynode.volume * sizeof(unsigned));

    ci.from_node =
        create_comm_node_vector(offset, ci.index, true);          // create receive end
    ci.to_node = create_comm_node_vector(offset, nullptr, false); // create sending end

    // set the total receive buffer size from the last vector
    const comm_node_struct &r = ci.from_node[ci.from_node.size() - 1];
    ci.receive_buf_size = r.buffer + r.sites;

    return ci;
}

#endif
