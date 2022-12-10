/// This routine uses internal MPI (or other) commands to remap the
/// node numbers so that the communications should be optimised.
///
/// Implements two calls,
/// int remap(int i) : given "computed" node index i, returns true mpi rank
///     the input i is the logical node index, runs fastest to directions x,y,z,(t)
///     must be compatible with node_rank layout!
///
/// int inverse_remap(int i) : inverse of above
///
/// Used especially in transformation coordinate -> node_rank, in
/// lattice_struct::node_rank in


#include "plumbing/defs.h"
#include "plumbing/lattice.h"

#if defined(NODE_LAYOUT_TRIVIAL)

////////////////////////////////////////////////////////////////////
// This is the trivial version, so that the order is unmodified
////////////////////////////////////////////////////////////////////

void lattice_struct::allnodes::create_remap() {
    hila::out0 << "Node remapping: NODE_LAYOUT_TRIVIAL (no reordering)\n";
    lattice.nodes.map_array = nullptr;
    lattice.nodes.map_inverse = nullptr;
}

unsigned lattice_struct::allnodes::remap(unsigned i) const {
    return i;
}

unsigned lattice_struct::allnodes::inverse_remap(unsigned i) const {
    return i;
}


#elif defined(NODE_LAYOUT_BLOCK)

////////////////////////////////////////////////////////////////////
// Arrange nodes so that NODE_LAYOUT_BLOCK nodes are "close"
//
// For example, in 2d and NODE_LAYOUT_BLOCK = 4 the MPI indices run as
//
//  1  2 |  5  6 | 9  10
//  3  4 |  7  8 | 11 12
//  -----+-------+------
// 13 14 | 17 18 | 21 22
// 15 16 | 19 20 | 23 24
//
////////////////////////////////////////////////////////////////////

void lattice_struct::allnodes::create_remap() {

    hila::out0 << "Node remapping: NODE_LAYOUT_BLOCK with blocksize " << NODE_LAYOUT_BLOCK
            << '\n';

    // let us allow only factors of 5,3 and 2 in NODE_LAYOUT_BLOCK
    CoordinateVector blocksize;
    CoordinateVector blockdivs = lattice.nodes.n_divisions;
    int nblocks = NODE_LAYOUT_BLOCK;
    blocksize.fill(1);

    bool found = true;
    while (found && nblocks > 1) {

        found = false; // if not divisible stop trying
        foralldir (d) {
            // divide directions one by one
            int div;
            for (div = 2; div <= 5; div++)
                if (nblocks % div == 0 && blockdivs[d] % div == 0)
                    break;
            if (div <= 5) {
                blocksize[d] *= div;
                blockdivs[d] /= div;
                nblocks /= div;
                found = true;
            }
        }
    }

    hila::out0 << "Node block size " << blocksize << "  block division " << blockdivs << '\n';

    // now blocksize tells us how many nodes to each direction in block
    // these are taken in order

    lattice.nodes.map_array =
        (unsigned *)memalloc(lattice.nodes.number * sizeof(unsigned));
    // we don't need the inverse at the moment?
    // lattice.nodes.map_inverse = (unsigned *)memalloc(lattice.nodes.number *
    // sizeof(unsigned));
    lattice.nodes.map_inverse = nullptr;

    nblocks = 1;
    foralldir (d)
        nblocks *= blocksize[d];

    // Loop over the "logical" node indices (i.e.)

    for (int i = 0; i < lattice.nodes.number; i++) {
        // lcoord is the coordinate of the logical node,
        // bcoord the block coord and icoord coord inside block
        CoordinateVector lcoord, bcoord, icoord;
        int idiv = i;
        foralldir (d) {
            lcoord[d] = idiv % lattice.nodes.n_divisions[d];
            idiv /= lattice.nodes.n_divisions[d];

            bcoord[d] = lcoord[d] / blocksize[d];
            icoord[d] = lcoord[d] % blocksize[d];
        }

        // ii - index within a block
        // bi - index of a block
        int ii, bi, im, bm;
        ii = bi = 0;
        im = bm = 1;
        foralldir(d) {
            ii += icoord[d] * im;
            im *= blocksize[d];

            bi += bcoord[d] * bm;
            bm *= blockdivs[d];
        }

        // hila::out0 << lcoord << bcoord << icoord << '\n';
        // hila::out0 << "ii " << ii << " bi " << bi << '\n';

        lattice.nodes.map_array[i] = bi * nblocks + ii;
    }
}

/// And the call interface for remapping

unsigned lattice_struct::allnodes::remap(unsigned i) const {
    return lattice.nodes.map_array[i];
}

unsigned lattice_struct::allnodes::inverse_remap(unsigned idx) const {
    for (int i = 0; i < lattice.nodes.number; i++)
        if (lattice.nodes.map_array[i] == idx)
            return i;

    return 1 << 30; // big number, crash
}


#else

NODE_LAYOUT_BLOCK or NODE_LAYOUT_TRIVIAL must be defined

#endif
