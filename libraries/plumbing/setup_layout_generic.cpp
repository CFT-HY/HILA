/// Setup layout does the node division.  This version
/// first tries an even distribution, with equally sized
/// nodes, and if that fails allows slightly different
/// node sizes.

#include "plumbing/defs.h"
#include "plumbing/lattice.h"

/***************************************************************/

/* number of primes to be used in factorization */
#define NPRIMES 12
const static int prime[NPRIMES] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};

/* Set up now squaresize and nsquares - arrays
 * Print info to outf as we proceed
 */

void lattice_struct::setup_layout() {
    int nfactors[NPRIMES];
    CoordinateVector nodesiz;

    hila::print_dashed_line();
    hila::out0 << "LAYOUT: lattice size  ";
    foralldir(d) {
        if (d != 0)
            hila::out0 << " x ";
        hila::out0 << l_size[d];
    }
    hila::out0 << "  =  " << l_volume << " sites\n";
    hila::out0 << "Dividing to " << hila::number_of_nodes() << " nodes\n";

    // foralldir(d) if (l_size[d] % 2 != 0) {
    //     hila::out0 << "Lattice must be even to all directions (odd size:TODO)\n";
    //     hila::finishrun();
    // }

    // Factorize the node number in primes
    // These factors must be used in slicing the lattice!
    int nn = hila::number_of_nodes();

    int i = nn;
    for (int n = 0; n < NPRIMES; n++) {
        nfactors[n] = 0;
        while (i % prime[n] == 0) {
            i /= prime[n];
            nfactors[n]++;
        }
    }
    if (i != 1) {
        hila::out0 << "Cannot factorize " << nn << " nodes with primes up to " << prime[NPRIMES - 1]
                   << '\n';
        hila::finishrun();
    }

    // strategy: try to increase the box size to one of the directions until rem = 0
    // find the optimal direction to do it
    // Use simple heuristic: take the dim with the least amount of added "ghost sites"

    CoordinateVector nsize;
    int64_t ghosts[NDIM]; // int is too small

    foralldir(d) {
        int64_t cosize = l_volume / l_size[d];
        int64_t n = l_size[d];
        while ((n * cosize) % nn != 0)
            n++; // virtual size can be odd
        // now nsize is the new would-be size
        ghosts[d] = (n - l_size[d]) * cosize;
        nsize[d] = n;
    }

    int mdir = 0;
    bool secondtime = false;
    do {
        // try the division a couple of times, if the 1st fails

        foralldir(j) if (ghosts[mdir] > ghosts[j]) mdir = j;
        // mdir is the direction where we do uneven division (if done)
        // hila::out0 << "MDIR " << mdir << " ghosts mdir " << ghosts[mdir] << " nsize " <<
        // nsize[mdir] << '\n';

        foralldir(i) {
            nodesiz[i] = (i == mdir) ? nsize[i] : l_size[i]; // start with ghosted lattice size
            nodes.n_divisions[i] = 1;
        }

        for (int n = NPRIMES - 1; n >= 0; n--)
            for (int i = 0; i < nfactors[n]; i++) {
                // figure out which direction to divide -- start from the largest prime,
                // because we don't want this to be last divisor! (would probably wind up
                // with size 1)

                // find largest divisible dimension of h-cubes - start from last, because
                int msize = 1, dir;
                for (dir = 0; dir < NDIM; dir++)
                    if (nodesiz[dir] > msize && nodesiz[dir] % prime[n] == 0)
                        msize = nodesiz[dir];

                // if one direction with largest dimension has already been
                // divided, divide it again.  Otherwise divide first direction
                // with largest dimension.

                // Switch here to first divide along t-direction, in
                // order to
                // a) minimize spatial blocks
                // b) In sf t-division is cheaper (1 non-communicating slice)

                for (dir = NDIM - 1; dir >= 0; dir--)
                    if (nodesiz[dir] == msize && nodes.n_divisions[dir] > 1 &&
                        nodesiz[dir] % prime[n] == 0)
                        break;

                // If not previously sliced, take one direction to slice
                if (dir < 0)
                    for (dir = NDIM - 1; dir >= 0; dir--)
                        if (nodesiz[dir] == msize && nodesiz[dir] % prime[n] == 0)
                            break;

                if (dir < 0) {
                    // This cannot happen
                    hila::out0 << "CANNOT HAPPEN! in setup_layout_generic.c\n";
                    hila::finishrun();
                }

                // Now slice it
                nodesiz[dir] /= prime[n];
                nodes.n_divisions[dir] *= prime[n];
            }

        // now check that the div makes sens
        bool fail = false;
        foralldir(dir) if (nodesiz[dir] < 2) fail = true; // don't allow nodes of size 1
        if (fail && !secondtime) {
            secondtime = true;
            ghosts[mdir] =
                (1ULL << 62); // this short-circuits direction mdir, some other taken next
        } else if (fail) {
            hila::out0 << "Could not successfully lay out the lattice with "
                       << hila::number_of_nodes() << " nodes\n";
            hila::finishrun();
        }

    } while (secondtime);

    // set up struct nodes variables
    nodes.number = hila::number_of_nodes();
    setup_node_divisors();

    // Now division done - check how good it is
    int ghost_slices = nsize[mdir] - l_size[mdir];
    if (ghost_slices > 0) {
        hila::out0 << "\nUsing uneven node division to direction " << mdir << ":\n";
        hila::out0 << "Lengths: " << nodes.n_divisions[mdir] - ghost_slices << " * ("
                   << nodesiz[mdir] << " sites) + " << ghost_slices << " * (" << nodesiz[mdir] - 1
                   << " sites)\n";
        hila::out0 << "Divisions: ";
        for (int i = 0; i < nodes.n_divisions[mdir]; i++) {
            if (i > 0)
                hila::out0 << " - ";
            hila::out0 << nodes.divisors[mdir][i + 1] - nodes.divisors[mdir][i];
        }
        hila::out0 << "\nFilling efficiency: " << (100.0 * l_size[mdir]) / nsize[mdir] << "%\n";

        if (ghost_slices > nodes.n_divisions[mdir] / 2)
            hila::out0 << "NOTE: number of smaller nodes > large nodes \n";
    }

    // this was hila::number_of_nodes() > 1
    if (1) {
        hila::out0 << "\nSites on node: ";
        foralldir(dir) {
            if (dir > 0)
                hila::out0 << " x ";
            if (dir == mdir && ghost_slices > 0)
                hila::out0 << '(' << nodesiz[dir] - 1 << '-' << nodesiz[dir] << ')';
            else
                hila::out0 << nodesiz[dir];
        }
        int ns = 1;
        foralldir(dir) ns *= nodesiz[dir];
        if (ghost_slices > 0) {
            int ns2 = ns * (nodesiz[mdir] - 1) / nodesiz[mdir];
            hila::out0 << "  =  " << ns2 << " - " << ns << '\n';
        } else {
            hila::out0 << "  =  " << ns << '\n';
        }

        hila::out0 << "Processor layout: ";
        foralldir(dir) {
            if (dir > 0)
                hila::out0 << " x ";
            hila::out0 << nodes.n_divisions[dir];
        }
        hila::out0 << "  =  " << hila::number_of_nodes() << " nodes\n";
    }

    // For MPI, remap the nodes for periodic torus
    // in the desired manner
    // we have at least 2 options:
    // map_node_layout_trivial.c
    // map_node_layout_block2.c - for 2^n n.n. blocks

    nodes.create_remap();


    hila::print_dashed_line();
}
