/// This routine decides how the lattice is distributed among
/// the nodes.  This is the generic version, which just
/// throws the system to as "square"" blocks as possible.
/// TODO:
/// First tries to divide the lattice so that every node has
/// the same size, but if that fails then allow different
/// division to one direction

#include "plumbing/defs.h"
#include "plumbing/lattice.h"

/***************************************************************/

/* number of primes to be used in factorization */
#define NPRIMES 9
static int prime[NPRIMES] = {2, 3, 5, 7, 11, 13, 17, 19, 23};

/* Set up now squaresize and nsquares - arrays
 * Print info to outf as we proceed
 */

void lattice_struct::setup_layout() {
    int i, msize, dir, nfactors[NPRIMES], nodesiz[NDIM];

    hila::out0 << "Standard lattice layout:\n " << NDIM << " dimensions\n";

    /* Figure out dimensions of hyperrectangle - this version ensures all are same size */

    if (l_volume % hila::number_of_nodes()) {
        hila::out0 << " No hope of laying out the lattice using " << hila::number_of_nodes()
                << " nodes\n";
        hila::finishrun();
    }

    /* Factorize the node number in primes
     * These factors must be used in slicing the lattice!
     */
    i = hila::number_of_nodes();
    for (int n = 0; n < NPRIMES; n++) {
        nfactors[n] = 0;
        while (i % prime[n] == 0) {
            i /= prime[n];
            nfactors[n]++;
        }
    }
    if (i != 1) {
        hila::out0 << " Cannot factorize " << hila::number_of_nodes() << " nodes with primes up to "
                << prime[NPRIMES - 1] << '\n';
        hila::finishrun();
    }

    for (i = 0; i < NDIM; i++) {
        nodesiz[i] = size(i);
        nodes.n_divisions[i] = 1;
    }

    for (int n = NPRIMES - 1; n >= 0; n--)
        for (i = 0; i < nfactors[n]; i++) {
            /* figure out which direction to divide -- start from the largest prime,
             * because we don't want this to be last divisor! (would probably wind up with
             * size 1)
             */

            // find largest divisible dimension of h-cubes - start from last, because
            // SF and spatial FFT.
            for (msize = 1, dir = 0; dir < NDIM; dir++)
                if (nodesiz[dir] > msize && nodesiz[dir] % prime[n] == 0)
                    msize = nodesiz[dir];

            // if one direction with largest dimension has already been
            // divided, divide it again.  Otherwise divide first direction
            // with largest dimension.

            // Switch here to first divide along t-direction, in
            // order to
            // a) minimize spatial blocks, for FFT
            // b) In sf t-division is cheaper (1 non-communicating slice)

            for (dir = NDIM - 1; dir >= 0; dir--)
                if (nodesiz[dir] == msize && nodes.n_divisions[dir] > 1 &&
                    nodesiz[dir] % prime[n] == 0)
                    break;

            /* If not previously sliced, take one direction to slice */
            if (dir < 0)
                for (dir = NDIM - 1; dir >= 0; dir--)
                    if (nodesiz[dir] == msize && nodesiz[dir] % prime[n] == 0)
                        break;

            if (dir < 0) {
                /* This cannot happen */
                hila::out0 << "CANNOT HAPPEN! in setup_layout_simple\n";
                hila::finishrun();
            }

            /* Now slice it */
            nodesiz[dir] /= prime[n];
            nodes.n_divisions[dir] *= prime[n];
        }

    // set up struct nodes variables
    nodes.number = hila::number_of_nodes();
    foralldir(dir) {
        nodes.divisors[dir].resize(nodes.n_divisions[dir] + 1);
        // trivial, evenly spaced divisors -- note: last element == size(dir)
        for (int i = 0; i <= nodes.n_divisions[dir]; i++)
            nodes.divisors[dir].at(i) = i * nodesiz[dir];
    }


    // For MPI, remap the nodes for periodic torus
    // in the desired manner
    // we have at least 2 options:
    // map_node_layout_trivial.c
    // map_node_layout_block2.c - for 2^n n.n. blocks

    nodes.create_remap();

    if (hila::myrank() == 0) {
        hila::out0 << "\n Sites on node: ";
        foralldir(dir) {
            if (dir > 0) {
                hila::out0 << " x ";
            }
            hila::out0 << nodesiz[dir];
        }
        hila::out0 << "\n Processor layout: ";
        foralldir(dir) {
            if (dir > 0) {
                hila::out0 << " x ";
            }
            hila::out0 << nodes.n_divisions[dir];
        }
        hila::out0 << '\n';
    }
}
