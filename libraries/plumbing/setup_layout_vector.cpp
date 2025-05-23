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

// Set up now squaresize and nsquares - arrays
// Print info to outf as we proceed

void lattice_struct::setup_layout() {
    int nfactors[NPRIMES];
    CoordinateVector nodesiz;

    hila::print_dashed_line();
    hila::out0 << "LAYOUT: subnode lattice, with " << VECTOR_SIZE / sizeof(float)
            << " subnodes\n"
            << "Enabling vector length " << VECTOR_SIZE * 8
            << " bits = " << VECTOR_SIZE / sizeof(double) << " doubles or "
            << VECTOR_SIZE / sizeof(float) << " floats/ints\n";
    hila::out0 << "Lattice size ";
    foralldir (d) {
        if (d != 0)
            hila::out0 << " x ";
        hila::out0 << size(d);
    }
    hila::out0 << "  =  " << l_volume << " sites\n";
    hila::out0 << "Dividing to " << hila::number_of_nodes() << " nodes\n";
    hila::out0 << "Layout using vector of " << number_of_subnodes << " elements\n";

    foralldir (d)
        if (size(d) % 2 != 0) {
            hila::out0 << "Lattice must be even to all directions (odd size:TODO)\n";
            hila::finishrun();
        }

    // we want to divide up to numnode * vector_size virtual nodes
    // use the float vector size to divide to hila::number_of_nodes() * vector_size
    // nodes, this is the most demanding. The directions where the extra divisions have
    // been done the node size must be even, so that these can be handled by vectors

    // Factorize the node number in primes
    // These factors must be used in slicing the lattice!

    // number of virtual nodes
    int nn = hila::number_of_nodes() * number_of_subnodes;

    int i = nn;
    for (int n = 0; n < NPRIMES; n++) {
        nfactors[n] = 0;
        while (i % prime[n] == 0) {
            i /= prime[n];
            nfactors[n]++;
        }
    }
    if (i != 1) {
        hila::out0 << "Cannot factorize " << hila::number_of_nodes()
                << " nodes with primes up to " << prime[NPRIMES - 1] << '\n';
        hila::finishrun();
    }

    // strategy: try to increase the box size to one of the directions until rem = 0
    // find the optimal direction to do it
    // Use simple heuristic: take the dim with the least amount of added "ghost sites"

    CoordinateVector nsize;
    int64_t ghosts[NDIM];
    foralldir (d) {
        int64_t cosize = l_volume / size(d);
        int64_t n = size(d);
        while ((n * cosize) % nn != 0)
            n++; // virtual size can be odd
        // now nsize is the new would-be size
        ghosts[d] = (n - size(d)) * cosize;
        nsize[d] = n;
    }

    int64_t ghost_volume = 1;
    foralldir (d)
        ghost_volume *= nsize[d];

    // if the division goes even nsize = size() and ghost_volume = volume

    // now try to divide the nsize-volume to subnodes.  We don't try to do
    // the division to direction where there are ghost sites

    CoordinateVector divisions, subdiv;

    int gdir;
    bool secondtime = false;
    do {
        // try the division a couple of times, if the 1st fails

        if (ghost_volume > l_volume) {
            gdir = NDIM - 1;
            foralldir (j)
                if (ghosts[gdir] > ghosts[j])
                    gdir = j;
            // gdir is the direction where we do uneven division (if done)
            // hila::out0 << "gdir " << gdir << " ghosts gdir " << ghosts[gdir] << " nsize
            // "
            // << nsize[gdir] << '\n';
        } else
            gdir = -1;

        foralldir (i) {
            nodesiz[i] =
                (i == gdir) ? nsize[i] : size(i); // start with ghosted lattice size
            divisions[i] = 1;
        }

        for (int n = NPRIMES - 1; n >= 0; n--)
            for (i = 0; i < nfactors[n]; i++) {
                // figure out which direction to divide -- start from the largest prime,
                // because we don't want this to be last divisor! (would probably wind
                // up with size 1)

                // Try to keep even node sizes, these are needed for vectors
                // We don't worry about evenness to gdir

                // find largest divisible dimension of h-cubes
                int msize = 1;
                int dir, mdir;
                for (dir = 0; dir < NDIM; dir++) {
                    if (nodesiz[dir] > msize &&
                        ((dir == gdir && nodesiz[dir] % prime[n] == 0) ||
                         (dir != gdir && nodesiz[dir] % (2 * prime[n]) == 0)) &&
                        nodesiz[dir] / prime[n] > 3) {
                        msize = nodesiz[dir];
                        mdir = dir;
                    }
                }

                // even divide failed, divide to odd
                if (msize == 1) {
                    for (dir = 0; dir < NDIM; dir++) {
                        if (nodesiz[dir] > msize && dir != gdir &&
                            nodesiz[dir] % (prime[n]) == 0) {
                            msize = nodesiz[dir];
                            mdir = dir;
                        }
                    }
                }

                if (msize == 1) {
                    // This cannot happen
                    hila::out0 << "CANNOT HAPPEN! in setup_layout_vector.c\n";
                    hila::finishrun();
                }

                // Now slice it
                nodesiz[mdir] /= prime[n];
                divisions[mdir] *= prime[n];

                //  hila::out0 << nodesiz << ' ' << divisions << '\n';
            }

        // Division done, now check that the div makes sense

        bool fail = false;

        foralldir (dir)
            if (nodesiz[dir] < 3)
                fail = true; // don't allow nodes of size 1 or 2

        if (!fail) {

            // check here that this can be used for vectorized division

            subdiv.fill(1);
            bool div_done;
            int n_subn = 1;
            do {
                div_done = false;
                foralldir (dir) {
                    // the direction where the vector subnodes are must not be
                    // an uneven direction, node size to the direction should be
                    // divisible by 2 and the number of nodes to this dir should also be
                    // a multiple of subdivs

                    int sd = subdiv[dir] * 2;
                    if (dir != gdir && nodesiz[dir] % 2 == 0 &&
                        divisions[dir] % sd == 0 && n_subn < number_of_subnodes) {
                        subdiv[dir] = sd;
                        n_subn *= 2;
                        div_done = true;
                        mynode.subnodes.merged_subnodes_dir = dir;
                    }
                }
            } while (div_done && n_subn < number_of_subnodes);

            if (n_subn != number_of_subnodes)
                fail = true;
        }

        if (fail && !secondtime && gdir >= 0) {
            secondtime = true;
            ghosts[gdir] =
                (1ULL
                 << 62); // this short-circuits direction gdir, some other taken next
        } else if (fail) {
            hila::out0 << "Could not successfully lay out the lattice with "
                    << hila::number_of_nodes() << " nodes!\n";
            hila::out0 << "  The division of ";
            foralldir (d) {
                hila::out0 << lattice.size(d);
                if (d < NDIM - 1)
                    hila::out0 << '*';
            }
            hila::out0 << " lattice using " << hila::number_of_nodes()
                    << " nodes with vector layout can be done\n";
            hila::out0 << "  if the lattice can be divided into ";
            hila::out0 << hila::number_of_nodes() << '*' << number_of_subnodes
                    << " virtual nodes so that the virtual node size is\n";
            hila::out0 << "  even to directions where the extra divisions are done, "
                       "and the node size is > 2.\n";

            hila::finishrun();
        }

    } while (secondtime);

    // set up struct nodes variables
    nodes.number = hila::number_of_nodes();
    foralldir (dir) {
        nodesiz[dir] *= subdiv[dir];
        nodes.n_divisions[dir] = divisions[dir] / subdiv[dir];
        nodes.divisors[dir].resize(nodes.n_divisions[dir] + 1);
        // Node divisors: note, this MUST BE compatible with
        // node_rank in lattice.cpp
        // to be sure, we use naively the same method than in node_rank
        // last element will be size(dir), for convenience
        int n = -1;
        for (int i = 0; i <= size(dir); i++)
            if ((i * nodes.n_divisions[dir]) / size(dir) != n) {
                ++n;
                nodes.divisors[dir][n] = i;
            }
    }

    // set up the subnode divisions here -- rest is set in setup_node
    foralldir (d)
        mynode.subnodes.divisions[d] = subdiv[d];

    // mynode is set up in setup_node

    // Now division done - check how good it is
    int ghost_slices = 0;
    if (gdir >= 0) {
        ghost_slices = nsize[gdir] - size(gdir);

        hila::out0 << "\nUsing uneven node division to direction " << gdir << ":\n";
        hila::out0 << "Lengths: " << nodes.n_divisions[gdir] - ghost_slices << " * ("
                << nodesiz[gdir] << " sites) + " << ghost_slices << " * ("
                << nodesiz[gdir] - 1 << " sites)\n";
        hila::out0 << "Divisions: ";
        for (int i = 0; i < nodes.n_divisions[gdir]; i++) {
            if (i > 0)
                hila::out0 << " - ";
            hila::out0 << nodes.divisors[gdir][i + 1] - nodes.divisors[gdir][i];
        }
        hila::out0 << "\nFilling efficiency: " << (100.0 * size(gdir)) / nsize[gdir]
                << "%\n";

        if (ghost_slices > nodes.n_divisions[gdir] / 2)
            hila::out0 << "NOTE: number of smaller nodes > large nodes \n";
    }

    // this was hila::number_of_nodes() > 1
    if (1) {
        hila::out0 << "\nSites on node: ";
        foralldir (dir) {
            if (dir > 0)
                hila::out0 << " x ";
            if (dir == gdir)
                hila::out0 << '(' << nodesiz[dir] - 1 << '-' << nodesiz[dir] << ')';
            else
                hila::out0 << nodesiz[dir];
        }
        int64_t ns = 1;
        foralldir (dir)
            ns *= nodesiz[dir];
        if (ghost_slices > 0) {
            int64_t ns2 = ns * (nodesiz[gdir] - 1) / nodesiz[gdir];
            hila::out0 << "  =  " << ns2 << " - " << ns << '\n';
        } else {
            hila::out0 << "  =  " << ns << '\n';
        }

        hila::out0 << "Node layout: ";
        foralldir (dir) {
            if (dir > 0)
                hila::out0 << " x ";
            hila::out0 << nodes.n_divisions[dir];
        }
        hila::out0 << "  =  " << hila::number_of_nodes() << " nodes\n";

#ifdef VECTORIZED

        hila::out0 << "Node subdivision to 32bit elems: ";
        foralldir (dir) {
            if (dir > 0)
                hila::out0 << " x ";
            hila::out0 << subdiv[dir];
        }
        hila::out0 << "  =  " << number_of_subnodes << " subnodes\n";

        hila::out0 << "Sites on subnodes: ";
        foralldir (dir) {
            if (dir > 0)
                hila::out0 << " x ";
            if (dir == gdir)
                hila::out0 << '(' << nodesiz[dir] - 1 << '-' << nodesiz[dir] << ')';
            else
                hila::out0 << nodesiz[dir] / subdiv[dir];
        }
        hila::out0 << '\n';

        Direction dmerge = mynode.subnodes.merged_subnodes_dir;

        hila::out0 << "Node subdivision to 64bit elems: ";
        foralldir (dir) {
            if (dir > 0)
                hila::out0 << " x ";
            hila::out0 << ((dir == dmerge) ? subdiv[dir] / 2 : subdiv[dir]);
        }
        hila::out0 << "  =  " << number_of_subnodes / 2 << " subnodes\n";

        hila::out0 << "Sites on subnodes: ";
        foralldir (dir) {
            if (dir > 0)
                hila::out0 << " x ";
            if (dir == gdir)
                hila::out0 << '(' << nodesiz[dir] - 1 << '-' << nodesiz[dir] << ')';
            else
                hila::out0 << ((dir == dmerge) ? 2 * nodesiz[dir] / subdiv[dir]
                                            : nodesiz[dir] / subdiv[dir]);
        }
        hila::out0 << '\n';

#endif
    }

    // For MPI, remap the nodes for periodic torus
    // in the desired manner
    // we have at least 2 options:
    // map_node_layout_trivial.c
    // map_node_layout_block2.c - for 2^n n.n. blocks

    nodes.create_remap();


    hila::print_dashed_line();
}
