/* This routine decides how the lattice is distributed among
 * the nodes.  This is the generic version, which just
 * throws the system to as square blocks as possible
 */

#include "../plumbing/lattice.h"


/***************************************************************/

/* number of primes to be used in factorization */
#define NPRIMES 8
static int prime[NPRIMES] = {2,3,5,7,11,13,17,19};

/* Set up now squaresize and nsquares - arrays 
 * Print info to outf as we proceed 
 */

void setup_layout( int siz[NDIM], 
		   int nsquares[NDIM], 
		   int squaresize[NDIM],
		   int *map_node_list )
{
  int i,j,dir,nfactors[NPRIMES];

  output0 << "LATTICE LAYOUT (GENERIC):\n " << NDIM << " dimensions, layout options: ";
  #ifdef EVENFIRST
  output0 << "EVENFIRST ";
  #endif
  output0 << '\n';

  /* Figure out dimensions of hyperrectangle - this version ensures all are same size */

  foralldir(dir) {
    nsquares[dir] = 1;
    squaresize[dir] = lattice.size[dir];
  }

  if (lattice.volume % numnodes()) {
    output0 << " No hope of laying out the lattice using " << numnodes() << " nodes\n";
    finishrun();
  }

  /* Factorize the node number in primes
   * These factors must be used in slicing the lattice!
   */
  i = numnodes(); 
  for (int n=0; n<NPRIMES; n++) {
    nfactors[n] = 0;
    while (i%prime[n] == 0) { i/=prime[n]; nfactors[n]++; }
  }
  if (i != 1) {
    output0 << " Cannot factorize " << numnodes() << " nodes with primes up to " 
            << prime[NPRIMES-1] << '\n';
    finishrun();
  }
  
  for (int n=NPRIMES-1; n>=0; n--) for(i=0; i<nfactors[n]; i++) {
    /* figure out which direction to divide -- start from the largest prime, because
     * we don't want this to be last divisor! (would probably wind up with size 1) 
     */

    /* find largest divisible dimension of h-cubes 
     *   - start the division from direction 0, because
     * these are likely to be within the same node!  */
    for(j=1,dir=0; dir<NDIM; dir++)
      if( squaresize[dir]>j && squaresize[dir]%prime[n] == 0 ) j=squaresize[dir];
    /* if one direction with largest dimension has already been
       divided, divide it again.  Otherwise divide first direction
       with largest dimension. */

    /* NEW-Switch here to first divide along t-direction, in
     * order to 
     * a) minimize spatial blocks, for FFT
     * b) In sf t-division is cheaper (1 non-communicating slice)
     */

    for (dir=NDIM-1; dir>=0; dir--)
      if( squaresize[dir]==j && nsquares[dir]>1 && 
	  squaresize[dir]%prime[n] == 0) break;

    /* If not previously sliced, take one direction to slice */
    if (dir < 0) for (dir=NDIM-1; dir>=0; dir--)
      if( squaresize[dir]==j && squaresize[dir]%prime[n] == 0) break;

    if (dir < 0) {
      /* This cannot happen */
      output0 << "CANNOT HAPPEN! in setup_layout_generic.c\n";
      finishrun();
    }

    /* Now slice it */
    squaresize[dir] /= prime[n]; nsquares[dir] *= prime[n];

  }
  
  if (mynode() == 0) {
    output0 << "\n Sites on node: ";
    foralldir(dir) {
      if (dir > 0) output0 << " x ";
      output0 << squaresize[dir];
    }
    output0 << "\n Processor layout: ";
    foralldir(dir) {
      if (dir > 0) output0 << " x ";
      output0 << nsquares[dir];
    }
    output0 << '\n';
  }

#ifdef USE_MPI
  /* Then, for MPI, remap the nodes for periodic torus
   * in the desired manner 
   * we have at least 2 options:
   * map_node_layout_trivial.c
   * map_node_layout_block2.c - for 2^n n.n. blocks
   */
  map_node_layout( nsquares, squaresize, map_node_list );
#endif

}
