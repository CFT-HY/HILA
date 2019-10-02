/// This routine decides how the lattice is distributed among
/// the nodes.  This is the generic version, which just
/// throws the system to as "square"" blocks as possible.
/// TODO:
/// First tries to divide the lattice so that every node has 
/// the same size, but if that fails then allow different
/// division to one direction

#include "../plumbing/globals.h"


/***************************************************************/

/* number of primes to be used in factorization */
#define NPRIMES 8
static int prime[NPRIMES] = {2,3,5,7,11,13,17,19};

/* Set up now squaresize and nsquares - arrays 
 * Print info to outf as we proceed 
 */

void lattice_struct::setup_layout( )
{
  int i,msize,dir,nfactors[NPRIMES], nodesiz[NDIM];

  output0 << "Standard lattice layout:\n " << NDIM << " dimensions\n";

  /* Figure out dimensions of hyperrectangle - this version ensures all are same size */

  if (volume % numnodes()) {
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
  
  for (i=0; i<NDIM; i++) {
    nodesiz[i] = size[i];
    nodes.ndir[i] = 1;
  }
  
  for (int n=NPRIMES-1; n>=0; n--) for(i=0; i<nfactors[n]; i++) {
    /* figure out which direction to divide -- start from the largest prime, because
     * we don't want this to be last divisor! (would probably wind up with size 1) 
     */

    // find largest divisible dimension of h-cubes - start from last, because 
    // SF and spatial FFT.
    for (msize=1,dir=0; dir<NDIM; dir++)
      if (nodesiz[dir] > msize && nodesiz[dir]%prime[n] == 0 ) msize = nodesiz[dir];

    // if one direction with largest dimension has already been
    // divided, divide it again.  Otherwise divide first direction
    // with largest dimension. 

    // Switch here to first divide along t-direction, in
    // order to 
    // a) minimize spatial blocks, for FFT
    // b) In sf t-division is cheaper (1 non-communicating slice)
    
    for (dir=NDIM-1; dir>=0; dir--)
      if (nodesiz[dir]==msize && nodes.ndir[dir]>1 && 
          nodesiz[dir]%prime[n] == 0) break;

    /* If not previously sliced, take one direction to slice */
    if (dir < 0) for (dir=NDIM-1; dir>=0; dir--)
      if( nodesiz[dir]==msize && nodesiz[dir]%prime[n] == 0) break;

    if (dir < 0) {
      /* This cannot happen */
      output0 << "CANNOT HAPPEN! in setup_layout_generic.c\n";
      finishrun();
    }

    /* Now slice it */
    nodesiz[dir] /= prime[n]; nodes.ndir[dir] *= prime[n];

  }
  
  // set up struct nodes variables
  nodes.number = numnodes();
  foralldir(dir) {
    nodes.divisors[dir].resize(nodes.ndir[dir]+1);
    // trivial, evenly spaced divisors -- note: last element == size[dir]
    for (int i=0; i<=nodes.ndir[dir]; i++) 
      nodes.divisors[dir].at(i) = i*nodesiz[dir];
  }

  #ifdef USE_MPI
  // For MPI, remap the nodes for periodic torus
  // in the desired manner 
  // we have at least 2 options:
  // map_node_layout_trivial.c
  // map_node_layout_block2.c - for 2^n n.n. blocks
  
  nodes.create_remap();
  
  if (mynode() == 0) {
    output0 << "\n Sites on node: ";
    foralldir(dir) {
      if (dir > 0) { output0 << " x "; }
      output0 << nodesiz[dir];
    }
    output0 << "\n Processor layout: ";
    foralldir(dir) {
      if (dir > 0) { output0 << " x "; }
      output0 << nodes.ndir[dir];
    }
    output0 << '\n';
  }
  #endif
}
