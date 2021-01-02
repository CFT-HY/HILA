/// Setup layout does the node division.  This version
/// first tries an even distribution, with equally sized
/// nodes, and if that fails allows slightly different
/// node sizes.  

#include "plumbing/globals.h"

/***************************************************************/

/* number of primes to be used in factorization */
#define NPRIMES 12
const static int prime[NPRIMES] = {2,3,5,7,11,13,17,19,23,29,31,37};

/* Set up now squaresize and nsquares - arrays 
 * Print info to outf as we proceed 
 */

void lattice_struct::setup_layout()
{
  int nfactors[NPRIMES];
  coordinate_vector nodesiz;

  print_dashed_line();
  output0 << "LAYOUT: lattice size  ";
  foralldir(d) {
    if (d != 0) output0 << " x ";
    output0 << size(d);
  }
  output0 << "  =  " << l_volume << " sites\n";
  output0 << "Dividing to " << numnodes() << " nodes\n";

  //Factorize the node number in primes
  // These factors must be used in slicing the lattice!
  int nn = numnodes();

  int i = nn; 
  for (int n=0; n<NPRIMES; n++) {
    nfactors[n] = 0;
    while (i%prime[n] == 0) { i/=prime[n]; nfactors[n]++; }
  }
  if (i != 1) {
    output0 << "Cannot factorize " << nn << " nodes with primes up to " 
            << prime[NPRIMES-1] << '\n';
    hila::finishrun();
  }

  int remainder = l_volume % nn;   // remainder = 0 even division

  // strategy: try to increase the box size to one of the directions until rem = 0
  // find the optimal direction to do it
  // Use simple heuristic: take the dim with the least amount of added "ghost sites"

  coordinate_vector ghosts, nsize;
  foralldir(d) {
    int cosize = l_volume / size(d);
    int n = size(d);
    while ( (n*cosize) % nn != 0 ) n++;  // virtual size can be odd
    // now nsize is the new would-be size 
    ghosts[d] = (n-size(d))*cosize;
    nsize[d] = n;
  }

  int mdir = 0;
  bool secondtime = false;
  do {
    // try the division a couple of times, if the 1st fails

    foralldir(j) if (ghosts[mdir] > ghosts[j]) mdir = j;
    // mdir is the direction where we do uneven division (if done)
    // output0 << "MDIR " << mdir << " ghosts mdir " << ghosts[mdir] << " nsize " << nsize[mdir] << '\n';

    foralldir(i) {
      nodesiz[i] = (i == mdir) ? nsize[i] : size(i);   // start with ghosted lattice size
      nodes.n_divisions[i] = 1;
    }
    
    for (int n=NPRIMES-1; n>=0; n--) for (int i=0; i<nfactors[n]; i++) {
      // figure out which direction to divide -- start from the largest prime, because
      // we don't want this to be last divisor! (would probably wind up with size 1) 
        
      // find largest divisible dimension of h-cubes - start from last, because 
      int msize = 1, dir;
      for (dir=0; dir<NDIM; dir++)
        if (nodesiz[dir] > msize && nodesiz[dir]%prime[n] == 0 ) msize = nodesiz[dir];

      // if one direction with largest dimension has already been
      // divided, divide it again.  Otherwise divide first direction
      // with largest dimension. 

      // Switch here to first divide along t-direction, in
      // order to 
      // a) minimize spatial blocks
      // b) In sf t-division is cheaper (1 non-communicating slice)
      
      for (dir=NDIM-1; dir>=0; dir--)
        if (nodesiz[dir]==msize && nodes.n_divisions[dir]>1 && 
            nodesiz[dir]%prime[n] == 0) break;

      // If not previously sliced, take one direction to slice
      if (dir < 0) for (dir=NDIM-1; dir>=0; dir--)
        if( nodesiz[dir]==msize && nodesiz[dir]%prime[n] == 0) break;

      if (dir < 0) {
        // This cannot happen
        output0 << "CANNOT HAPPEN! in setup_layout_generic.c\n";
        hila::finishrun();
      }

      // Now slice it
      nodesiz[dir] /= prime[n]; 
      nodes.n_divisions[dir] *= prime[n];
    } 

    // now check that the div makes sens
    bool fail = false;
    foralldir(dir) 
      if (nodesiz[dir] < 3) fail = true;  // don't allow nodes of size 1 or 2
    if (fail && !secondtime) {
      secondtime = true;
      ghosts[mdir] = 1<<28;   // this short-circuits direction mdir, some other taken next
    } else if (fail) {
      output0 << "Could not successfully lay out the lattice with " << numnodes() << " nodes\n";
      hila::finishrun();
    }

  } while (secondtime);
  
  // set up struct nodes variables
  nodes.number = numnodes();
  foralldir(dir) {
    nodes.divisors[dir].resize(nodes.n_divisions[dir]+1);
    // Node divisors: note, this MUST BE compatible with
    // node_rank in lattice.cpp
    // to be sure, we use naively the same method than in node_rank
    // last element will be size(dir), for convenience
    int n=-1;
    for (int i=0; i<=size(dir); i++) if ((i*nodes.n_divisions[dir])/size(dir) != n) {
      ++n;
      nodes.divisors[dir][n] = i;
    }
    // output0 << "Divisors ";
    // for (int i=0;i<nodes.n_divisions[dir]; i++) output0 << nodes.divisors[dir][i] << " ";
    // output0 << '\n';
  }


  // Now division done - check how good it is
  int ghost_slices = nsize[mdir] - size(mdir);
  if (ghost_slices > 0) {
    output0 << "\nUsing uneven node division to direction " << mdir << ":\n";
    output0 << "Lengths: " << nodes.n_divisions[mdir]-ghost_slices << " * (" 
            << nodesiz[mdir] << " sites) + " << ghost_slices << " * ("
            << nodesiz[mdir]-1 << " sites)\n";
    output0 << "Divisions: ";
    for (int i=0; i<nodes.n_divisions[mdir]; i++) {
      if (i>0) output0 << " - ";
      output0 << nodes.divisors[mdir][i+1] - nodes.divisors[mdir][i];
    }
    output0 << "\nFilling efficiency: " << (100.0*size(mdir))/nsize[mdir] << "%\n";

    if (ghost_slices > nodes.n_divisions[mdir]/2) 
      output0 << "NOTE: number of smaller nodes > large nodes \n";
  }
  
  // this was numnodes() > 1
  if (1) {   
    output0 << "\nSites on node: ";
    foralldir(dir) {
      if (dir > 0) output0 << " x ";
      if (dir == mdir && ghost_slices > 0)
        output0 << '(' << nodesiz[dir]-1 << '-' << nodesiz[dir] << ')';
      else
        output0 << nodesiz[dir];
    }
    int ns = 1;
    foralldir(dir) ns *= nodesiz[dir];
    if (ghost_slices > 0) {
      int ns2 = ns*(nodesiz[mdir]-1)/nodesiz[mdir];
      output0 << "  =  " << ns2 << " - " << ns << '\n';
    } else {
      output0 << "  =  " << ns << '\n';
    }

    output0 << "Processor layout: ";
    foralldir(dir) {
      if (dir > 0) output0 << " x ";
      output0 << nodes.n_divisions[dir];
    }
    output0 << "  =  " << numnodes() << " nodes\n";
  }
    
  #ifdef USE_MPI
  // For MPI, remap the nodes for periodic torus
  // in the desired manner 
  // we have at least 2 options:
  // map_node_layout_trivial.c
  // map_node_layout_block2.c - for 2^n n.n. blocks
  
  nodes.create_remap();

  #endif  

  print_dashed_line();

}
