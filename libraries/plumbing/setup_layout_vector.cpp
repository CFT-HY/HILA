/// Setup layout does the node division.  This version
/// first tries an even distribution, with equally sized
/// nodes, and if that fails allows slightly different
/// node sizes.  

#include "plumbing/globals.h"
#include "plumbing/lattice.h"

/***************************************************************/

/* number of primes to be used in factorization */
#define NPRIMES 12
const static int prime[NPRIMES] = {2,3,5,7,11,13,17,19,23,29,31,37};

// Set up now squaresize and nsquares - arrays 
// Print info to outf as we proceed 

void lattice_struct::setup_layout()
{
  int nfactors[NPRIMES];
  coordinate_vector nodesiz;

  print_dashed_line();
  output0 << "LAYOUT: subnode lattice, with " << VECTOR_SIZE/sizeof(float) << " subnodes\n"
              << "Enabling vector length " << VECTOR_SIZE*8 << " bits = "
              << VECTOR_SIZE/sizeof(double) << " doubles or " 
              << VECTOR_SIZE/sizeof(float) << " floats/ints\n";
  output0 << "Lattice size ";
  foralldir(d) {
    if (d != 0) output0 << " x ";
    output0 << size(d);
  }
  output0 << "  =  " << l_volume << " sites\n";
  output0 << "Dividing to " << numnodes() << " nodes\n";
  output0 << "Layout using vector of " << number_of_subnodes << " elements\n";

  // we want to divide up to numnode * vector_size virtual nodes
  // use the float vector size to divide to numnodes() * vector_size nodes, this is the most
  // demanding.
  // The directions where the extra divisions have been done the node size
  // must be even, so that these can be handled by vectors

  int nn = numnodes();

  //Factorize the node number in primes
  // These factors must be used in slicing the lattice!

  // try divisions:
  // number of virtual nodes
  nn *= number_of_subnodes;

  int i = nn; 
  for (int n=0; n<NPRIMES; n++) {
    nfactors[n] = 0;
    while (i%prime[n] == 0) { i/=prime[n]; nfactors[n]++; }
  }
  if (i != 1) {
    output0 << "Cannot factorize " << numnodes() << " nodes with primes up to " 
            << prime[NPRIMES-1] << '\n';
    finishrun();
  }
  

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

  long long ghost_volume = 1;
  foralldir(d) ghost_volume *= nsize[d];

  // if the division goes even nsize = size() and ghost_volume = volume

  // now try to divide the nsize-volume to subnodes.  We don't try to do
  // the division to direction where there are ghost sites

 
  coordinate_vector divisions,subdiv;

  int mdir;
  bool secondtime = false;
  do {
    // try the division a couple of times, if the 1st fails

    if (ghost_volume > l_volume) {
      mdir = NDIM-1;
      foralldir(j) if (ghosts[mdir] > ghosts[j]) mdir = j;
      // mdir is the direction where we do uneven division (if done)
      // output0 << "MDIR " << mdir << " ghosts mdir " << ghosts[mdir] << " nsize " << nsize[mdir] << '\n';
    } else mdir = -1;

    foralldir(i) {
      nodesiz[i] = (i == mdir) ? nsize[i] : size(i);   // start with ghosted lattice size
      divisions[i] = 1;
    }
    
    for (int n=NPRIMES-1; n>=0; n--) for(i=0; i<nfactors[n]; i++) {
      // figure out which direction to divide -- start from the largest prime, because
      // we don't want this to be last divisor! (would probably wind up with size 1) 
        
      // find largest divisible dimension of h-cubes
      int msize = 1, dir;
      for (dir=0; dir<NDIM; dir++)
        if (nodesiz[dir] > msize && nodesiz[dir]%prime[n] == 0 ) msize = nodesiz[dir];

      // if one direction with largest dimension has already been
      // divided, divide it again.  Otherwise divide first direction
      // with the largest dimension. 

      // Try to keep even node sizes, these are needed for vectors
      // We don't worry about evenness to mdir
      // for the 1st iteration this loop wants even division, 2nd time lets uneven through

      dir = NDIM;
      for (int rem=0; rem < 2 && dir == NDIM; rem++) {

        for (dir=0; dir<NDIM; dir++)
          if (nodesiz[dir]==msize && divisions[dir]>1 && nodesiz[dir]%prime[n] == 0 && 
              (dir == mdir || (nodesiz[dir]/prime[n])%2 <= rem) ) break;

        // If not previously sliced, take one direction to slice - try to keep even
        if (dir == NDIM) for (dir=0; dir<NDIM; dir++)
          if( nodesiz[dir]==msize && nodesiz[dir]%prime[n] == 0 &&
              (dir == mdir || (nodesiz[dir]/prime[n])%2 <= rem) ) break;

      }

      if (dir == NDIM) {
        // This cannot happen
        output0 << "CANNOT HAPPEN! in setup_layout_vector.c\n";
        finishrun();
      }

      // Now slice it
      nodesiz[dir] /= prime[n]; 
      divisions[dir] *= prime[n];

      //  output0 << nodesiz << ' ' << divisions << '\n';
    } 

    // Division done, now check that the div makes sense

    bool fail = false;

    foralldir(dir) 
      if (nodesiz[dir] < 3) fail = true;  // don't allow nodes of size 1 or 2

    if (!fail) {

      // check here that this can be used for vectorized division

      subdiv = 1;
      bool div_done;
      int n_subn = 1;
      do {
        div_done = false;
        foralldir(dir) {
          // the direction where the vector subnodes are must not be 
          // the uneven direction, node size to the direction should be divisible by 2
          // and the number of nodes to this dir should also be a multiple of subdivs

          int sd = subdiv[dir] * 2;
          if (dir != mdir && nodesiz[dir] % 2 == 0 && divisions[dir] % sd == 0 
              && n_subn < number_of_subnodes) {
            subdiv[dir] = sd;
            n_subn *= 2;
            div_done = true;
            this_node.subnodes.merged_subnodes_dir = dir;
          }
        }
      } while (div_done && n_subn < number_of_subnodes);

      if (n_subn != number_of_subnodes) fail = true;

    }

    if (fail && !secondtime && mdir >= 0) {
      secondtime = true;
      ghosts[mdir] = 1<<28;   // this short-circuits direction mdir, some other taken next
    } else if (fail) {
      output0 << "Could not successfully lay out the lattice with " << numnodes() << " nodes\n";
      finishrun();
    }

  } while (secondtime);


  // set up struct nodes variables
  nodes.number = numnodes();
  foralldir(dir) {
    nodesiz[dir] *= subdiv[dir];
    nodes.n_divisions[dir] = divisions[dir]/subdiv[dir];
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
  }

  // set up the subnode divisions here -- rest is set in setup_node
  foralldir(d) this_node.subnodes.divisions[d] = subdiv[d];

  // this_node is set up in setup_node

  // Now division done - check how good it is
  int ghost_slices = 0;
  if (mdir >= 0) {
    ghost_slices = nsize[mdir] - size(mdir);

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
      if (dir == mdir)
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

    output0 << "Node layout: ";
    foralldir(dir) {
      if (dir > 0) output0 << " x ";
      output0 << nodes.n_divisions[dir];
    }
    output0 << "  =  " << numnodes() << " nodes\n";

#ifdef VECTORIZED

    output0 << "Node subdivision to 32bit elems: ";
    foralldir(dir) {
      if (dir > 0) output0 << " x ";
      output0 << subdiv[dir];
    }
    output0 << "  =  " << number_of_subnodes << " subnodes\n";
    
    output0 << "Sites on subnodes: ";
    foralldir(dir) {
      if (dir > 0) output0 << " x ";
      if (dir == mdir)
        output0 << '(' << nodesiz[dir]-1 << '-' << nodesiz[dir] << ')';
      else
        output0 << nodesiz[dir]/subdiv[dir];
    }
    output0 << '\n';

    direction dmerge = this_node.subnodes.merged_subnodes_dir;

    output0 << "Node subdivision to 64bit elems: ";
    foralldir(dir) {
      if (dir > 0) output0 << " x ";
      output0 << ( (dir == dmerge) ? subdiv[dir]/2 : subdiv[dir] );
    }
    output0 << "  =  " << number_of_subnodes/2 << " subnodes\n";
    
    output0 << "Sites on subnodes: ";
    foralldir(dir) {
      if (dir > 0) output0 << " x ";
      if (dir == mdir)
        output0 << '(' << nodesiz[dir]-1 << '-' << nodesiz[dir] << ')';
      else
        output0 << ( (dir == dmerge) ? 2*nodesiz[dir]/subdiv[dir] : nodesiz[dir]/subdiv[dir] );
    }
    output0 << '\n';
    
#endif


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
