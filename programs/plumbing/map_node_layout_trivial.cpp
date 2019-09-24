// This routine uses internal MPI (or other) commands to remap the
// node numbers so that the communications should be optimised.  
// This is the trivial version, so that the order is straight
// here.

#include "../plumbing/lattice.h"

void lattice::map_node_layout( int size[NDIM], 
                               int g[NDIM], 
                               int *map_arr )
{
  output0 << " map_node_layout: TRIVIAL (no effort made to reorder)\n";
}

unsigned lattice::remap_node(unsigned i) { return i; }
