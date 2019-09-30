// This routine uses internal MPI (or other) commands to remap the
// node numbers so that the communications should be optimised.  
// This is the trivial version, so that the order is straight
// here.

#include "../plumbing/lattice.h"

void lattice::nodes.create_remap()
{
  output0 << " Node remapping: TRIVIAL (no effort made to reorder)\n";
  nodes.map_array = nodes.map_inverse = nullptr;
}

unsigned lattice::nodes.remap(unsigned i) { return i; }

unsigned lattice::nodes.inverse_remap(unsigned i) { return i; }
