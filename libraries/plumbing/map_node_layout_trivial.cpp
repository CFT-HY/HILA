// This routine uses internal MPI (or other) commands to remap the
// node numbers so that the communications should be optimised.
// This is the trivial version, so that the order is straight
// here.

#include "plumbing/globals.h"

void lattice_struct::allnodes::create_remap() {
    output0 << "Node remapping: TRIVIAL (no effort made to reorder)\n";
    map_array = map_inverse = nullptr;
}

unsigned lattice_struct::allnodes::remap(unsigned i) { return i; }

unsigned lattice_struct::allnodes::inverse_remap(unsigned i) { return i; }
