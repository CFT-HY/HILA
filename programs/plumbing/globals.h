#ifndef GLOBALS_H
#define GLOBALS_H

/// This file contains general global state of the simulation programs.
/// Probably this could be done much more elegantly

#include <iostream>
#include "../plumbing/lattice.h"

// text output section -- defines also output0, which writes from node 0 only

namespace hila {
  // this is our default output file stream
  extern std::ostream &output;
  // this is just a hook to store output file, if it is in use
  extern std::ofstream output_file;
};

// this is pretty hacky but easy.  Probably could do without #define too
// do this through else-branch in order to avoid if-statement problems
#define output0 if (mynode() != 0) {} else hila::output

#endif
