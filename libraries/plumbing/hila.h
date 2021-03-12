#ifndef HILA_H_
#define HILA_H_

///////////////////////////////////////////////////////////////////////
/// Catch-(almost)all include to get in most of the hila-system .h -files

#include "plumbing/defs.h"
#include "datatypes/cmplx.h"
#include "datatypes/matrix.h"
#include "plumbing/coordinates.h"
#include "plumbing/lattice.h"
#include "plumbing/field.h"
#include "plumbing/reduction.h"
#include "plumbing/vectorreduction.h"
#include "plumbing/input.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

#endif