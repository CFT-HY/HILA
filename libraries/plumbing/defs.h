#ifndef DEFS_H_
#define DEFS_H_

// Useful global definitions here -- this file should be included by (almost) all others

#include <iostream>
#include <array>
#include <vector>
#include <assert.h>
#include <sstream>

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef HILAPP
// The compiler is hilapp
#define __device__
#define __host__
#endif

#include "plumbing/mersenne.h"
#include "plumbing/memalloc.h" // memory allocator
#include "plumbing/timing.h"

/// Allow other than periodic boundary conditions
// #define SPECIAL_BOUNDARY_CONDITIONS


/// Define __restrict__?  It is non-standard but supported by most (all?) compilers.
/// ADD HERE GUARD FOR THOSE WHICH DO not HAVE IT
#define RESTRICT __restrict__
// #ifndef CUDA
// #define RESTRICT __restrict__
// #else
// #define RESTRICT // disabled here
// #endif

#define EVEN_SITES_FIRST

// TODO: default type real_t definition somewhere (makefile?)
using real_t = double;

// move these somewhere - use consts?
// Have this defined in the program?
#ifndef NDIM
#define NDIM 4
#endif

// This below declares "output_only" -qualifier.  It is empty on purpose. Do not remove!
#define output_only

// text output section -- defines also output0, which writes from node 0 only
namespace hila {

/// this is our default output file stream
extern std::ostream output;
/// this is just a hook to store output file, if it is in use
extern std::ofstream output_file;

///  rank of this node
int myrank();

/// Seed random generators
void seed_random(unsigned long seed);

// about_to_finish becomes true at the end.  Signals that
// better not rely on MPI or existence of objects any more.
extern bool about_to_finish;

// check_input is used to notify that we're just checking the
// input values and will exit before fields are allocated.
extern bool check_input;
extern int check_with_nodes;

void initialize(int argc, char **argv);
void finishrun();
void terminate(int status);
void error(const std::string &msg);
void error(const char *msg);

} // namespace hila

// The logger uses hila::myrank, so it cannot be included on top
#include "plumbing/logger.h"
namespace hila {
/// Now declare the logger
extern logger_class log;
} // namespace hila

/// We want to define ostream
///     "output0 << stuff;"
/// which is done only by rank 0.
/// This is hacky but easy.  Probably should be done without #define.
/// Do this through else-branch in order to avoid if-statement problems.
/// #define output0 if (hila::myrank() != 0) {} else hila::output
///
/// Above #define can trigger "dangling-else" warning.  Let us
/// try to avoid it with the following a bit more ugly trick:
#define output0                                                                        \
    for (int _dummy_i_ = 1; hila::myrank() == 0 && _dummy_i_; --_dummy_i_)             \
    hila::output

// The following disables the "dangling-else" warning, but not needed now
//#if defined(__clang__) || defined(__GNUC__)
//#pragma GCC diagnostic ignored "-Wdangling-else"
//#endif

/// define a class for FFT direction
enum class fft_direction { forward, inverse };



// Backend defs-headers

#if defined(CUDA)
#include "plumbing/backend_cuda/defs.h"
#elif defined(AVX)
#include "plumbing/backend_vector/defs.h"
#else
#include "plumbing/backend_cpu/defs.h"
#endif

// this include has to be after the backend defs, because those define hila::random()
#include "plumbing/random.h"

// This contains useful template tools
#include "plumbing/type_tools.h"

#if defined(CUDA)
#include "plumbing/backend_cuda/cuda_templated_ops.h"
#endif


// MPI Related functions and definitions
#define MAX_GATHERS 1000

#define DEFAULT_OUTPUT_NAME "output"

#ifndef USE_MPI

// broadcast does nothing if not MPI
template <typename T> void broadcast(T &v) {}

template <typename T> void broadcast_array(T *var, int n) {}

#endif



int numnodes();
void initialize_communications(int &argc, char ***argv);
void split_into_sublattices(int rank);
void synchronize();
bool is_comm_initialized(void);
void finish_communications();
void abort_communications(int status);

// and print a dashed line
void print_dashed_line(const std::string &txt = {});


#endif
