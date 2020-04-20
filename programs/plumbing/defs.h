#ifndef DEFS_H
#define DEFS_H

// Useful global definitions here -- this file should be included by (almost) all others

#include <iostream>
#include <array>
#include <vector>
#include <assert.h>
#include <sstream>

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "../plumbing/mersenne.h"
#include "../plumbing/memalloc.h"   // memory allocator


/// Define __restrict__?  It is non-standard but supported by most (all?) compilers.
/// ADD HERE GUARD FOR THOSE WHICH DO not HAVE IT
#if 1
#define RESTRICT __restrict__
#else
#define RESTRICT   // disabled here
#endif

#define EVEN_SITES_FIRST

// TODO: default type real_t definition somewhere (makefile?)
using real_t = float;

// move these somewhere - use consts?
// Have this defined in the program?
#ifndef NDIM
#define NDIM 4
#endif

#define pure_output     // nothing here - keyword for loop function parameters

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




// Global functions: setup
void initial_setup(int argc, char **argv);





// Backend defs-headers

#if defined(CUDA)
#include "../plumbing/backend_cuda/defs.h"
#elif defined(AVX)
#include "../plumbing/backend_vector/defs.h"
#else
#include "../plumbing/backend_cpu/defs.h"
#endif





// MPI Related functions and definitions
#define MAX_GATHERS 1000

#ifndef USE_MPI

// Trivial, no MPI
#define mynode() 0
#define numnodes() 1
inline void initialize_machine(int &argc, char ***argv) {}
inline void finishrun() {
  exit(0);
}

#else

int mynode();
int numnodes();
void finishrun();
void initialize_machine(int &argc, char ***argv);

#endif


// Synchronization
#ifndef USE_MPI

inline void synchronize(){
  synchronize_threads();
}

#else

inline void synchronize(){
  static int n=1;
  synchronize_threads();
  MPI_Barrier(MPI_COMM_WORLD); 
  n++;
}

#endif




// Useful c++14 template missing in Puhti compilation of transformer
#if defined(PUHTI) && defined(TRANSFORMER)
namespace std {
  template< bool B, class T = void >
  using enable_if_t = typename std::enable_if<B,T>::type;
}
#endif


/// Utility for selecting the numeric base type of a class
template<class T, class Enable = void>
struct base_type_struct {
  using type = typename T::base_type;
};

template<typename T>
struct base_type_struct< T, typename std::enable_if_t<is_arithmetic<T>::value>> {
  using type = T;
};




#endif
