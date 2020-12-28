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
#include "plumbing/memalloc.h"   // memory allocator
#include "plumbing/timing.h"

/// Define __restrict__?  It is non-standard but supported by most (all?) compilers.
/// ADD HERE GUARD FOR THOSE WHICH DO not HAVE IT
#ifndef CUDA
#define RESTRICT __restrict__
#else
#define RESTRICT   // disabled here
#endif

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
  // this is our default output file stream
  extern std::ostream output;
  // this is just a hook to store output file, if it is in use
  extern std::ofstream output_file;

  // store the rank of this process to a global variable - this will not vanish during 
  // object destruction at the end!
  extern int my_rank_n;
  inline int myrank() { return my_rank_n; }

  extern bool about_to_finish;

  void initialize(int argc, char **argv);
  void finishrun();
  void terminate(int status);
  void error(const std::string & msg);
  void error(const char * msg);

}

// We want to define ostream
//     "output0 << stuff;"
// which is done only by rank 0.  
// This is hacky but easy.  Probably should be done without #define.
// Do this through else-branch in order to avoid if-statement problems.
// #define output0 if (hila::myrank() != 0) {} else hila::output
//
// Above #define can trigger "dangling-else" warning.  Let us
// try to avoid it with the following a bit more ugly trick:
#define output0 for(int _dummy_i_=1; hila::myrank()==0 && _dummy_i_; --_dummy_i_) hila::output

// The following disables the "dangling-else" warning, but not needed now
//#if defined(__clang__) || defined(__GNUC__)
//#pragma GCC diagnostic ignored "-Wdangling-else"
//#endif

// define a class for FFT direction
enum class fft_direction { forward, backward, back };

// Allow other than periodic boundary conditions
#define SPECIAL_BOUNDARY_CONDITIONS

// Backend defs-headers

#if defined(CUDA)
#include "plumbing/backend_cuda/defs.h"
#elif defined(AVX)
#include "plumbing/backend_vector/defs.h"
#else
#include "plumbing/backend_cpu/defs.h"
#endif

// this include has to be after the backend defs, because those define hila_random()
#include "plumbing/random.h"

// MPI Related functions and definitions
#define MAX_GATHERS 1000

#define DEFAULT_OUTPUT_NAME "output"

#ifndef USE_MPI

// broadcast does nothing if not MPI
template <typename T>
void broadcast(T & v) {}

template <typename T>
void broadcast_array(T * var, int n) {}

#endif

int mynode();
int numnodes();
void initialize_communications(int &argc, char ***argv);
void split_into_sublattices( int rank );
void synchronize();
bool is_comm_initialized(void);
void finish_communications();
void abort_communications(int status);


// and print a dashed line
void print_dashed_line();

// Useful c++14 template missing in Puhti compilation of hilapp
#if defined(PUHTI) && defined(HILAPP)
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

/// Utility for selecting the numeric base type of a class
template<typename T>
struct base_type_struct< T, typename std::enable_if_t<is_arithmetic<T>::value>> {
  using type = T;
};

/// Utility for selecting the numeric base type of a class
template<typename T>
using number_type = typename base_type_struct<T>::type;
 

// These are helpers, to make generic templates
// e.g. type_plus<A,B> gives the type of the operator a + b, where a is of type A and b B.
template<typename A, typename B>
using type_plus = decltype(std::declval<A>() + std::declval<B>());
template<typename A, typename B>
using type_minus= decltype(std::declval<A>() - std::declval<B>());
template<typename A, typename B>
using type_mul  = decltype(std::declval<A>() * std::declval<B>());
template<typename A, typename B>
using type_div  = decltype(std::declval<A>() / std::declval<B>());



#endif
