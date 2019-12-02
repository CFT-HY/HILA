#ifndef DEFS_H
#define DEFS_H

// Useful global definitions here -- this file should be included by (almost) all others

#include <array>
#include <vector>
#include <assert.h> 
#include "../plumbing/mersenne.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

#if !defined(CUDA) && !defined(AVX512)
#define VANILLA
#endif

#ifdef CUDA
#define layout_SOA
#endif

#define EVENFIRST

// TODO: default type real_t definition somewhere (makefile?)
using real_t = float;

// move these somewhere - use consts?
// Have this defined in the program?
#ifndef NDIM
  #define NDIM 4
#endif


// HACK  -- this is needed for pragma handlin, do not change!
// #pragma transformer _transformer_cmd_dump_ast_

// HACK
#ifdef TRANSFORMER
#define transformer_ctl(a) extern int _transformer_ctl_##a
#else
#define transformer_ctl(a)
#endif
//void transformer_control(const char *);


// Direction and parity

#if NDIM==4
enum direction :unsigned { XUP, YUP, ZUP, TUP, TDOWN, ZDOWN, YDOWN, XDOWN, NDIRS };
#elif NDIM==3
enum direction { XUP, YUP, ZUP, ZDOWN, YDOWN, XDOWN, NDIRS };
#elif NDIM==2
enum direction { XUP, YUP, YDOWN, XDOWN, NDIRS };
#elif NDIM==1
enum direction { XUP, XDOWN, NDIRS };
#endif

/**
 * Increment op for directions
 * */
#pragma transformer loop_function
inline direction & operator++(direction & dir, int dummy){
  const int i = static_cast<int>(dir);
  return dir=static_cast<direction>((i + 1)%NDIRS);
}

static inline direction opp_dir(const direction d) { return static_cast<direction>(NDIRS - 1 - static_cast<int>(d)); }
static inline direction opp_dir(const int d) { return static_cast<direction>(NDIRS - 1 - d); }

enum class parity : unsigned { none, even, odd, all, x };
// use here #define instead of const parity. Makes EVEN a protected symbol
const parity EVEN = parity::even;
const parity ODD  = parity::odd;
const parity ALL  = parity::all;
const parity X    = parity::x;

// turns EVEN <-> ODD, ALL remains.  X->none, none->none
static inline parity opp_parity(const parity p) {
  unsigned u = 0x3 & static_cast<unsigned>(p);
  return static_cast<parity>(0x3 & ((u<<1)|(u>>1)));
}

/// Return a vector for iterating over  parities included in par
/// If par is ALL, this returns vector of EVEN and ODD, otherwise
/// just par
static std::vector<parity> loop_parities(parity par){
  std::vector<parity> parities;
  if( par == ALL){
    parities.insert(parities.end(), { EVEN, ODD });
  } else {
    parities.insert(parities.end(), { par });
  }
  return parities;
}



#define foralldir(d) for(direction d=XUP; d<NDIM; d++) 

static inline int is_up_dir(const int d) { return d<NDIM; }


// location type

struct location {
    int r[NDIM];
    int& operator[] (const int i) { return r[i]; }
    int& operator[] (const direction d) { return r[(int)d]; }
    const int& operator[] (const int i) const { return r[i]; }
    const int& operator[] (const direction d) const { return r[(int)d]; }
};

inline location operator+(const location & a, const location & b) {
  location r;
  foralldir(d) r[d] = a[d] + b[d];
  return r;
}

inline location operator-(const location & a, const location & b) {
  location r;
  foralldir(d) r[d] = a[d] - b[d];
  return r;
}

inline location operator-(const location & a) {
  location r;
  foralldir(d) r[d] = -a[d];
  return r;
}

inline parity location_parity(const location & a) {
  int s = 0;
  foralldir(d) s += a[d];
  if (s % 2 == 0) return parity::even;
  else return parity::odd;
}

// Replaced by transformer
location coordinates(parity X);




// Global functions: setup
void initial_setup(int argc, char **argv);

// Communication functions
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

#define MAX_GATHERS 1000


inline void assert_even_odd_parity( parity p ) {
    assert(p == EVEN || p == ODD || p == ALL);
}



#ifdef CUDA
#include "../plumbing/hila_cuda.h"

#elif defined(openacc)

#define seed_random(seed) seed_mersenne(seed)
#define hila_random() mersenne()
inline void synchronize_threads(){}

#else

#define seed_random(seed) seed_mersenne(seed)
#define hila_random() mersenne()
inline void synchronize_threads(){}

#endif


#if defined(PUHTI) && defined(TRANSFORMER)
namespace std {
  // This is missing in c++11, which appears to be what we have on Puhti
  template< bool B, class T = void >
  using enable_if_t = typename std::enable_if<B,T>::type;
}
#endif



// Synchronization
#ifndef USE_MPI

inline void synchronize(){
  synchronize_threads();
}

#else

inline void synchronize(){
  static int n=1;
  //printf("node %d, in barrier %d\n", mynode(), n);
  synchronize_threads();
  //printf("node %d, waiting for mpi in barrier %d\n", mynode(), n);
  MPI_Barrier(MPI_COMM_WORLD); 
  //printf("node %d, barrier cleared %d\n", mynode(), n);
  n++;
}

#endif





#endif
