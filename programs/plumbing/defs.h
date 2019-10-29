#ifndef DEFS_H
#define DEFS_H

#ifdef PUHTI_TRANSFORMER_CUDA
#define PUHTI_TRANSFOMER
#define CUDA
#endif

#ifdef PUHTI_TRANSFORMER_openacc
#define PUHTI_TRANSFOMER
#define openacc
#endif

// Useful global definitions here -- this file should be included by (almost) all others

#include <array>
#include <vector>
#include <assert.h> 
#include "../plumbing/mersenne.h"

#define EVENFIRST
#define layout_SOA

// TODO: default type real_t definition somewhere (makefile?)
using real_t = float;

// move these somewhere - use consts?
// Have this defined in the program?
#ifndef NDIM
  #define NDIM 4
#endif

#ifndef USE_MPI
  #define VANILLA
#endif

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

#define foralldir(d) for(direction d=XUP; d<NDIM; d++) 

static inline int is_up_dir(const int d) { return d<NDIM; }


// location type


using location = std::array<int,NDIM>;

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
inline location coordinates(parity X){location l; return l;};



#ifndef USE_MPI
#define mynode() 0
#define numnodes() 1
inline void finishrun() {
  exit(1);
}
#endif

#define MAX_GATHERS 1000


inline void assert_even_odd_parity( parity p ) {
    assert(p == EVEN || p == ODD || p == ALL);
}



#ifdef CUDA
#include "../plumbing/hila_cuda.h"

#elif defined(openacc)

//#include <openacc.h>

#define loop_callable _Pragma("acc routine seq")
#define seed_random(seed) seed_mersenne(seed)
#define hila_random() mersenne()

#else
#define seed_random(seed) seed_mersenne(seed)
#define hila_random() mersenne()
#define loop_callable
#endif


#ifdef PUHTI_TRANSFOMER
namespace std {
  // This is missing in c++11, which appears to be what we have on Puhti
  template< bool B, class T = void >
  using enable_if_t = typename std::enable_if<B,T>::type;
}
#endif


#endif
