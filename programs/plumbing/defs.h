#ifndef DEFS_H
#define DEFS_H

// Useful global definitions here -- this file should be included by (almost) all others

#include <array>
#include <vector>
#include <assert.h> 
#include "../plumbing/mersenne.h"

#ifdef AVX
#include "../plumbing/AVX.h"
#endif


#ifdef USE_MPI
#include <mpi.h>
#endif

#if !defined(CUDA) && !defined(AVX)
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
enum direction :unsigned { XUP = 0, YUP, ZUP, TUP, TDOWN, ZDOWN, YDOWN, XDOWN, NDIRS };
#elif NDIM==3
enum direction :unsigned { XUP = 0, YUP, ZUP, ZDOWN, YDOWN, XDOWN, NDIRS };
#elif NDIM==2
enum direction :unsigned { XUP = 0, YUP, YDOWN, XDOWN, NDIRS };
#elif NDIM==1
enum direction :unsigned { XUP = 0, XDOWN, NDIRS };
#endif

// Increment/decrement ops for directions
// Post-increment 
#pragma transformer loop_function
inline direction operator++(direction & dir, int dummy){
  const unsigned d = dir;
  dir = static_cast<direction>(d + 1);
  return static_cast<direction>(d);
}

// Pre-increment
#pragma transformer loop_function
inline direction & operator++(direction & dir) {
  dir = static_cast<direction>(dir + 1);
  return dir;
}

static inline direction opp_dir(const direction d) { return static_cast<direction>(NDIRS - 1 - static_cast<int>(d)); }
static inline direction opp_dir(const int d) { return static_cast<direction>(NDIRS - 1 - d); }

enum class parity : unsigned { none = 0, even, odd, all, x };
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

#define foralldir(d) for(direction d=XUP; d<NDIM; ++d)

static inline int is_up_dir(const int d) { return d<NDIM; }


// location type

class coordinate_vector {
 private:
  int r[NDIM];

 public:
  int& operator[] (const int i) { return r[i]; }
  int& operator[] (const direction d) { return r[(int)d]; }
  const int& operator[] (const int i) const { return r[i]; }
  const int& operator[] (const direction d) const { return r[(int)d]; }

  // Parity of this coordinate
  parity coordinate_parity() {
    int s = 0;
    foralldir(d) s += r[d];
    if (s % 2 == 0) return parity::even;
    else return parity::odd;
  }

  // cast to std::array
  operator std::array<int,NDIM>(){std::array<int,NDIM> a; for(int d=0; d<NDIM;d++) a[d] = r[d]; return a;}
};

inline coordinate_vector operator+(const coordinate_vector & a, const coordinate_vector & b) {
  coordinate_vector r;
  foralldir(d) r[d] = a[d] + b[d];
  return r;
}

inline coordinate_vector operator-(const coordinate_vector & a, const coordinate_vector & b) {
  coordinate_vector r;
  foralldir(d) r[d] = a[d] - b[d];
  return r;
}

inline coordinate_vector operator-(const coordinate_vector & a) {
  coordinate_vector r;
  foralldir(d) r[d] = -a[d];
  return r;
}

inline coordinate_vector operator*(const int i, const coordinate_vector & a) {
  coordinate_vector r;
  foralldir(d) r[d] = i*a[d];
  return r;
}

inline coordinate_vector operator*(const coordinate_vector & a, const int i) {
  return i*a;
}

inline coordinate_vector operator/(const coordinate_vector & a, const int i) {
  coordinate_vector r;
  foralldir(d) r[d] = a[d]/i;
  return r;
}

// Replaced by transformer
coordinate_vector coordinates(parity X);


/// Parity + dir -type: used in expressions of type f[X+dir]
/// It's a dummy type, will be removed by transformer
struct parity_plus_direction {
  parity p;
  direction d;
};

/// Declarations, no need to implement these (type removed by transformer)
const parity_plus_direction operator+(const parity par, const direction d);
const parity_plus_direction operator-(const parity par, const direction d);

/// Parity + coordinate offset, used in f[X+coordinate_vector] or f[X+dir1+dir2] etc.
struct parity_plus_offset {
  parity p;
  coordinate_vector cv;
};

const parity_plus_offset operator+(const parity par, const coordinate_vector & cv);
const parity_plus_offset operator+(const parity_plus_direction, const direction d);


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
inline double hila_random(){ return mersenne(); }
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




///Implements test for arithmetic operators in types, similar to 
///std::is_arithmetic but allows vector types

#ifndef VECTORIZED
template< class T >
struct is_arithmetic : std::integral_constant<
  bool,
  std::is_arithmetic<T>::value
> {};
#else
template< class T >
struct is_arithmetic : std::integral_constant<
  bool,
  std::is_arithmetic<T>::value ||
  std::is_same<T,Vec4d>::value ||
  std::is_same<T,Vec8f>::value ||
  std::is_same<T,Vec8i>::value ||
  std::is_same<T,Vec8d>::value ||
  std::is_same<T,Vec16f>::value ||
  std::is_same<T,Vec16i>::value 
> {};
#endif




#endif
