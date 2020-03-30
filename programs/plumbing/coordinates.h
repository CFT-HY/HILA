#ifndef COORDINATES_H_
#define COORDINATES_H_
/// This header file defines:
///   enum direction
///   enum class parity
///   class coordinate_vector
/// These are used to traverse the lattice coordinate systems

#include "../plumbing/defs.h"

/// enum direction - includes the opposite direction
/// defined as unsigned, but note that direction + int is not defined

#if NDIM==4
enum direction : unsigned { XUP = 0, YUP, ZUP, TUP, TDOWN, ZDOWN, YDOWN, XDOWN, NDIRECTIONS };
#elif NDIM==3
enum direction : unsigned { XUP = 0, YUP, ZUP, ZDOWN, YDOWN, XDOWN, NDIRECTIONS };
#elif NDIM==2
enum direction : unsigned { XUP = 0, YUP, YDOWN, XDOWN, NDIRECTIONS };
#elif NDIM==1
enum direction : unsigned { XUP = 0, XDOWN, NDIRECTIONS };
#endif

constexpr unsigned NDIRS = NDIRECTIONS;    // 

// Increment for directions:  ++dir,  dir++  does the obvious 
// dir-- not defined, should we?

#pragma transformer loop_function
static inline direction next_direction(direction dir) {
  return static_cast<direction>(static_cast<unsigned>(dir)+1);
}
#pragma transformer loop_function
static inline direction operator++(direction & dir) {
  return dir = static_cast<direction>(static_cast<unsigned>(dir)+1);
}
#pragma transformer loop_function
static inline direction operator++(direction & dir, int) {
  direction d = dir;
  ++dir;
  return d;
}

/// Basic direction looper, which defines the direction as you go.  
#define foralldir(d) for(direction d=XUP; d<NDIM; ++d)

static inline direction opp_dir(const direction d) { return static_cast<direction>(NDIRS - 1 - static_cast<int>(d)); }
static inline direction opp_dir(const int d) { return static_cast<direction>(NDIRS - 1 - d); }
/// unary + and -
static inline direction operator-(const direction d) { return opp_dir(d); }
static inline direction operator+(const direction d) { return d; }

/// is_up_dir is true if the dir is "up" to coord direction
static inline int is_up_dir(const direction d) { return d<NDIM; }
static inline int is_up_dir(const int d) { return d<NDIM; }

inline int dir_dot_product(direction d1, direction d2) {
  if (d1 == d2) return 1;
  else if (d1 == opp_dir(d2)) return -1;
  else return 0;
}


/// enum class parity type definition

enum class parity : unsigned { none = 0, even, odd, all, x };
// should use here #define instead of const parity? Makes EVEN a protected symbol
constexpr parity EVEN = parity::even;      // bit pattern:  001
constexpr parity ODD  = parity::odd;       //               010
constexpr parity ALL  = parity::all;       //               011
constexpr parity X    = parity::x;         //               100

// utilities for getting the bit patterns
static inline unsigned parity_bits(parity p) {
  return 0x3 & static_cast<unsigned>(p);
}
static inline unsigned parity_bits_inverse(parity p) {
  return 0x3 & ~static_cast<unsigned>(p);
}

// turns EVEN <-> ODD, ALL remains.  X->none, none->none
static inline parity opp_parity(const parity p) {
  unsigned u = parity_bits(p);
  return static_cast<parity>(0x3 & ((u<<1)|(u>>1)));
}

static inline bool is_even_odd_parity( parity p ) { return (p == EVEN || p == ODD); }

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

// COORDINATE_VECTOR

class coordinate_vector {
 private:
  int r[NDIM];

 public:
  coordinate_vector() = default;

  #pragma transformer loop_function
  coordinate_vector(const coordinate_vector & v) = default;
  // coordinate_vector(const coordinate_vector & v) {
  //   foralldir(d) r[d] = v[d];
  // }

  ~coordinate_vector() = default;

  #pragma transformer loop_function
  coordinate_vector(const int i) { set(i); }

  // initialize with direction -- useful for automatic conversion
  #pragma transformer loop_function
  coordinate_vector(const direction & dir) {
    foralldir(d) r[d] = dir_dot_product(d,dir);
  }

  #pragma transformer loop_function
  int& operator[] (const int i) { return r[i]; }
  #pragma transformer loop_function
  int& operator[] (const direction d) { return r[(int)d]; }
  #pragma transformer loop_function
  const int& operator[] (const int i) const { return r[i]; }
  #pragma transformer loop_function
  const int& operator[] (const direction d) const { return r[(int)d]; }

  // Parity of this coordinate
  #pragma transformer loop_function
  parity coordinate_parity() {
    int s = 0;
    foralldir(d) s += r[d];
    if (s % 2 == 0) return parity::even;
    else return parity::odd;
  }

  // cast to std::array
  #pragma transformer loop_function
  operator std::array<int,NDIM>(){std::array<int,NDIM> a; for(int d=0; d<NDIM;d++) a[d] = r[d]; return a;}

  #pragma transformer loop_function
  coordinate_vector & set(int i) {
    foralldir(d) r[d] = i;
    return *this;
  }

  #pragma transformer loop_function
  coordinate_vector & operator+=(const coordinate_vector &v) {
    foralldir(d) r[d] += v[d];
    return *this;
  }

  #pragma transformer loop_function
  coordinate_vector & operator-=(const coordinate_vector &v) {
    foralldir(d) r[d] -= v[d];
    return *this;
  }

  #pragma transformer loop_function
  coordinate_vector & operator*=(const int i) {
    foralldir(d) r[d] *= i;
    return *this;
  }

  #pragma transformer loop_function
  coordinate_vector & operator/=(const int i) {
    foralldir(d) r[d] /= i;
    return *this;
  }

  // and also additions for direction -- dir acts like a unit vector
  #pragma transformer loop_function
  coordinate_vector & operator+=(const direction dir) {
    if (is_up_dir(dir)) ++r[dir]; 
    else --r[dir];
    return *this;
  }

  #pragma transformer loop_function
  coordinate_vector & operator-=(const direction dir) {
    if (is_up_dir(dir)) --r[dir]; 
    else ++r[dir];
    return *this;
  }

  // and unary -
  #pragma transformer loop_function
  coordinate_vector operator-() {
    coordinate_vector v;
    foralldir(d) v[d] = -r[d];
    return v;
  }

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

/// Positive mod(): we define here the positive mod utility function mod(a,b).
/// It mods the arguments 0 <= a < m.  This is not the standard
/// for integer operator% in c++, which gives negative results if a<0.  This is useful in 
/// calculating lattice vector additions on a periodic box

static inline int mod(const int a, const int b) {
  int r = a % b;
  if (r < 0) r += b;
  return r;
}

/// Positive mod for coordinate vector, see  int mod(int a, int b).  If 
/// 2nd argument m is lattice.size(), this mods the vector a to periodic lattice.

inline coordinate_vector mod(const coordinate_vector & a, const coordinate_vector & m) {
  coordinate_vector r;
  foralldir(d) {
    r[d] = mod(a[d], m[d]);
  }
  return r;
}

// dot product, just in case...
inline int dot(const coordinate_vector & d1, const coordinate_vector & d2) {
  int r = 1;
  foralldir(d) r += d1[d]*d2[d];
  return r;
}


// Replaced by transformer
coordinate_vector coordinates(parity X);

/// Special direction operators: dir + dir -> coordinate_vector
inline coordinate_vector operator+(const direction d1, const direction d2) {
  coordinate_vector r;
  foralldir(d) {
    r[d]  = dir_dot_product(d1,d);
    r[d] += dir_dot_product(d2,d);
  }
  return r;
}

inline coordinate_vector operator-(const direction d1, const direction d2) {
  coordinate_vector r;
  foralldir(d) {
    r[d]  = dir_dot_product(d1,d);
    r[d] -= dir_dot_product(d2,d);
  }
  return r;
}

/// Special operators: int*direction -> coordinate_vector
inline coordinate_vector operator*(const int i, const direction dir) {
  coordinate_vector r;
  foralldir(d) r[d] = i*dir_dot_product(d,dir);
  return r;
}

inline coordinate_vector operator*(const direction d, const int i) {
  return i*d;
}

// coordinate vector + direction -- dir is a unit vector
inline coordinate_vector operator+( coordinate_vector cv, const direction dir ) {
  cv += dir;
  return cv;
}
inline coordinate_vector operator-( coordinate_vector cv, const direction dir ) {
  cv -= dir;
  return cv;
}

inline coordinate_vector operator+( const direction dir, coordinate_vector cv ) {
  cv += dir;
  return cv;
}
inline coordinate_vector operator-( const direction dir, coordinate_vector cv ) {
  foralldir(d) cv[d] = dir_dot_product(dir,d) - cv[d];
  return cv;
}


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
const parity_plus_offset operator-(const parity par, const coordinate_vector & cv);
const parity_plus_offset operator+(const parity_plus_direction, const direction d);
const parity_plus_offset operator-(const parity_plus_direction, const direction d);
const parity_plus_offset operator+(const parity_plus_direction, const coordinate_vector & cv);
const parity_plus_offset operator-(const parity_plus_direction, const coordinate_vector & cv);
const parity_plus_offset operator+(const parity_plus_offset, const direction d);
const parity_plus_offset operator-(const parity_plus_offset, const direction d);
const parity_plus_offset operator+(const parity_plus_offset, const coordinate_vector & cv);
const parity_plus_offset operator-(const parity_plus_offset, const coordinate_vector & cv);


#endif
