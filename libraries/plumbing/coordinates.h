#ifndef COORDINATES_H_
#define COORDINATES_H_
/// This header file defines:
///   enum direction
///   enum class parity
///   class coordinate_vector
/// These are used to traverse the lattice coordinate systems

#include "plumbing/defs.h"

///////////////////////////////////////////////////////////////////////////
/// enum direction - includes the opposite direction
/// defined as unsigned, but note that direction + int is not defined
/// Direction can be used as an array index (interchangably with int)
///////////////////////////////////////////////////////////////////////////

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

#pragma hila loop_function
static inline direction next_direction(direction dir) {
  return static_cast<direction>(static_cast<unsigned>(dir)+1);
}
#pragma hila loop_function
static inline direction operator++(direction & dir) {
  return dir = static_cast<direction>(static_cast<unsigned>(dir)+1);
}
#pragma hila loop_function
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

static inline direction abs(direction dir) { if (is_up_dir(dir)) return dir; else return -dir; }

inline int dir_dot_product(direction d1, direction d2) {
  if (d1 == d2) return 1;
  else if (d1 == opp_dir(d2)) return -1;
  else return 0;
}

/// dir_mask_t  type  used in masking directions
/// unsigned char is ok up to 4 dim (2*4 bits)
using dir_mask_t = unsigned char;

inline dir_mask_t get_dir_mask(const direction d) { return (dir_mask_t)(1<<d); }


////////////////////////////////////////////////////////////////////
/// enum class parity type definition - stronger protection than std enum
/// 
////////////////////////////////////////////////////////////////////

enum class parity : unsigned { none = 0, even, odd, all };
// 
// should use here #define instead of const parity? Makes EVEN a protected symbol
constexpr parity EVEN = parity::even;      // bit pattern:  001
constexpr parity ODD  = parity::odd;       //               010
constexpr parity ALL  = parity::all;       //               011

// this is used in diagnostics - make static inline so can be defd here
inline const char * parity_name(parity p) {
  const char * parity_name_s[4] = {"parity::none", "EVEN", "ODD", "ALL"};
  return parity_name_s[(int)p]; 
}

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


//////////////////////////////////////////////////////////////////////
/// COORDINATE_VECTOR type 
//////////////////////////////////////////////////////////////////////

class coordinate_vector {
 private:
  int r[NDIM];

 public:
  coordinate_vector() = default;

  #pragma hila loop_function
  coordinate_vector(const coordinate_vector & v) = default;
  // coordinate_vector(const coordinate_vector & v) {
  //   foralldir(d) r[d] = v[d];
  // }

  ~coordinate_vector() = default;

  #pragma hila loop_function
  coordinate_vector(const int i) { set(i); }

  // initialize with direction -- useful for automatic conversion
  #pragma hila loop_function
  coordinate_vector(const direction & dir) {
    foralldir(d) r[d] = dir_dot_product(d,dir);
  }

  #pragma hila loop_function
  int& operator[] (const int i) { return r[i]; }
  #pragma hila loop_function
  int& operator[] (const direction d) { return r[(int)d]; }
  #pragma hila loop_function
  const int& operator[] (const int i) const { return r[i]; }
  #pragma hila loop_function
  const int& operator[] (const direction d) const { return r[(int)d]; }

  // Parity of this coordinate
  #pragma hila loop_function
  ::parity parity() {
    int s = 0;
    foralldir(d) s += r[d];
    if (s % 2 == 0) return parity::even;
    else return parity::odd;
  }

  // cast to std::array
  #pragma hila loop_function
  operator std::array<int,NDIM>(){std::array<int,NDIM> a; for(int d=0; d<NDIM;d++) a[d] = r[d]; return a;}

  #pragma hila loop_function
  coordinate_vector & set(int i) {
    foralldir(d) r[d] = i;
    return *this;
  }

  #pragma hila loop_function
  coordinate_vector & operator+=(const coordinate_vector &v) {
    foralldir(d) r[d] += v[d];
    return *this;
  }

  #pragma hila loop_function
  coordinate_vector & operator-=(const coordinate_vector &v) {
    foralldir(d) r[d] -= v[d];
    return *this;
  }

  #pragma hila loop_function
  coordinate_vector & operator*=(const int i) {
    foralldir(d) r[d] *= i;
    return *this;
  }

  #pragma hila loop_function
  coordinate_vector & operator/=(const int i) {
    foralldir(d) r[d] /= i;
    return *this;
  }

  // and also additions for direction -- dir acts like a unit vector
  #pragma hila loop_function
  coordinate_vector & operator+=(const direction dir) {
    if (is_up_dir(dir)) ++r[dir]; 
    else --r[-dir];
    return *this;
  }

  #pragma hila loop_function
  coordinate_vector & operator-=(const direction dir) {
    if (is_up_dir(dir)) --r[dir]; 
    else ++r[-dir];
    return *this;
  }

  // and unary -
  #pragma hila loop_function
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

/// Somewhat unorthodox vector * vector and vector / vector -operators; elem by elem

inline coordinate_vector operator*(const coordinate_vector & a, const coordinate_vector & b) {
  coordinate_vector r;
  foralldir(d) r[d] = a[d]*b[d];
  return r;
}

inline coordinate_vector operator/(const coordinate_vector & a, const coordinate_vector & b) {
  coordinate_vector r;
  foralldir(d) r[d] = a[d]/b[d];
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
  int r = 0;
  foralldir(d) r += d1[d]*d2[d];
  return r;
}

inline int norm_sq(const coordinate_vector & d) {
  return dot(d,d);
}


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

// finally, output
inline std::ostream& operator<<(std::ostream &strm, const coordinate_vector & c) {
  foralldir(d) strm << c[d] << ' ';
  return strm;
}



/////////////////////////////////////////////////////////////////////
///  X-coordinate type - "dummy" class, used only
///  in loop index and is removed by code analysis
///  Has nothing inside except some methods
/////////////////////////////////////////////////////////////////////

class X_index_type {
  public:
    const coordinate_vector & coordinates() const;
    int coordinate(direction d) const;
    ::parity parity() const;
};

/// this defines the "point" dummy variable!
static const X_index_type X;


/// X + dir -type: used in expressions of type f[X+dir]
/// It's a dummy type, will be removed by hilapp
struct X_plus_direction {};

/// Declarations, no need to implement these (type removed by hilapp)
const X_plus_direction operator+(const X_index_type x, const direction d);
const X_plus_direction operator-(const X_index_type x, const direction d);

/// X + coordinate offset, used in f[X+coordinate_vector] or f[X+dir1+dir2] etc.
struct X_plus_offset {};

// and prototypes, these operators are not defined anywhere but needed for types
const X_plus_offset operator+(const X_index_type x, const coordinate_vector & cv);
const X_plus_offset operator-(const X_index_type x, const coordinate_vector & cv);
const X_plus_offset operator+(const X_plus_direction, const direction d);
const X_plus_offset operator-(const X_plus_direction, const direction d);
const X_plus_offset operator+(const X_plus_direction, const coordinate_vector & cv);
const X_plus_offset operator-(const X_plus_direction, const coordinate_vector & cv);
const X_plus_offset operator+(const X_plus_offset, const direction d);
const X_plus_offset operator-(const X_plus_offset, const direction d);
const X_plus_offset operator+(const X_plus_offset, const coordinate_vector & cv);
const X_plus_offset operator-(const X_plus_offset, const coordinate_vector & cv);


/// Prototypes for coordinates, replaced by hilapp -- TO BE REMOVED!
const coordinate_vector & coordinates(X_index_type x);
// const coordinate_vector & coordinates(X_plus_direction xd);
// const coordinate_vector & coordinates(X_plus_offset xo);



#endif
