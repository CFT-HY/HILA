#ifndef COORDINATES_H_
#define COORDINATES_H_
/// This header file defines:
///   enum direction
///   enum class parity
///   class coordinate_vector
/// These are used to traverse the lattice coordinate systems

#include "plumbing/defs.h"
#include "datatypes/matrix.h"

///////////////////////////////////////////////////////////////////////////
/// enum direction - includes the opposite direction
/// defined as unsigned, but note that direction + int is not defined
/// Direction can be used as an array index (interchangably with int)
///////////////////////////////////////////////////////////////////////////

#if NDIM==4
enum direction : unsigned { e_x = 0, e_y, e_z, e_t, e_t_down, e_z_down, e_y_down, e_x_down, NDIRECTIONS };
#elif NDIM==3
enum direction : unsigned { e_x = 0, e_y, e_z, e_z_down, e_y_down, e_x_down, NDIRECTIONS };
#elif NDIM==2
enum direction : unsigned { e_x = 0, e_y, e_y_down, e_x_down, NDIRECTIONS };
#elif NDIM==1
enum direction : unsigned { e_x = 0, e_x_down, NDIRECTIONS };
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
#define foralldir(d) for(direction d=e_x; d<NDIM; ++d)

#pragma hila loop_function
static inline direction opp_dir(const direction d) { return static_cast<direction>(NDIRS - 1 - static_cast<int>(d)); }
#pragma hila loop_function
static inline direction opp_dir(const int d) { return static_cast<direction>(NDIRS - 1 - d); }
/// unary + and -
#pragma hila loop_function
static inline direction operator-(const direction d) { return opp_dir(d); }
#pragma hila loop_function
static inline direction operator+(const direction d) { return d; }

/// is_up_dir is true if the dir is "up" to coord direction
#pragma hila loop_function
static inline bool is_up_dir(const direction d) { return d<NDIM; }
#pragma hila loop_function
static inline bool is_up_dir(const int d) { return d<NDIM; }

#pragma hila loop_function
static inline direction abs(direction dir) { if (is_up_dir(dir)) return dir; else return -dir; }

#pragma hila loop_function
inline int dir_dot_product(direction d1, direction d2) {
  if (d1 == d2) return 1;
  else if (d1 == -d2) return -1;
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

// This should not be used from loops ...
inline std::ostream & operator<<(std::ostream& os, const parity p) {
  os << parity_name(p);
  return os;
}

// utilities for getting the bit patterns
static inline unsigned parity_bits(parity p) {
  return 0x3 & static_cast<unsigned>(p);
}
static inline unsigned parity_bits_inverse(parity p) {
  return 0x3 & ~static_cast<unsigned>(p);
}

// turns EVEN <-> ODD, ALL remains.  X->none, none->none
#pragma hila loop_function
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

template <typename T>
class coordinate_vector_t : public Vector<NDIM,T> {

 public:

    // std incantation for field types
    using base_type = typename base_type_struct<T>::type;

    // define these to ensure std::is_trivial
    coordinate_vector_t() = default;
    ~coordinate_vector_t() = default;
    coordinate_vector_t(const coordinate_vector_t & v) = default;

    // initialize with direction -- useful for automatic conversion
    #pragma hila loop_function
    inline coordinate_vector_t(const direction dir) {
      foralldir(d) this->e(d) = dir_dot_product(d,dir);
    }

    // construct from vector
    #pragma hila loop_function
    coordinate_vector_t(const Vector<NDIM,T> & v) {
      foralldir(d) this->e(d) = v.e(d);
    }

    // assign and construct from zero
    #pragma hila loop_function
    coordinate_vector_t(const Zero z) {
      foralldir(d) this->e(d) = static_cast<T>(0);
    }

    #pragma hila loop_function
    inline coordinate_vector_t & operator= (const Zero z) {
      static_cast<coordinate_vector_t &>( *this = z );
    }


    #pragma hila loop_function
    T& operator[] (const int i) { return this->e(i); }
    #pragma hila loop_function
    T& operator[] (const direction d) { return this->e((int)d); }
    #pragma hila loop_function
    T operator[] (const int i) const { return this->e(i); }
    #pragma hila loop_function
    T operator[] (const direction d) const { return this->e((int)d); }



    // Parity of this coordinate
    #pragma hila loop_function
    ::parity parity() {
      int s = 0;
      foralldir(d) s += this->e(d);
      if (s % 2 == 0) return parity::even;
      else return parity::odd;
    }

    // cast to std::array
    #pragma hila loop_function
    operator std::array<int,NDIM>(){
      std::array<int,NDIM> a;
      for(int d=0; d<NDIM;d++) a[d] = this->e(d);
      return a;
    }

    // cast to Vector
    #pragma hila loop_function
    operator Vector<NDIM,T>(){
      Vector<NDIM,T> a;
      for(int d=0; d<NDIM;d++) a[d] = this->e(d);
      return a;
    }



    // add coordinate vector
    #pragma hila loop_function
    coordinate_vector_t & operator+=(const coordinate_vector_t & rhs) {
      foralldir(d) this->e(d) += rhs.e(d);
      return *this;
    }


    // and also additions for direction -- dir acts like a unit vector
    #pragma hila loop_function
    coordinate_vector_t & operator+=(const direction dir) {
      if (is_up_dir(dir)) ++this->e(dir);
      else --this->e(-dir);
      return *this;
    }

    #pragma hila loop_function
    coordinate_vector_t & operator-=(const direction dir) {
      if (is_up_dir(dir)) --this->e(dir);
      else ++this->e(-dir);
      return *this;
    }

    // unary -
    #pragma hila loop_function
    inline coordinate_vector_t operator-() const {
      coordinate_vector_t res;
      foralldir(d) res.e(d) = -this->e(d);
      return res;
    }




};

/// Define the std coordinate_vector type here
using coordinate_vector = coordinate_vector_t<int>;

/// Positive mod(): we define here the positive mod utility function mod(a,b).
/// It mods the arguments 0 <= a < m.  This is not the standard
/// for integer operator% in c++, which gives negative results if a<0.  This is useful in 
/// calculating lattice vector additions on a periodic box

#pragma hila loop_function
static inline int mod(const int a, const int b) {
  int r = a % b;
  if (r < 0) r += b;
  return r;
}

/// Positive mod for coordinate vector, see  int mod(int a, int b).  If 
/// 2nd argument m is lattice.size(), this mods the vector a to periodic lattice.

template <typename T>
#pragma hila loop_function
inline coordinate_vector_t<T> mod(const coordinate_vector_t<T> & a, const coordinate_vector_t<T> & m) {
  coordinate_vector_t<T> r;
  foralldir(d) {
    r.e(d) = mod(a[d], m[d]);
  }
  return r;
}

template <typename T>
#pragma hila loop_function
inline coordinate_vector_t<T> operator+(coordinate_vector_t<T> cv1, const coordinate_vector_t<T> & cv2) {
  coordinate_vector_t<T> res;
  foralldir(d) res.c[d] = cv1.c[d] + cv2.c[d];
  return res;
}

template <typename T>
#pragma hila loop_function
inline coordinate_vector_t<T> operator-(coordinate_vector_t<T> cv1, const coordinate_vector_t<T> & cv2) {
  coordinate_vector_t<T> res;
  foralldir(d) res.c[d] = cv1.c[d] - cv2.c[d];
  return res;
}

/// Special direction operators: dir + dir -> coordinate_vector
template <typename T>
#pragma hila loop_function
inline coordinate_vector_t<T> operator+(const direction d1, const direction d2) {
  coordinate_vector_t<T> r;
  foralldir(d) {
    r.e(d)  = dir_dot_product(d1,d);
    r.e(d) += dir_dot_product(d2,d);
  }
  return r;
}

template <typename T>
#pragma hila loop_function
inline coordinate_vector_t<T> operator-(const direction d1, const direction d2) {
  coordinate_vector_t<T> r;
  foralldir(d) {
    r.e(d)  = dir_dot_product(d1,d);
    r.e(d) -= dir_dot_product(d2,d);
  }
  return r;
}

/// Special operators: int*direction -> coordinate_vector (of type int!)
#pragma hila loop_function
inline coordinate_vector operator*(const int i, const direction dir) {
  coordinate_vector r;
  foralldir(d) r.e(d) = i*dir_dot_product(d,dir);
  return r;
}

#pragma hila loop_function
inline coordinate_vector operator*(const direction d, const int i) {
  return i*d;
}

// coordinate vector + direction -- dir is a unit vector
template <typename T>
#pragma hila loop_function
inline coordinate_vector_t<T> operator+( coordinate_vector_t<T> cv, const direction dir ) {
  cv += dir;
  return cv;
}
template <typename T>
#pragma hila loop_function
inline coordinate_vector_t<T> operator-( coordinate_vector_t<T> cv, const direction dir ) {
  cv -= dir;
  return cv;
}

template <typename T>
#pragma hila loop_function
inline coordinate_vector_t<T> operator+( const direction dir, coordinate_vector_t<T> cv ) {
  cv += dir;
  return cv;
}
template <typename T>
#pragma hila loop_function
inline coordinate_vector_t<T> operator-( const direction dir, coordinate_vector_t<T> cv ) {
  foralldir(d) cv.e(d) = dir_dot_product(dir,d) - cv.e(d);
  return cv;
}




/////////////////////////////////////////////////////////////////////
///  X-coordinate type - "dummy" class, used only
///  in loop index and is removed by code analysis
///  Has nothing inside except some methods
/////////////////////////////////////////////////////////////////////

class X_index_type {
  public:
    #pragma hila loop_function
    const coordinate_vector & coordinates() const;

    #pragma hila loop_function
    int coordinate(direction d) const;

    #pragma hila loop_function
    ::parity parity() const;
};

/// this defines the "point" dummy variable!
static const X_index_type X;


/// X + dir -type: used in expressions of type f[X+dir]
/// It's a dummy type, will be removed by hilapp
struct X_plus_direction {};

/// Declarations, no need to implement these (type removed by hilapp)
#pragma hila loop_function
const X_plus_direction operator+(const X_index_type x, const direction d);
#pragma hila loop_function
const X_plus_direction operator-(const X_index_type x, const direction d);

/// X + coordinate offset, used in f[X+coordinate_vector] or f[X+dir1+dir2] etc.
struct X_plus_offset {};

// and prototypes, these operators are not defined anywhere but needed for types
#pragma hila loop_function
const X_plus_offset operator+(const X_index_type x, const coordinate_vector & cv);
#pragma hila loop_function
const X_plus_offset operator-(const X_index_type x, const coordinate_vector & cv);
#pragma hila loop_function
const X_plus_offset operator+(const X_plus_direction, const direction d);
#pragma hila loop_function
const X_plus_offset operator-(const X_plus_direction, const direction d);
#pragma hila loop_function
const X_plus_offset operator+(const X_plus_direction, const coordinate_vector & cv);
#pragma hila loop_function
const X_plus_offset operator-(const X_plus_direction, const coordinate_vector & cv);
#pragma hila loop_function
const X_plus_offset operator+(const X_plus_offset, const direction d);
#pragma hila loop_function
const X_plus_offset operator-(const X_plus_offset, const direction d);
#pragma hila loop_function
const X_plus_offset operator+(const X_plus_offset, const coordinate_vector & cv);
#pragma hila loop_function
const X_plus_offset operator-(const X_plus_offset, const coordinate_vector & cv);


/// Prototypes for coordinates, replaced by hilapp -- TO BE REMOVED!
// const coordinate_vector & coordinates(X_index_type x);
// const coordinate_vector & coordinates(X_plus_direction xd);
// const coordinate_vector & coordinates(X_plus_offset xo);



#endif
