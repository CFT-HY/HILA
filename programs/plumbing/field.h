// -*- mode: c++ -*-
#ifndef FIELD_H
#define FIELD_H
#include <iostream>
#include <string>
#include <math.h>

// using namespace std;

#include "../infrastructure/lattice.h"


// HACK  -- this is needed for pragma handlin, do not change!
// #pragma transformer _transformer_cmd_dump_ast_

// HACK
#define transformer_ctl(a) extern int _transformer_ctl_##a
//void transformer_control(const char *);

// TODO: default type real_t definition somewhere
typedef double real_t;

enum class direction : unsigned { xup, yup, zup, tup, tdown, zdown, ydown, xdown, NONE };

static inline direction opp_dir(const direction d) { 
  return static_cast<direction>(NDIRS - static_cast<unsigned>(d)); }

enum class parity : unsigned { none, even, odd, all, x };
const parity EVEN = parity::even;
const parity ODD  = parity::odd;
const parity ALL  = parity::all;
const parity X    = parity::x;

// turns EVEN <-> ODD, ALL remains.  X->none, none->none
static inline parity opp_parity(const parity p) {
  unsigned u = 0x3 & static_cast<unsigned>(p);
  return static_cast<parity>(0x3 & ((u<<1)|(u>>1)));
}


struct parity_plus_direction {
  parity p;
  direction d;
};

const parity_plus_direction operator+(const parity par, const direction d);
const parity_plus_direction operator-(const parity par, const direction d);

// This is a marker for transformer -- for does not survive as it is
#define onsites(p) for(parity parity_type_var_(p);;)

//class parity {
//  parity() {}
//  parity( parity &par ) { p = par.p; }
//  parity( parity_enum pare ) { p = pare; }
//  
//  parity & operator=(const parity &par) { p = par.p; return *this; }
//  parity & operator=(const parity_enum par) { p = par; return *this; }
//  const parity_plus_direction operator+(const direction &d) { parity_plus_direction pd; return(pd); }
//  const parity_plus_direction operator-(const direction &d) { parity_plus_direction pd; return(pd); }
//};

// #define onallsites(i) for (int i=0; i<N; i++) 


template <typename R> struct field_vector {
  R vec[32/sizeof(R)];   // TODO: size should be the size of the vector
};

// fwd definition
// template <typename T> class field;

// field_element class: virtual class, no storage allocated,
// wiped out by the transformer
template <typename T>
class field_element  {
 private:
  T v;   // TODO: should rather be a short vector type?
  
 public:
  // the following are just placemarkers, and are substituted by the transformer

  // implicit conversion to type T: this works for storage types where field is an std
  // array of type T - this is again for type propagation
  // TODO: IS THIS NEEDED?  WE WANT TO AVOID CONVERSIONS FROM field<T> v -> T
  // operator T() { return v; }
      
  // The type is important for ensuring correctness
  // Possibility: write these so that they work without the transformer
  field_element<T>& operator= (const T &d) {
    v = d; return *this;}
  
  // field_element = field_element
  field_element<T>& operator=  (const field_element<T>& rhs) {
    v  = rhs.v; return *this;}
  field_element<T>& operator+= (const field_element<T>& rhs) {
    v += rhs.v; return *this;}
  field_element<T>& operator-= (const field_element<T>& rhs) {
    v -= rhs.v; return *this;}
  field_element<T>& operator*= (const field_element<T>& rhs) {
    v *= rhs.v; return *this;}
  field_element<T>& operator/= (const field_element<T>& rhs) {
    v /= rhs.v; return *this;}

  //field_element<T>& operator+= (const T rhs) {
  //  v += rhs; return *this;}
  //field_element<T>& operator-= (const T rhs) {
  //  v -= rhs; return *this;}
  //field_element<T>& operator*= (const T rhs) {
  //  v *= rhs; return *this;}
  //field_element<T>& operator/= (const T rhs) {
  //  v /= rhs; return *this;}

  field_element<T>& operator+= (const double rhs) {
    v += rhs; return *this;}
  field_element<T>& operator-= (const double rhs) {
    v -= rhs; return *this;}
  field_element<T>& operator*= (const double rhs) {
    v *= rhs; return *this;}
  field_element<T>& operator/= (const double rhs) {
    v /= rhs; return *this;}


  // access the raw value - TODO:short vectors 
  T get_value() { return v; }

  T reduce_plus() {
    return v;   // TODO: short vector!
  }

  T reduce_mult() {
    return v;   // TODO: short vector!
  }
  
};

// declarations, implemented by transformer -- not necessarily defined anywhere
// +
template <typename T>
field_element<T> operator+( const field_element<T> &lhs, const field_element<T> &rhs);

//template <typename T>
//field_element<T> operator+( const T &lhs, const field_element<T> &rhs);

//template <typename T>
//field_element<T> operator+( const field_element<T> &lhs,  const T &rhs);

template <typename T,typename L>
field_element<T> operator+( const L &lhs, const field_element<T> &rhs);

template <typename T,typename R>
field_element<T> operator+( const field_element<T> &lhs, const R &rhs);

// -
template <typename T>
field_element<T> operator-( const field_element<T> &lhs, const field_element<T> &rhs);

template <typename T,typename L>
field_element<T> operator-( const L &lhs, const field_element<T> &rhs);

template <typename T,typename R>
field_element<T> operator-( const field_element<T> &lhs,  const R &rhs);

// template <typename T>
// field_element<T> operator-( double lhs, const field_element<T> &rhs);

// template <typename T>
// field_element<T> operator-( const field_element<T> &lhs, double rhs);

// *
template <typename T>
field_element<T> operator*( const field_element<T> &lhs, const field_element<T> &rhs);

template <typename T,typename L>
field_element<T> operator*( const L &lhs, const field_element<T> &rhs);

template <typename T,typename R>
field_element<T> operator*( const field_element<T> &lhs,  const R &rhs);

// template <typename T>
// field_element<T> operator*( double lhs, const field_element<T> &rhs);

// template <typename T>
// field_element<T> operator*( const field_element<T> &lhs, double rhs);


// /    Division is not implemented for all types, but define generically here
template <typename T>
field_element<T> operator/( const field_element<T> &lhs, const field_element<T> &rhs);

template <typename T,typename L>
field_element<T> operator/( const L &lhs, const field_element<T> &rhs);

template <typename T,typename R>
field_element<T> operator/( const field_element<T> &lhs,  const R &rhs);

// template <typename T>
// field_element<T> operator/( const field_element<T> &lhs, double rhs);



// a function
template <typename T>
field_element<T> exp( field_element<T> &arg) {
  field_element<T> res;
  res = exp(arg.get_value());
  return res;
}


// TRY NOW AUTOMATIC REDUCTION IDENTIFICATION
// Overload operator  res += expr, where
// res is type T and expr is field_element<T>
// Make these void, because these cannot be assigned from
// These will be modified by transformer

template <typename T>
void operator += (T& lhs, field_element<T>& rhs) {
  lhs += rhs.reduce_plus();
}

template <typename T>
void operator *= (T& lhs, field_element<T>& rhs) {
  lhs *= rhs.reduce_mult();
}



// Type alias would be nice, but one cannot specialize those!
// Placemarker, will be specialized by transformer
template <typename T>
struct field_storage_type {
  T c;
};

/// The following struct holds the data + information about the field
/// TODO: field-specific boundary conditions?
template <typename T>
class field_struct {
  field_storage_type<T> * payload; // TODO: must be maximally aligned, modifiers - never null
  lattice_struct * lattice;
  unsigned is_fetched[NDIRS];
  
};


// ** class field

template <typename T>
class field {
private:
  // here correct data
  field_struct<T> * fs;
  
public:
  
  field<T>() {
    // std::cout << "In constructor 1\n";
    fs = nullptr;             // lazy allocation on 1st use
  }
  
  // copy constructor  -- let the latfield_element<> type find if this is OK
  template <typename R>
  field<T>(const field<R>& other) {
    fs = nullptr;
    (*this)[ALL] = other[X];
  }

  template <typename R>
  field<T>(const R& val) {
    fs = nullptr;
    (*this)[ALL] = val;
  }
  
  // move constructor - steal the content
  field<T>(field<T>&& rhs) {
    // std::cout << "in move constructor\n";
    fs = rhs.fs;
    rhs.fs = nullptr;
  }

  ~field<T>() {
    free();
  }
    
  void alloc() {
    assert(fs == nullptr);
    if (current_lattice == nullptr) {
      // TODO: write to some named stream
      std::cout << "Can not use field variables before lattice is initialized\n";
      exit(1);  // TODO - more ordered exit?
    }
    fs = new field_struct<T>;
    fs->lattice = &current_lattice;
    fs->payload = (field_storage_type<T> *)alloc_field_memory(sizeof(T)*current_lattice->node_fullsize());

    mark_changed(ALL);
  }

  void free() {
    if (fs != nullptr) {
      free_field_memory( fs->payload );
      delete fs;
      fs = nullptr;
    }
  }

  // call this BEFORE the var is written to
  void mark_changed(const parity p) {
    if (fs == nullptr) alloc();
    else {
      char pc = static_cast<char>(p);
      assert(p == EVEN || p == ODD || p == ALL);
      unsigned up = 0x3 & (!(static_cast<unsigned>(opp_parity(p))));
      for (int i=0; i<NDIRS; i++) fs->is_fetched[i] &= up;
    }
  }
  
  void assert_is_initialized() {
    if (fs == nullptr) {
      std::cout << "field variable used before it is assigned to\n";
      exit(1);
    }
  }
  
  // Overloading [] 
  // placemarker, should not be here
  // T& operator[] (const int i) { return data[i]; }

  // these give the field_element -- WILL BE modified by transformer
  field_element<T>& operator[] (const parity p) const;
  field_element<T>& operator[] (const parity_plus_direction p);
  //{ 
  //  return (field_element<T>) *this;
  //}

  // Overloading =  
  // it is left to "field_element<>" to decide whether the assignment is OK
  template <typename R>
  field<T>& operator= (const field<R>& rhs) {
    (*this)[ALL] = rhs[X];
    return *this;
  }
  
  template <typename R>
  field<T>& operator= (const R& d) {
    (*this)[ALL] = d;
    return *this;
  }
  
  // Do also move assignment
  field<T>& operator= (field<T> && rhs) {
    if (this != &rhs) {
      free();
      fs = rhs.fs;
      rhs.fs = nullptr;
    }
    return *this;
  }
  
  template <typename R>
  field<T>& operator+= (const field<R>& rhs) { (*this)[ALL] += rhs[X]; return *this;}
  template <typename R>
  field<T>& operator-= (const field<R>& rhs) { (*this)[ALL] -= rhs[X]; return *this;}
  template <typename R>
  field<T>& operator*= (const field<R>& rhs) { (*this)[ALL] *= rhs[X]; return *this;}
  template <typename R>
  field<T>& operator/= (const field<R>& rhs) { (*this)[ALL] /= rhs[X]; return *this;}

  template <typename R>
  field<T>& operator+= (const R & rhs) { (*this)[ALL] += rhs; return *this;}
  template <typename R>
  field<T>& operator-= (const R & rhs) { (*this)[ALL] -= rhs; return *this;}

  template <typename R>
  field<T>& operator*= (const R & rhs) { (*this)[ALL] *= rhs; return *this;}
  template <typename R>
  field<T>& operator/= (const R & rhs) { (*this)[ALL] /= rhs; return *this;}

  
};


template <typename T, typename L, typename R>
field<T> operator+( const field<L> &lhs, const field<R> &rhs) {
  field<T> tmp;
  tmp[ALL] = lhs[X] + rhs[X];
  return tmp;
}

template <typename T,typename L>
field<T> operator+( const L &lhs, const field<T> &rhs) {
  field<T> tmp;
  tmp[ALL] = lhs + rhs[X];
  return tmp;
}

template <typename T,typename R>
field<T> operator+( const field<T> &lhs,  const R &rhs) {
  return operator+(rhs,lhs);
}

// template <typename T>
// field<T> operator+( const double lhs, const field<T> &rhs) {
//   field<T> tmp;
//   tmp[ALL] = lhs + rhs[X];
//   return tmp;
// }

// template <typename T>
// field<T> operator+( const field<T> &lhs, const double rhs) {
  
//   return lhs;
// }

template <typename T, typename R>
field<T> operator*( const field<T> &lhs, const field<T> &rhs) {
  return lhs;
}



// template <typename T>
// class reduction {
// private:
//   field_element<T> value;
  
// public:
//   void sum(const field_element<T> & rhs) { value += rhs; }
//   void operator+=(const field_element<T> & rhs) { value += rhs; }
//   void product(const field_element<T> & rhs) { value *= rhs; }
//   void operator*=(const field_element<T> & rhs) { value *= rhs; }

//   // TODO - short vectors!
//   void max(const field_element<T> & rhs) { if (rhs > value) value = rhs; }
//   void min(const field_element<T> & rhs) { if (rhs < value) value = rhs; }

//   // TODO - short vectors!
//   T get() { return value.get_value(); }

//   // implicit conversion
//   operator T() { return get(); }
  
// };


#endif

