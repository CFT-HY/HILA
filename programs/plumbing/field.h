// -*- mode: c++ -*-
#ifndef FIELD_H
#define FIELD_H
#include <iostream>
#include <string>
#include <math.h>
#include <type_traits>

// using namespace std;

#include "../plumbing/lattice.h"


// HACK  -- this is needed for pragma handlin, do not change!
// #pragma transformer _transformer_cmd_dump_ast_

// HACK
#define transformer_ctl(a) extern int _transformer_ctl_##a
//void transformer_control(const char *);


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


// fwd definition
// template <typename T> class field;

#if 0
// field_element class: virtual class, no storage allocated,
// wiped out by the transformer
template <typename T>
class field_element  {
 private:
  T v;   // TODO: this must be set appropriately?
  
 public:
  // the following are just placemarkers, and are substituted by the transformer

  // implicit conversion to type T: this works for storage types where field is an std
  // array of type T - this is again for type propagation
  // TODO: IS THIS NEEDED?  WE WANT TO AVOID CONVERSIONS FROM field<T> v -> T
  // operator T() { return v; }
      
  // The type is important for ensuring correctness
  // Possibility: write these so that they work without the transformer
  template <typename A,
            std::enable_if_t<std::is_assignable<T&,A>::value, int> = 0 >
  field_element<T>& operator= (const A &d) {
    v = d; 
    return *this;
  }
  
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

#endif

template <typename T>
using field_element = T;


// Type alias would be nice, but one cannot specialize those!
// Placemarker, will be specialized by transformer
template <typename T>
struct field_storage_type {
  T c;
};

/// The following struct holds the data + information about the field
// template <typename T>
// class field_struct {
//   friend field<T>;
// 
//   field_storage_type<T> * payload; // TODO: must be maximally aligned, modifiers - never null
//   lattice_struct * lattice;
//   unsigned is_fetched[NDIRS];
//   
// };


// These are helpers, to make generic templates
// e.g. t_plus(A,B) gives the type of the operator a + b, where a is of type A and b B.
#define t_plus(A,B)  decltype(std::declval<A>() + std::declval<B>())
#define t_minus(A,B) decltype(std::declval<A>() - std::declval<B>())
#define t_mul(A,B)   decltype(std::declval<A>() * std::declval<B>())
#define t_div(A,B)   decltype(std::declval<A>() / std::declval<B>())


// ** class field

template <typename T>
class field {
private:

  static_assert( std::is_trivial<T>::value, "Field expects only trivial elements");
  
  /// TODO: field-specific boundary conditions?
  class field_struct {
  private:
    field_storage_type<T> * payload; // TODO: must be maximally aligned, modifiers - never null
    
  public:
    void allocate_payload() {
      payload = (field_storage_type<T> *) 
                   allocate_field_mem(sizeof(T) * lattice->mynode.field_alloc_size);
      if (payload == nullptr) {
        std::cout << "Failure in field memory allocation\n";
        exit(1);
      }
    }
    
    void free_payload() {
      free_field_mem((void *)payload);
      payload = nullptr;
    }

    unsigned is_fetched[NDIRS];
    lattice_struct * lattice;

    T& operator[] (int i)
    {
      return (T&) (payload[i].c);
    }
    
  };
  
  field_struct * fs;
    
public:
  
  field<T>() {
    // std::cout << "In constructor 1\n";
    fs = nullptr;             // lazy allocation on 1st use
  }
  
  // TODO: for some reason this straightforward copy constructor seems to be necessary, the
  // one below it does not catch implicit copying.  Try to understand why
  field<T>(const field<T>& other) {
    fs = nullptr;  // this is probably unnecessary
    (*this)[ALL] = other[X];
  }
    
  // copy constructor - from fields which can be assigned
  template <typename A,
            std::enable_if_t<std::is_convertible<A,T>::value, int> = 0 >  
  field<T>(const field<A>& other) {
    fs = nullptr;  // this is probably unnecessary
    (*this)[ALL] = other[X];
  }


  template <typename A,
            std::enable_if_t<std::is_convertible<A,T>::value, int> = 0 >  
  field<T>(const A& val) {
    fs = nullptr;
    // static_assert(!std::is_same<A,int>::value, "in int constructor");
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
    
  void allocate() {
    assert(fs == nullptr);
    if (lattice == nullptr) {
      // TODO: write to some named stream
      std::cout << "Can not allocate field variables before lattice.setup()\n";
      exit(1);  // TODO - more ordered exit?
    }
    fs = new field_struct;
    fs->lattice = lattice;
    fs->allocate_payload();

    mark_changed(ALL);
  }

  void free() {
    if (fs != nullptr) {
      fs->free_payload();
      delete fs;
      fs = nullptr;
    }
  }

  bool is_allocated() { return (fs != nullptr); }
  
  // call this BEFORE the var is written to
  void mark_changed(const parity p) {
    if (fs == nullptr) allocate();
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

  T& operator[] (int i)
  {
    return (T&) (this->fs[0][i]);
  }

  // fetch the element at this loc
  // T get(int i) const;
  
  // Overloading = - possible only if T = A is OK
  template <typename A, 
            std::enable_if_t<std::is_assignable<T&,A>::value, int> = 0 >
  field<T>& operator= (const field<A>& rhs) {
    (*this)[ALL] = rhs[X];
    return *this;
  }

  // same but without the field
  template <typename A, 
            std::enable_if_t<std::is_assignable<T&,A>::value, int> = 0 >
  field<T>& operator= (const A& d) {
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
  
  // is OK if T+A can be converted to type T
  template <typename A,
            std::enable_if_t<std::is_convertible<t_plus(T,A),T>::value, int> = 0>
  field<T>& operator+= (const field<A>& rhs) { 
    (*this)[ALL] += rhs[X]; return *this;
  }
  
  template <typename A,
            std::enable_if_t<std::is_convertible<t_minus(T,A),T>::value, int> = 0>  
  field<T>& operator-= (const field<A>& rhs) { 
    (*this)[ALL] -= rhs[X];
    return *this;
  }
  
  template <typename A,
            std::enable_if_t<std::is_convertible<t_mul(T,A),T>::value, int> = 0>
  field<T>& operator*= (const field<A>& rhs) {
    (*this)[ALL] *= rhs[X]; 
    return *this;
  }

  template <typename A,
            std::enable_if_t<std::is_convertible<t_div(T,A),T>::value, int> = 0>
  field<T>& operator/= (const field<A>& rhs) {
    (*this)[ALL] /= rhs[X];
    return *this;
  }

  template <typename A,
            std::enable_if_t<std::is_convertible<t_plus(T,A),T>::value, int> = 0>
  field<T>& operator+= (const A & rhs) { (*this)[ALL] += rhs; return *this;}

  template <typename A,
            std::enable_if_t<std::is_convertible<t_minus(T,A),T>::value, int> = 0>  
  field<T>& operator-= (const A & rhs) { (*this)[ALL] -= rhs; return *this;}

  template <typename A,
            std::enable_if_t<std::is_convertible<t_mul(T,A),T>::value, int> = 0>
  field<T>& operator*= (const A & rhs) { (*this)[ALL] *= rhs; return *this;}
  
  template <typename A,
            std::enable_if_t<std::is_convertible<t_div(T,A),T>::value, int> = 0>
  field<T>& operator/= (const A & rhs) { (*this)[ALL] /= rhs; return *this;}

  
};


// these operators rely on SFINAE, OK if field_t_plus(A,B) exists i.e. A+B is OK
/// operator +
template <typename A, typename B>
auto operator+( const field<A> &lhs, const field<B> &rhs) -> field<t_plus(A,B)>
{
  field <t_plus(A,B)> tmp;
  tmp[ALL] = lhs[X] + rhs[X];
  return tmp;
}

template <typename A, typename B>
auto operator+( const A &lhs, const field<B> &rhs) -> field<t_plus(A,B)>
{
  field<t_plus(A,B)> tmp;
  tmp[ALL] = lhs + rhs[X];
  return tmp;
}

template <typename A, typename B>
auto operator+( const field<A> &lhs, const B &rhs) -> field<t_plus(A,B)>
{
  field<t_plus(A,B)> tmp;
  tmp[ALL] = lhs[X] + rhs;
  return tmp;
}

/// operator -
template <typename A, typename B>
auto operator-( const field<A> &lhs, const field<B> &rhs) -> field<t_minus(A,B)>
{
  field <t_minus(A,B)> tmp;
  tmp[ALL] = lhs[X] - rhs[X];
  return tmp;
}

template <typename A, typename B>
auto operator-( const A &lhs, const field<B> &rhs) -> field<t_minus(A,B)>
{
  field<t_minus(A,B)> tmp;
  tmp[ALL] = lhs - rhs[X];
  return tmp;
}

template <typename A, typename B>
auto operator-( const field<A> &lhs, const B &rhs) -> field<t_minus(A,B)>
{
  field<t_minus(A,B)> tmp;
  tmp[ALL] = lhs[X] - rhs;
  return tmp;
}


/// operator *
template <typename A, typename B>
auto operator*( const field<A> &lhs, const field<B> &rhs) -> field<t_mul(A,B)>
{
  field <t_mul(A,B)> tmp;
  tmp[ALL] = lhs[X] * rhs[X];
  return tmp;
}

template <typename A, typename B>
auto operator*( const A &lhs, const field<B> &rhs) -> field<t_mul(A,B)>
{
  field<t_mul(A,B)> tmp;
  tmp[ALL] = lhs * rhs[X];
  return tmp;
}

template <typename A, typename B>
auto operator*( const field<A> &lhs, const B &rhs) -> field<t_mul(A,B)>
{
  field<t_mul(A,B)> tmp;
  tmp[ALL] = lhs[X] * rhs;
  return tmp;
}

/// operator /
template <typename A, typename B>
auto operator/( const field<A> &lhs, const field<B> &rhs) -> field<t_div(A,B)>
{
  field <t_div(A,B)> tmp;
  tmp[ALL] = lhs[X] / rhs[X];
  return tmp;
}

template <typename A, typename B>
auto operator/( const A &lhs, const field<B> &rhs) -> field<t_div(A,B)>
{
  field<t_div(A,B)> tmp;
  tmp[ALL] = lhs / rhs[X];
  return tmp;
}

template <typename A, typename B>
auto operator/( const field<A> &lhs, const B &rhs) -> field<t_div(A,B)>
{
  field<t_div(A,B)> tmp;
  tmp[ALL] = lhs[X] / rhs;
  return tmp;
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

//field<T> operator*( const field<T> &lhs, const field<T> &rhs) {
//  return lhs;
//}


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

