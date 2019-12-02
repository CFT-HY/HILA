// -*- mode: c++ -*-
#ifndef FIELD_H
#define FIELD_H
#include <iostream>
#include <string>
#include <cstring> //Memcpy is here...
#include <math.h>
#include <type_traits>

#include "../plumbing/globals.h"

#ifdef USE_MPI
#include "../plumbing/comm_mpi.h"
#endif

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




// Pointer to field data and accessors. Only this is passed to the
// CUDA kernels and other accelerators and it only contains a minimal
// amount of data.
template <typename T>
class field_storage {
  public:
    #ifndef layout_SOA
      // Array of structures implementation
    T * fieldbuf;

    void allocate_field( const int field_alloc_size ) {
      fieldbuf = (T*) allocate_field_mem( sizeof(T) * field_alloc_size);
      #pragma acc enter data create(fieldbuf)
    }

    void free_field() {
      #pragma acc exit data delete(fieldbuf)
      free_field_mem((void *)fieldbuf);
      fieldbuf = nullptr;
    }

    #pragma transformer loop_function
    T get(const int i, const int field_alloc_size) const
    {
      return ((T *) fieldbuf)[i];
    }

    #pragma transformer loop_function
    void set(T value, const int i, const int field_alloc_size) 
    {
      ((T *) fieldbuf)[i] = value;
    }


    #else
    // Structure of arrays implementation
    constexpr static int t_elements = sizeof(T) / sizeof(real_t);
    real_t * fieldbuf;

    void allocate_field( const int field_alloc_size ) {
      fieldbuf = (real_t*) allocate_field_mem( t_elements*sizeof(real_t) * field_alloc_size );
      #pragma acc enter data create(fieldbuf)
    }

    void free_field() {
      #pragma acc exit data delete(fieldbuf)
      free_field_mem((void *)fieldbuf);
      fieldbuf = nullptr;
    }

    /// Get a single element in a field
    /// With CUDA this only works in a loop
    #pragma transformer loop_function
    inline T get(const int idx, const int field_alloc_size) const {
      assert( idx < field_alloc_size);
      T value;
      real_t *value_f = static_cast<real_t *>(static_cast<void *>(&value));
      for (int i=0; i<(sizeof(T)/sizeof(real_t)); i++) {
         value_f[i] = fieldbuf[i*field_alloc_size + idx];
      }
      return value; 
    }

    /// Set a single element in a field
    /// With CUDA this only works in a loop  
    #pragma transformer loop_function
    inline void set(T value, const int idx, const int field_alloc_size) {
      assert( idx < field_alloc_size);
      real_t *value_f = static_cast<real_t *>(static_cast<void *>(&value));
      for (int i=0; i<(sizeof(T)/sizeof(real_t)); i++) {
        fieldbuf[i*field_alloc_size + idx] = value_f[i];
      }
    }
    #endif
};





// These are helpers, to make generic templates
// e.g. t_plus<A,B> gives the type of the operator a + b, where a is of type A and b B.
template<typename A, typename B>
using t_plus = decltype(std::declval<A>() + std::declval<B>());
template<typename A, typename B>
using t_minus= decltype(std::declval<A>() - std::declval<B>());
template<typename A, typename B>
using t_mul  = decltype(std::declval<A>() * std::declval<B>());
template<typename A, typename B>
using t_div  = decltype(std::declval<A>() / std::declval<B>());




// ** class field

template <typename T>
class field {
private:

  /// The following struct holds the data + information about the field
  /// TODO: field-specific boundary conditions?
  class field_struct {
    private:
      constexpr static int t_elements = sizeof(T) / sizeof(real_t);
    public:
      field_storage<T> payload; // TODO: must be maximally aligned, modifiers - never null
      lattice_struct * lattice;
      unsigned is_fetched[NDIRS];
      bool move_started[3*NDIRS];
#ifdef USE_MPI
      std::vector<MPI_Request> receive_request[3*NDIRS];
      std::vector<MPI_Request> send_request[3*NDIRS];
      std::vector<char *> receive_buffer[3*NDIRS];
      std::vector<char *> send_buffer[3*NDIRS];
      void set_buffer_sizes(){
        for(int d=0; d<NDIRS; d++) for(parity par: {EVEN,ODD}) {
          int tag = d + NDIRS*(int)par;
          lattice_struct::comminfo_struct ci = lattice->comminfo[d];
          receive_buffer[tag].resize(ci.from_node.size());
          send_buffer[tag].resize(ci.to_node.size());
          receive_request[tag].resize(ci.from_node.size());
          send_request[tag].resize(ci.to_node.size());
        }
      }
#else
      void set_buffer_sizes(){};
#endif

      void allocate_payload() { 
        payload.allocate_field(lattice->field_alloc_size());
        set_buffer_sizes();
      }
      void free_payload() { payload.free_field(); }
      
      /// Getter for an individual elements. Will not work in CUDA host code,
      /// but must be defined
      T get(const int i) const { return  payload.get( i, lattice->field_alloc_size() ); }
      /// Getter for an individual elements. Will not work in CUDA host code,
      /// but must be defined
      void set(T value, const int i) { payload.set( value, i, lattice->field_alloc_size() ); }

      /// Gather boundary elements for communication
      void gather_comm_elements(char * buffer, lattice_struct::comm_node_struct to_node, parity par) const;
      /// Place boundary elements from neighbour
      void scatter_comm_elements(char * buffer, lattice_struct::comm_node_struct from_node, parity par);
  };

  static_assert( std::is_trivial<T>::value, "Field expects only trivial elements");
  
public:

  field_struct * fs;
  
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

  bool is_allocated() const { return (fs != nullptr); }
  
  /// call this BEFORE the var is written to
  void mark_changed(const parity p) {
    if (fs == nullptr) allocate();
    else {
      char pc = static_cast<char>(p);
      assert(p == EVEN || p == ODD || p == ALL);
      unsigned up = 0x3 & (!(static_cast<unsigned>(opp_parity(p))));
      for (int i=0; i<NDIRS; i++) fs->is_fetched[i] &= up;
      for (int i=0; i<3*NDIRS; i++) fs->move_started[i] = false;
    }
  }

  void mark_changed(const parity p) const {
    assert(is_allocated());
    char pc = static_cast<char>(p);
    assert(p == EVEN || p == ODD || p == ALL);
    unsigned up = 0x3 & (!(static_cast<unsigned>(opp_parity(p))));
    for (int i=0; i<NDIRS; i++) fs->is_fetched[i] &= up;
    for (int i=0; i<3*NDIRS; i++) fs->move_started[i] = false;
  }

  /// Mark the field parity fetched from direction
  void mark_fetched(int dir, const parity p) const {
    char pc = static_cast<char>(p);
    assert(p == EVEN || p == ODD || p == ALL);
    unsigned up = 0x3 & (static_cast<unsigned>(opp_parity(p)));
    fs->is_fetched[dir] |= up;
  }

  /// Check if the field has been changed since the previous communication
  bool is_fetched( int dir, parity par) const{
    assert(dir < NDIRS);
    int par_int = static_cast<unsigned>(opp_parity(par));
    return (fs->is_fetched[dir] ^ par_int) == 0;
  }

  /// Check if communication has started
  bool is_move_started( int dir, parity par) const{
    assert(dir < NDIRS);
    return fs->move_started[(int)par + 2*dir];
  }

  /* Mark communication started */
  void mark_move_started( int dir, parity par) const{
    assert(dir < NDIRS);
    fs->move_started[(int)par + 2*dir] = true;
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
  field_element<T>& operator[] (const parity_plus_direction p) const;
  //{ 
  //  return (field_element<T>) *this;
  //}


  /// Get an individual element outside a loop. This is also used as a getter in the vanilla code.
  T get_value_at(int i) const { return this->fs->get(i); }

  /// Set an individual element outside a loop. This is also used as a setter in the vanilla code.
  void set_value_at(T value, int i) { this->fs->set(value, i); }


  // fetch the element at this loc
  // T get(int i) const;
  
  // NOTE: THIS SHOULD BE INCLUDED IN TEMPLATE BELOW; SEEMS NOT???  
  field<T>& operator= (const field<T>& rhs) {
   (*this)[ALL] = rhs[X];
   return *this;
  }

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
            std::enable_if_t<std::is_convertible<t_plus<T,A>,T>::value, int> = 0>
  field<T>& operator+= (const field<A>& rhs) { 
    (*this)[ALL] += rhs[X]; return *this;
  }
  
  template <typename A,
            std::enable_if_t<std::is_convertible<t_minus<T,A>,T>::value, int> = 0>  
  field<T>& operator-= (const field<A>& rhs) { 
    (*this)[ALL] -= rhs[X];
    return *this;
  }
  
  template <typename A,
            std::enable_if_t<std::is_convertible<t_mul<T,A>,T>::value, int> = 0>
  field<T>& operator*= (const field<A>& rhs) {
    (*this)[ALL] *= rhs[X]; 
    return *this;
  }

  template <typename A,
            std::enable_if_t<std::is_convertible<t_div<T,A>,T>::value, int> = 0>
  field<T>& operator/= (const field<A>& rhs) {
    (*this)[ALL] /= rhs[X];
    return *this;
  }

  template <typename A,
            std::enable_if_t<std::is_convertible<t_plus<T,A>,T>::value, int> = 0>
  field<T>& operator+= (const A & rhs) { (*this)[ALL] += rhs; return *this;}

  template <typename A,
            std::enable_if_t<std::is_convertible<t_minus<T,A>,T>::value, int> = 0>  
  field<T>& operator-= (const A & rhs) { (*this)[ALL] -= rhs; return *this;}

  template <typename A,
            std::enable_if_t<std::is_convertible<t_mul<T,A>,T>::value, int> = 0>
  field<T>& operator*= (const A & rhs) { (*this)[ALL] *= rhs; return *this;}
  
  template <typename A,
            std::enable_if_t<std::is_convertible<t_div<T,A>,T>::value, int> = 0>
  field<T>& operator/= (const A & rhs) { (*this)[ALL] /= rhs; return *this;}


  // Communication routines
  void start_move(direction d, parity p) const;
  void wait_move(direction d, parity p) const;
};


// these operators rely on SFINAE, OK if field_t_plus<A,B> exists i.e. A+B is OK
/// operator +
template <typename A, typename B>
auto operator+( field<A> &lhs, field<B> &rhs) -> field<t_plus<A,B>>
{
  field <t_plus<A,B>> tmp;
  tmp[ALL] = lhs[X] + rhs[X];
  return tmp;
}

template <typename A, typename B>
auto operator+( const A &lhs, const field<B> &rhs) -> field<t_plus<A,B>>
{
  field<t_plus<A,B>> tmp;
  tmp[ALL] = lhs + rhs[X];
  return tmp;
}

template <typename A, typename B>
auto operator+( const field<A> &lhs, const B &rhs) -> field<t_plus<A,B>>
{
  field<t_plus<A,B>> tmp;
  tmp[ALL] = lhs[X] + rhs;
  return tmp;
}

/// operator -
template <typename A, typename B>
auto operator-( const field<A> &lhs, const field<B> &rhs) -> field<t_minus<A,B>>
{
  field <t_minus<A,B>> tmp;
  tmp[ALL] = lhs[X] - rhs[X];
  return tmp;
}

template <typename A, typename B>
auto operator-( const A &lhs, const field<B> &rhs) -> field<t_minus<A,B>>
{
  field<t_minus<A,B>> tmp;
  tmp[ALL] = lhs - rhs[X];
  return tmp;
}

template <typename A, typename B>
auto operator-( const field<A> &lhs, const B &rhs) -> field<t_minus<A,B>>
{
  field<t_minus<A,B>> tmp;
  tmp[ALL] = lhs[X] - rhs;
  return tmp;
}


/// operator *
template <typename A, typename B>
auto operator*( const field<A> &lhs, const field<B> &rhs) -> field<t_mul<A,B>>
{
  field <t_mul<A,B>> tmp;
  tmp[ALL] = lhs[X] * rhs[X];
  return tmp;
}

template <typename A, typename B>
auto operator*( const A &lhs, const field<B> &rhs) -> field<t_mul<A,B>>
{
  field<t_mul<A,B>> tmp;
  tmp[ALL] = lhs * rhs[X];
  return tmp;
}

template <typename A, typename B>
auto operator*( const field<A> &lhs, const B &rhs) -> field<t_mul<A,B>>
{
  field<t_mul<A,B>> tmp;
  tmp[ALL] = lhs[X] * rhs;
  return tmp;
}

/// operator /
template <typename A, typename B>
auto operator/( const field<A> &lhs, const field<B> &rhs) -> field<t_div<A,B>>
{
  field <t_div<A,B>> tmp;
  tmp[ALL] = lhs[X] / rhs[X];
  return tmp;
}

template <typename A, typename B>
auto operator/( const A &lhs, const field<B> &rhs) -> field<t_div<A,B>>
{
  field<t_div<A,B>> tmp;
  tmp[ALL] = lhs / rhs[X];
  return tmp;
}

template <typename A, typename B>
auto operator/( const field<A> &lhs, const B &rhs) -> field<t_div<A,B>>
{
  field<t_div<A,B>> tmp;
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







/* Communication routines for fields */
#ifndef USE_MPI
/* No MPI, trivial implementations */

#include "../plumbing/comm_vanilla.h"

///* start_move(): Trivial implementation when no MPI is used
template<typename T>
void field<T>::start_move(direction d, parity p) const {}
template<typename T>
void field<T>::wait_move(direction d, parity p) const {}


#else
/* MPI implementations
 * For simplicity, these functions do not use field expressions and
 * can be ignored by the transformer. Since the transformer does not
 * have access to mpi.h, it cannot process this branch.
 */


#if !defined(CUDA) || defined(TRANSFORMER)
/// The standard (Non-CUDA) implementation of gather_comm_elements
/* Gathers sites at the boundary that need to be communicated to neighbours */
template<typename T>
void field<T>::field_struct::gather_comm_elements(char * buffer, lattice_struct::comm_node_struct to_node, parity par) const {
  for (int j=0; j<to_node.n_sites(par); j++) {
    T element = get(to_node.site_index(j, par));
    std::memcpy( buffer + j*sizeof(T), (char *) (&element), sizeof(T) );
  }
}
/// The standard (Non-CUDA) implementation of scatter_comm_elements
/* Sets the values the neighbour elements from the communication buffer */
template<typename T>
void field<T>::field_struct::scatter_comm_elements(char * buffer, lattice_struct::comm_node_struct from_node, parity par){
  for (int j=0; j<from_node.n_sites(par); j++) {
    T element = *((T*) ( buffer + j*sizeof(T) ));
    set(element, from_node.offset(par)+j);
  }
}

#else
/* CUDA implementations */

/// A kernel that gathers neighbour elements for communication (using the getter)
template <typename T>
__global__ void gather_comm_elements_kernel( field_storage<T> field, char *buffer, int * site_index, const int sites, const int field_alloc_size )
{
  int Index = threadIdx.x + blockIdx.x * blockDim.x;
  if( Index < sites ) {
    ((T*) buffer)[Index] = field.get(site_index[Index], field_alloc_size);
  }
}

/// CUDA implementation of gather_comm_elements without CUDA aware MPI
/// Gathers sites at the boundary that need to be communicated to neighbours
template<typename T>
void field<T>::field_struct::gather_comm_elements(char * buffer, lattice_struct::comm_node_struct to_node, parity par) const {
  int *site_index, *d_site_index;
  char * d_buffer;
  int sites = to_node.n_sites(par);
  
  // Copy the list of boundary site indexes to the device
  site_index = (int *)std::malloc( sites*sizeof(int) );
  cudaMalloc( (void **)&(d_site_index), sites*sizeof(int));
  for (int j=0; j<sites; j++) {
    site_index[j] = to_node.site_index(j, par);
  }
  cudaMemcpy( d_site_index, site_index, sites*sizeof(int), cudaMemcpyHostToDevice );
  std::free(site_index);

  // Call the kernel to build the list of elements
  cudaMalloc( (void **)&(d_buffer), sites*sizeof(T));
  int N_blocks = sites/N_threads + 1; 
  gather_comm_elements_kernel<<< N_blocks, N_threads >>>( payload, d_buffer, d_site_index, sites, lattice->field_alloc_size() );
  
  // Copy the result to the host
  cudaMemcpy( buffer, d_buffer, sites*sizeof(T), cudaMemcpyDeviceToHost );

  cudaFree(d_site_index);
  cudaFree(d_buffer);
}


/// A kernel that scatters the neighbour sites received from a neihbour into 
/// it's proper place (using the setter)
template <typename T>
__global__ void scatter_comm_elements_kernel( field_storage<T> field, char *buffer, const int offset, const int sites, const int field_alloc_size )
{
  int Index = threadIdx.x + blockIdx.x * blockDim.x;
  if( Index < sites ) {
    field.set( ((T*) buffer)[Index], offset + Index, field_alloc_size);
  }
}

/// CUDA implementation of gather_comm_elements without CUDA aware MPI
/// Sets the values the neighbour elements from the communication buffer 
template<typename T>
void field<T>::field_struct::scatter_comm_elements(char * buffer, lattice_struct::comm_node_struct from_node, parity par){
  char * d_buffer;
  int sites = from_node.n_sites(par);

  // Allocate space and copy the buffer to the device
  cudaMalloc( (void **)&(d_buffer), sites*sizeof(T));
  cudaMemcpy( d_buffer, buffer, sites*sizeof(T), cudaMemcpyHostToDevice );

  // Call the kernel to place the elements 
  int N_blocks = sites/N_threads + 1; 
  scatter_comm_elements_kernel<<< N_blocks, N_threads >>>( payload, d_buffer, from_node.offset(par), sites, lattice->field_alloc_size() );

  cudaFree(d_buffer);
}
#endif



/// wait_move(): Communicate the field at parity par from direction
///  d. Uses accessors to prevent dependency on the layout.
template<typename T>
void field<T>::start_move(direction d, parity p) const {

  for( parity par: loop_parities(p) ) {
    if( is_move_started(d, par) ){
      // Not changed, return directly
      // Keep count of gathers optimized away
      lattice->n_gather_avoided += 1;
      return;
    }

    // Communication hasn't been started yet, do it now
    int tag = d + NDIRS*(int)par;
    constexpr int size = sizeof(T);

    lattice_struct::comminfo_struct ci = lattice->comminfo[d];
    int n = 0;
    std::vector<MPI_Request> & receive_request = fs->receive_request[tag];
    std::vector<MPI_Request> & send_request = fs->send_request[tag];
    std::vector<char *> & receive_buffer = fs->receive_buffer[tag];
    std::vector<char *> & send_buffer = fs->send_buffer[tag];

    /* HANDLE RECEIVES: loop over nodes which will send here */
    for( lattice_struct::comm_node_struct from_node : ci.from_node ){
      unsigned sites = from_node.n_sites(par);
      if(receive_buffer[n] == NULL)
        receive_buffer[n] = (char *)malloc( sites*size );

      //printf("node %d, recv tag %d from %d\n", mynode(), tag, from_node.index);

      MPI_Irecv( receive_buffer[n], sites*size, MPI_BYTE, from_node.index, 
	             tag, lattice->mpi_comm_lat, &receive_request[n] );
      n++;
    }

    /* HANDLE SENDS: Copy field elements on the boundary to a send buffer and send */
    n=0;
    for( lattice_struct::comm_node_struct to_node : ci.to_node ){
      /* gather data into the buffer  */
      unsigned sites = to_node.n_sites(par);
      send_buffer[n] = (char *)malloc( sites*size );
      fs->gather_comm_elements(send_buffer[n], to_node, par);

      //printf("node %d, send tag %d to %d\n", mynode(), tag, to_node.index);
      /* And send */
      MPI_Send( send_buffer[n], sites*size, MPI_BYTE, to_node.index, 
               tag, lattice->mpi_comm_lat);
      //printf("node %d, sent tag %d\n", mynode(), tag);
      std::free(send_buffer[n]);
      n++;
    }

    mark_move_started(d, par);
  }
}


///* wait_move(): Wait for communication at parity par from
///  direction d completes the communication in the function.
///  If the communication has not started yet, also calls
///  start_move()
///
///  NOTE: This will be called even if the field is marked const.
///  Therefore this function is const, even though it does change
///  the internal content of the field at the boundaries. From the 
///  point of view of the user, the value of the field does not change.
template<typename T>
void field<T>::wait_move(direction d, parity p) const {

  // Loop over parities
  // (if p=ALL, do both EVEN and ODD otherwise just p);
  for( parity par: loop_parities(p) ) {
    int tag = d + NDIRS*(int)par;

    if( is_fetched(d, par) ){
      // Not changed, return directly
      // Keep count of gathers optimized away
      lattice->n_gather_avoided += 1;
      return;
    }

    // This will start the communication if it has not been started yet
    start_move(d, par);

    lattice_struct::comminfo_struct ci = lattice->comminfo[d];
    std::vector<MPI_Request> & receive_request = fs->receive_request[tag];
    std::vector<char *> & receive_buffer = fs->receive_buffer[tag];

    /* Wait for the data here */
    int n = 0;
    for( lattice_struct::comm_node_struct from_node : ci.from_node ){
      MPI_Status status;
      //printf("node %d, waiting for recv tag %d\n", mynode(), tag);
      MPI_Wait(&receive_request[n], &status);
      //printf("node %d, received tag %d\n", mynode(), tag);

      fs->scatter_comm_elements(receive_buffer[n], from_node, par);
      n++;
    }

    /* Mark the parity fetched from direction dir */
    mark_fetched(d, par);

    /* Keep count of communications */
    lattice->n_gather_done += 1;
  }
}

#endif



#endif


