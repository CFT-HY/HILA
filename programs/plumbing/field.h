#ifndef FIELD_H
#define FIELD_H
#include <iostream>
#include <string>
#include <cstring> //Memcpy is here...
#include <math.h>
#include <type_traits>

#include "../plumbing/globals.h"
#include "../plumbing/defs.h"
#include "../plumbing/field_storage.h"
#include "../plumbing/lattice.h"

#ifdef USE_MPI
#include "../plumbing/comm_mpi.h"
#endif

static int next_mpi_field_tag = 0;


// This is a marker for transformer -- for does not survive as it is
#define onsites(p) for(parity parity_type_var_(p);;)

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

template <typename T>
field_element<T> operator*( const field_element<T> &lhs, const field_element<T> &rhs);

template <typename T,typename L>
field_element<T> operator*( const L &lhs, const field_element<T> &rhs);

template <typename T,typename R>
field_element<T> operator*( const field_element<T> &lhs,  const R &rhs);

template <typename T>
field_element<T> operator/( const field_element<T> &lhs, const field_element<T> &rhs);

template <typename T,typename L>
field_element<T> operator/( const L &lhs, const field_element<T> &rhs);

template <typename T,typename R>
field_element<T> operator/( const field_element<T> &lhs,  const R &rhs);

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
using element = T;

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

// field class 
template <typename T>
class field {
 private:

  /// The following struct holds the data + information about the field
  /// TODO: field-specific boundary conditions?
  class field_struct {
    public:
      field_storage<T> payload; // TODO: must be maximally aligned, modifiers - never null
      lattice_struct * lattice;
      bool is_fetched[2*NDIRS];
      bool move_started[2*NDIRS];
#ifdef USE_MPI
      std::vector<MPI_Request> receive_request[3*NDIRS];
      std::vector<MPI_Request> send_request[3*NDIRS];
      std::vector<char *> receive_buffer[3*NDIRS];
      std::vector<char *> send_buffer[3*NDIRS];
      int mpi_tag;
      void initialize_communication(){
        for(int d=0; d<NDIRS; d++) for(parity par: {EVEN,ODD}) {
          int tag = d + NDIRS*(int)par;
          lattice_struct::comminfo_struct ci = lattice->get_comminfo(d);
          receive_buffer[tag].resize(ci.from_node.size());
          send_buffer[tag].resize(ci.to_node.size());
          receive_request[tag].resize(ci.from_node.size());
          send_request[tag].resize(ci.to_node.size());
        }
        mpi_tag = next_mpi_field_tag;
        next_mpi_field_tag++;
      }
#else
      void initialize_communication(){};
#endif

      void allocate_payload() { 
        payload.allocate_field(lattice);
        initialize_communication();
      }
      void free_payload() { payload.free_field(); }

      /// Getter for an individual elements. Will not work in CUDA host code,
      /// but must be defined
      inline auto get(const int i) const {
        return payload.get( i, lattice->field_alloc_size() );
      }

      template<typename A>
      inline void set(const A & value, const int i) {
        payload.set( value, i, lattice->field_alloc_size() );
      }

      /// Gather boundary elements for communication
      void gather_comm_elements(char * buffer, lattice_struct::comm_node_struct to_node, parity par) const {
        payload.gather_comm_elements(buffer, to_node, par, lattice);
      };

      /// Place boundary elements from neighbour
      void place_comm_elements(char * buffer, lattice_struct::comm_node_struct from_node, parity par){
        payload.place_comm_elements(buffer, from_node, par, lattice);
      };
      
      /// Place boundary elements from local lattice (used in vectorized version)
      void set_local_boundary_elements(direction dir, parity par){
        payload.set_local_boundary_elements(dir, par, lattice);
      };
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
      assert(p == EVEN || p == ODD || p == ALL);
      if(p==EVEN) for (int i=0; i<2*NDIRS; i+=2) fs->is_fetched[i] = false;
      if(p==ODD)  for (int i=1; i<2*NDIRS; i+=2) fs->is_fetched[i] = false;
      if(p==ALL)  for (int i=0; i<2*NDIRS; i+=1) fs->is_fetched[i] = false;
      if(p==EVEN) for (int i=0; i<2*NDIRS; i+=2) fs->move_started[i] = false;
      if(p==ODD)  for (int i=1; i<2*NDIRS; i+=2) fs->move_started[i] = false;
      if(p==ALL)  for (int i=0; i<2*NDIRS; i+=1) fs->move_started[i] = false;
    }
  }

  void mark_changed(const parity p) const {
    assert(is_allocated());
    assert(p == EVEN || p == ODD || p == ALL);
    if(p==EVEN) for (int i=0; i<2*NDIRS; i+=2) fs->is_fetched[i] = false;
    if(p==ODD)  for (int i=1; i<2*NDIRS; i+=2) fs->is_fetched[i] = false;
    if(p==ALL)  for (int i=0; i<2*NDIRS; i+=1) fs->is_fetched[i] = false;
    if(p==EVEN) for (int i=0; i<2*NDIRS; i+=2) fs->move_started[i] = false;
    if(p==ODD)  for (int i=1; i<2*NDIRS; i+=2) fs->move_started[i] = false;
    if(p==ALL)  for (int i=0; i<2*NDIRS; i+=1) fs->move_started[i] = false;
  }

  /// Mark the field parity fetched from direction
  void mark_fetched(int dir, const parity p) const {
    assert(p == EVEN || p == ODD || p == ALL);
    if(p==EVEN || p==ALL) fs->is_fetched[2*dir] = true;
    if(p==ODD || p==ALL ) fs->is_fetched[2*dir+1] = true;
    
  }

  /// Check if the field has been changed since the previous communication
  bool is_fetched( int dir, parity par) const{
    assert(dir < NDIRS);
    bool * fd = fs->is_fetched + 2*dir;
    return (par==ALL && fd[0] && fd[1] ) || (par==EVEN && fd[0]) || (par==ODD && fd[1]);
  }

  /* Mark communication started */
  void mark_move_started( int dir, parity p) const{
    assert(dir < NDIRS);
    if(p==EVEN || p==ALL) fs->move_started[2*dir] = true;
    if(p==ODD  || p==ALL) fs->move_started[2*dir+1] = true;
  }

  /// Check if communication has started
  bool is_move_started( int dir, parity par) const{
    assert(dir < NDIRS);
    bool * sd = fs->move_started + 2*dir;
    return (par==ALL && sd[0] && sd[1] ) || (par==EVEN && sd[0]) || (par==ODD && sd[1]);
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

  // declarations -- WILL BE implemented by transformer, not written here
  element<T>& operator[] (const parity p) const;
  element<T>& operator[] (const parity_plus_direction p) const;
  element<T>& operator[] (const parity_plus_offset p) const;

  #if defined(VANILLA)
  // TEMPORARY HACK: return ptr to bare array
  inline auto field_buffer() const { return this->fs->payload.get_buffer(); }
  #endif

  /// Get an individual element outside a loop. This is also used as a getter in the vanilla code.
  inline auto get_value_at(int i) const { return this->fs->get(i); }

  /// Set an individual element outside a loop. This is also used as a setter in the vanilla code.
  template<typename A>
  inline void set_value_at(const A & value, int i) { this->fs->set( value, i); }

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
  void start_move(direction d) const {start_move(d, ALL);}
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



#if defined(USE_MPI) && !defined(TRANSFORMER) 
/* MPI implementations
 * For simplicity, these functions do not use field expressions and
 * can be ignored by the transformer. Since the transformer does not
 * have access to mpi.h, it cannot process this branch.
 */

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
    int index = d + NDIRS*(int)par;
    int tag =  fs->mpi_tag*3*NDIRS + index;
    constexpr int size = sizeof(T);

    lattice_struct::comminfo_struct ci = lattice->comminfo[d];
    int n = 0;
    std::vector<MPI_Request> & receive_request = fs->receive_request[index];
    std::vector<MPI_Request> & send_request = fs->send_request[index];
    std::vector<char *> & receive_buffer = fs->receive_buffer[index];
    std::vector<char *> & send_buffer = fs->send_buffer[index];

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
       if(send_buffer[n] == NULL)
         send_buffer[n] = (char *)malloc( sites*size );

       fs->gather_comm_elements(send_buffer[n], to_node, par);
 
       //printf("node %d, send tag %d to %d\n", mynode(), tag, to_node.index);
       /* And send */
       MPI_Isend( send_buffer[n], sites*size, MPI_BYTE, to_node.index, 
               tag, lattice->mpi_comm_lat, &send_request[n]);
       //printf("node %d, sent tag %d\n", mynode(), tag);
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
///  the internal content of the field, the halo. From the point
///  of view of the user, the value of the field does not change.
template<typename T>
void field<T>::wait_move(direction d, parity p) const {

  // Loop over parities
  // (if p=ALL, do both EVEN and ODD otherwise just p);
  for( parity par: loop_parities(p) ) {
    int index = d + NDIRS*(int)par;
    int tag =  fs->mpi_tag*3*NDIRS + index;

    if( is_fetched(d, par) ){
      // Not changed, return directly
      // Keep count of gathers optimized away
      lattice->n_gather_avoided += 1;
      return;
    }

    // This will start the communication if it has not been started yet
    start_move(d, par);

    // Update local elements in the halo (necessary for vectorized version)
    fs->set_local_boundary_elements(d, par);

    lattice_struct::comminfo_struct ci = lattice->comminfo[d];
    std::vector<MPI_Request> & receive_request = fs->receive_request[index];
    std::vector<char *> & receive_buffer = fs->receive_buffer[index];

    /* Wait for the data here */
    int n = 0;
    for( lattice_struct::comm_node_struct from_node : ci.from_node ){
      MPI_Status status;
      //printf("node %d, waiting for recv tag %d\n", mynode(), tag);
      MPI_Wait(&receive_request[n], &status);
      //printf("node %d, received tag %d\n", mynode(), tag);

      fs->place_comm_elements(receive_buffer[n], from_node, par);
      n++;
    }

    /* Mark the parity fetched from direction dir */
    mark_fetched(d, par);

    /* Keep count of communications */
    lattice->n_gather_done += 1;
  }
}


#else

///* Trivial implementation when no MPI is used
#include "../plumbing/comm_vanilla.h"
template<typename T>
void field<T>::start_move(direction d, parity p) const {}
template<typename T>
void field<T>::wait_move(direction d, parity p) const {
  // Update local elements in the halo (necessary for vectorized version)
  // Does not need to happen every time; should use tracking like in MPI
  fs->set_local_boundary_elements(d, p);
}


#endif


#endif // FIELD_H


