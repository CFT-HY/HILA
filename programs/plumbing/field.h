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


// This is a marker for transformer -- for does not survive as it is
#define onsites(p) for(parity parity_type_var_(p);;)


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

 public:
  enum class status : unsigned { NOT_DONE, STARTED, DONE };

 private:

  /// The following struct holds the data + information about the field
  /// TODO: field-specific boundary conditions?
  class field_struct {
    public:
      field_storage<T> payload; // TODO: must be maximally aligned, modifiers - never null
      lattice_struct * lattice;
      unsigned assigned_to;           // keeps track of first assignment to parities
      status move_status[3][NDIRS];     // is communication done
#ifdef USE_MPI
      MPI_Request receive_request[3][NDIRS];
      MPI_Request send_request[3][NDIRS];
#ifndef VANILLA
      // vanilla needs no special receive buffers
      char * receive_buffer[NDIRS];
#endif
      char * send_buffer[NDIRS];

      void initialize_communication(){
        for (int d=0; d<NDIRS; d++) {
          for (int p=0; p<3; p++) move_status[p][d] = status::NOT_DONE;
          send_buffer[d] = nullptr;
#ifndef VANILLA
          receive_buffer[d] = nullptr;
#endif
        }
      
        // for(int d=0; d<NDIRS; d++) for(parity par: {EVEN,ODD}) {
        //   int tag = d + NDIRS*(int)par;
        //   lattice_struct::nn_comminfo_struct ci = lattice->get_comminfo(d);
        //   receive_buffer[tag].resize(ci.from_node.size());
        //   send_buffer[tag].resize(ci.to_node.size());
        //   receive_request[tag].resize(ci.from_node.size());
        //   send_request[tag].resize(ci.to_node.size());
        // }
        // mpi_tag_base = get_field_mpi
        // next_mpi_field_tag++;
      }

      void free_communication() {
        for (int d=0; d<NDIRS; d++) {
          if (send_buffer[d] != nullptr) std::free(send_buffer[d]);
#ifndef VANILLA
          if (receive_buffer[d] != nullptr) std::free(receive_buffer[d]);
#endif
        }
      }
#else
      // empty stubs
      void initialize_communication(){}
      void free_communication(){}
#endif

      void allocate_payload() { 
        payload.allocate_field(lattice);
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
      void gather_comm_elements(char * RESTRICT buffer, const lattice_struct::comm_node_struct & to_node, parity par) const {
        payload.gather_comm_elements(buffer, to_node, par, lattice);
      };

      /// Place boundary elements from neighbour
      void place_comm_elements(char * RESTRICT buffer, const lattice_struct::comm_node_struct & from_node, parity par){
        payload.place_comm_elements(buffer, from_node, par, lattice);
      };
      
      /// Place boundary elements from local lattice (used in vectorized version)
      void set_local_boundary_elements(direction dir, parity par){
        payload.set_local_boundary_elements(dir, par, lattice);
      };
  };

  static_assert( std::is_pod<T>::value, "Field expects only pod-type elements (plain data): default constructor, copy and delete");
  
 public:

  field_struct * RESTRICT fs;
  
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
    fs->initialize_communication();
    mark_changed(ALL);      // guarantees communications will be done
    fs->assigned_to = 0;    // and this means that it is not assigned
  }

  void free() {
    if (fs != nullptr) {
      fs->free_payload();
      fs->free_communication();
      delete fs;
      fs = nullptr;
    }
  }

  bool is_allocated() const { return (fs != nullptr); }

  bool is_initialized(parity p) const { 
    return fs != nullptr && ((fs->assigned_to & parity_bits(p)) != 0);
  }
  
  status move_status(parity p, int d) const { 
    assert(parity_bits(p) && d>=0 && d<NDIRS);
    return fs->move_status[(int)p - 1][d]; 
  }
  void set_move_status(parity p, int d, status stat) const { 
    assert(parity_bits(p) && d>=0 && d<NDIRS);
    fs->move_status[(int)p - 1][d] = stat;
  }


  // If ALL changes, both parities invalid; if p != ALL, then p and ALL.
  void mark_changed_inner(const parity p) const {
    for (int i=0; i<NDIRS; i++) {
      set_move_status(p,i,status::NOT_DONE);
      if (p != ALL) set_move_status(ALL,i,status::NOT_DONE );
      else {
        set_move_status(EVEN,i,status::NOT_DONE);
        set_move_status(ODD,i,status::NOT_DONE);
      }
    }
    fs->assigned_to |= parity_bits(p);
  }

  // call this BEFORE the var is actually written to
  void mark_changed(const parity p) {
    if (fs == nullptr) allocate();
    mark_changed_inner(p);
  }

  // const version needed for fetches of const fields
  void mark_changed(const parity p) const {
    assert(is_allocated());
    mark_changed_inner(p);
  }
  
  /// Mark the field parity fetched from direction
  // In case p=ALL we could mark everything fetched, but we'll be conservative here 
  // and mark only this parity, because there might be other parities on the fly and corresponding
  // waits should be done,  This should never happen in automatically generated loops.
  // In any case start_get, is_fetched, get_move_parity has intelligence to figure out the right thing to do
  //

  void mark_fetched( int dir, const parity p) const {
    set_move_status(p,dir,status::DONE);
  }

  // Check if the field has been fetched since the previous communication
  // par = ALL:   ALL or (EVEN+ODD) are OK
  // par != ALL:  ALL or par are OK
  bool is_fetched( int dir, parity par) const {
    if (par != ALL) {
      return move_status(par,dir) == status::DONE || move_status(ALL,dir) == status::DONE;
    } else {
      return move_status(ALL,dir) == status::DONE ||
           ( move_status(EVEN,dir) == status::DONE && move_status(ODD,dir) == status::DONE );
    }
  }
   
  // Mark communication started -- this must be just the one
  // going on with MPI
  void mark_move_started( int dir, parity p) const{
    set_move_status(p,dir,status::STARTED);
  }

  /// Check if communication has started.  This is strict, checks exactly this parity
  bool is_move_started( int dir, parity par) const{
    return move_status(par,dir) == status::STARTED;
  }
    
  bool move_not_done( int dir, parity par) const {
    return move_status(par,dir) == status::NOT_DONE;
  }
 
  
  // Overloading [] 
  // placemarker, should not be here
  // T& operator[] (const int i) { return data[i]; }

  // declarations -- WILL BE implemented by transformer, not written here
  element<T>& operator[] (const parity p) const;             // f[EVEN]
  element<T>& operator[] (const X_index_type) const;         // f[X]
  element<T>& operator[] (const X_plus_direction p) const;   // f[X+dir]
  element<T>& operator[] (const X_plus_offset p) const;      // f[X+dir1+dir2] and others

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
  dir_mask_t start_get(direction d, parity p) const;
  dir_mask_t start_get(direction d) const { return start_get(d, ALL);}
  void wait_get(direction d, parity p) const;
  void get(direction d, parity p) const;

  // and declaration of shift methods
  field<T> shift(const coordinate_vector &v, parity par) const;
  field<T> shift(const coordinate_vector &v) const { return shift(v,ALL); }

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


#define NAIVE_SHIFT
#if defined(NAIVE_SHIFT)

// Define shift method here too - this is a placeholder, very inefficient
// works by repeatedly nn-copying the field

template<typename T>
field<T> field<T>::shift(const coordinate_vector &v, const parity par) const {
  field<T> r1, r2;
  r2 = *this;
  foralldir(d) {
    if (abs(v[d]) > 0) {
      direction dir;
      if (v[d] > 0) dir = d; else dir = -d;
    
      for (int i=0; i<abs(v[d]); i++) {
        r1[ALL] = r2[X+dir];
        r2 = r1;
      }
    }
  }
  return r2;
}

#elif !defined(USE_MPI)

template<typename T>
field<T> field<T>::shift(const coordinate_vector &v, const parity par) const {
  field<T> result;

  onsites(par) {
    if 
  }
  r2 = *this;
  foralldir(d) {
    if (abs(v[d]) > 0) {
      direction dir;
      if (v[d] > 0) dir = d; else dir = -d;
    
      for (int i=0; i<abs(v[d]); i++) {
        r1[ALL] = r2[X+dir];
        r2 = r1;
      }
    }
  }
  return r2;
}

#endif


#if defined(USE_MPI)
/* MPI implementations
 * For simplicity, these functions do not use field expressions andb
 * can be ignored by the transformer. Since the transformer does not
 * have access to mpi.h, it cannot process this branch.
 */

/// start_get(): Communicate the field at parity par from direction
///  d. Uses accessors to prevent dependency on the layout.
/// return the direction mask bits where something is happening
template<typename T>
dir_mask_t field<T>::start_get(direction d, parity p) const {

  // get the mpi message tag right away, to ensure that we are always synchronized with the
  // mpi calls -- some nodes might not need comms, but the tags must be in sync

  int tag = get_next_msg_tag();

  lattice_struct::nn_comminfo_struct  & ci = lattice->nn_comminfo[d];
  lattice_struct::comm_node_struct & from_node = ci.from_node;
  lattice_struct::comm_node_struct & to_node = ci.to_node;

  // check if this is done - either fetched or no comm to be done in the 1st place
  if (is_fetched(d,p) || (from_node.rank == mynode() && to_node.rank == mynode()) ) {
    lattice->n_gather_avoided++; 
    return 0;   // nothing to wait for
  }

  // if this parity or ALL-type fetch is going on nothing to be done
  if (!move_not_done(d,p) || !move_not_done(d,ALL)) {
    lattice->n_gather_avoided++;
    return get_dir_mask(d);     // nothing to do, but still need to wait 
  }

  parity par = p;
  // if p is ALL but ODD or EVEN is going on/done, turn off parity which is not needed
  // corresponding wait must do the same thing
  if (p == ALL) {
    if (!move_not_done(d,EVEN) && !move_not_done(d,ODD)) {
      // even and odd are going on or ready, nothing to be done
      lattice->n_gather_avoided++;
      return get_dir_mask(d);
    }
    if (!move_not_done(d,EVEN)) par = ODD;
    else if (!move_not_done(d,ODD)) par = EVEN;
    // if neither is the case par = ALL
  }

  mark_move_started(d, par);

  // Communication hasn't been started yet, do it now

  int par_i = static_cast<int>(par)-1;  // index to dim-3 arrays

  constexpr int size = sizeof(T);

  char * receive_buffer;
  char * send_buffer;

  if (from_node.rank != mynode()) {

    // HANDLE RECEIVES: get node which will send here

#ifdef VANILLA
    // in vanilla code the receive buffer is the field buffer, set offsets
    // field_buffer gives the right type, so addition gives the right offset without size
    receive_buffer = ((char *)field_buffer()) + from_node.offset(par) * size;
#else
    if (fs->receive_buffer[d] == nullptr) {
      fs->receive_buffer[d] = (char *)memalloc( size*from_node.sites);
    }
    receive_buffer = fs->receive_buffer[d] + from_node.offset(par) * size;
#endif
  
    unsigned sites = from_node.n_sites(par);

    // c++ version does not return errors??
    MPI_Irecv( receive_buffer, sites*size, MPI_BYTE, from_node.rank,
	             tag, lattice->mpi_comm_lat, &fs->receive_request[par_i][d] );
  }

  if (to_node.rank != mynode()) {
    // HANDLE SENDS: Copy field elements on the boundary to a send buffer and send
    unsigned sites = to_node.n_sites(par);

    if(fs->send_buffer[d] == nullptr)
      fs->send_buffer[d] = (char *)memalloc( to_node.sites*size );
    send_buffer = fs->send_buffer[d] + to_node.offset(par) * size;

    fs->gather_comm_elements(send_buffer, to_node, par);
 
    MPI_Isend( send_buffer, sites*size, MPI_BYTE, to_node.rank, 
               tag, lattice->mpi_comm_lat, &fs->send_request[par_i][d]);
  }


  return get_dir_mask(d);

}

///* wait_get(): Wait for communication at parity par from
///  direction d completes the communication in the function.
///  If the communication has not started yet, also calls
///  start_get()
///
///  NOTE: This will be called even if the field is marked const.
///  Therefore this function is const, even though it does change
///  the internal content of the field, the halo. From the point
///  of view of the user, the value of the field does not change.
template<typename T>
void field<T>::wait_get(direction d, parity p) const {

  lattice_struct::nn_comminfo_struct  & ci = lattice->nn_comminfo[d];
  lattice_struct::comm_node_struct & from_node = ci.from_node;
  lattice_struct::comm_node_struct & to_node = ci.to_node;

  // check if this is done - either fetched or no comm to be done in the 1st place
  if (is_fetched(d,p) || 
      (from_node.rank == mynode() && to_node.rank == mynode()) ) {
    return;
  }

  // if (!is_move_started(d,p)) {
  //   output0 << "Wait move error - wait_get without corresponding start_get\n";
  //   exit(-1);
  // }

  // Note: the move can be parity p OR ALL -- need to wait for it in any case
  // set par to be the "sum" over both parities
  // There never should be ongoing ALL and other parity fetch -- start_get takes care

  // check here consistency, this should never happen
  if (p != ALL && is_move_started(d,p) && is_move_started(d,ALL)) {
    output0 << "wait_get move parity error!\n";
    exit(-1);
  }

  parity par;
  int n_wait = 1;
  // what par to wait for?
  if (is_move_started(d,p)) par = p;            // standard match
  else if (p != ALL) {
    if (is_move_started(d,ALL)) par = ALL;      // if all is running wait for it
    else {
      output0 << "wait_get error: no matching wait found for parity " << (int)p << '\n';
      exit(-1);
    }
  } else {
    // now p == ALL and ALL is not running
    if (is_fetched(d,EVEN) && is_move_started(d,ODD)) par = ODD;
    else if (is_fetched(d,ODD) && is_move_started(d,EVEN)) par = EVEN;
    else if (is_move_started(d,EVEN) && is_move_started(d,ODD)) {
      n_wait = 2;  // need to wait for both! 
      par = ALL;  
    } else {
      output0 << "wait_get error: no matching wait found for parity ALL\n";
      exit(-1);
    }
  }

  // Update local elements in the halo (necessary for vectorized version)
  fs->set_local_boundary_elements(d, par);

  if (n_wait == 2) par = EVEN; // we'll flip both

  for (int wait_i = 0; wait_i < n_wait; ++wait_i ) {

    int par_i = (int)par - 1;

    if (from_node.rank != mynode()) {
      MPI_Status status;
      MPI_Wait( &fs->receive_request[par_i][d], &status);

#ifndef VANILLA
      fs->place_comm_elements( fs->receive_buffer[d], from_node, par);
#endif
    }

    // then wait for the sends
    if (to_node.rank != mynode()) {
      MPI_Status status;
      MPI_Wait( &fs->send_request[par_i][d], &status );
    }

    // Mark the parity fetched from direction dir
    mark_fetched(d, par);
  
    // Keep count of communications
    lattice->n_gather_done += 1;

    par = opp_parity(par);  // flip if 2 loops
  }
}



#else

///* Trivial implementation when no MPI is used
#include "../plumbing/comm_vanilla.h"
template<typename T>
dir_mask_t field<T>::start_get(direction d, parity p) const {
  // Update local elements in the halo (necessary for vectorized version)
  // Does not need to happen every time; should use tracking like in MPI
  fs->set_local_boundary_elements(d, p);
  return 0;
}

template<typename T>
void field<T>::wait_get(direction d, parity p) const {}


#endif  // MPI

/// And a conveniece combi function
template<typename T>
void field<T>::get(direction d, parity p) const {
  start_get(d,p);
  wait_get(d,p);
}



#endif // FIELD_H


