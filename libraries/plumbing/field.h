#ifndef FIELD_H
#define FIELD_H
#include <sstream>
#include<iostream>
#include <string>
#include <cstring> //Memcpy is here...
#include <math.h>
#include <type_traits>


#include "plumbing/defs.h"
#include "plumbing/field_storage.h"
#include "plumbing/lattice.h"

#include "plumbing/backend_vector/vector_types.h"


#ifdef USE_MPI
#include "plumbing/com_mpi.h"
#endif


// This is a marker for hilapp -- for does not survive as it is
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
  enum class fetch_status : unsigned { NOT_DONE, STARTED, DONE };

 private:


  /// The following struct holds the data + information about the field
  /// TODO: field-specific boundary conditions?
  class field_struct {
    public:
      field_storage<T> payload; // TODO: must be maximally aligned, modifiers - never null
      lattice_struct * lattice;
#ifdef VECTORIZED
      // get a direct ptr from here too, ease access
      vectorized_lattice_struct< vector_info<T>::vector_size > * vector_lattice;
#endif
      unsigned assigned_to;           // keeps track of first assignment to parities
      fetch_status move_status[3][NDIRS];     // is communication done

      // neighbour pointers - because of boundary conditions, can be different for diff. fields
      const unsigned * RESTRICT neighbours[NDIRS];
      boundary_condition_t boundary_condition[NDIRS];

#ifdef USE_MPI
      MPI_Request receive_request[3][NDIRS];
      MPI_Request send_request[3][NDIRS];
#ifndef VANILLA
      // vanilla needs no special receive buffers
      T * receive_buffer[NDIRS];
#endif
      T * send_buffer[NDIRS];


      void initialize_communication(){
        for (int d=0; d<NDIRS; d++) {
          for (int p=0; p<3; p++) move_status[p][d] = fetch_status::NOT_DONE;
          send_buffer[d] = nullptr;
#ifndef VANILLA
          receive_buffer[d] = nullptr;
#endif
        }
      }

      void free_communication() {
        for (int d=0; d<NDIRS; d++) {
          if (send_buffer[d] != nullptr) payload.free_mpi_buffer(send_buffer[d]);
#ifndef VANILLA
          if (receive_buffer[d] != nullptr) payload.free_mpi_buffer(receive_buffer[d]);
#endif
        }
      }

#else  // Not MPI

      // empty stubs
      void initialize_communication(){}
      void free_communication(){}

#endif

      void allocate_payload() { 
        payload.allocate_field(lattice);
      }
      void free_payload() { payload.free_field(); }

#ifndef VECTORIZED
      /// Getter for an individual elements in a loop
      inline auto get(const int i) const {
        return payload.get( i, lattice->field_alloc_size() );
      }

      template<typename A>
      inline void set(const A & value, const int i) {
        payload.set( value, i, lattice->field_alloc_size() );
      }

      /// Getter for an element outside a loop. Used to manipulate the field directly outside loops.
      inline auto get_element(const int i) const {
        return payload.get_element( i, lattice );
      }

      template<typename A>
      inline void set_element(const A & value, const int i) {
        payload.set_element( value, i, lattice );
      }
#else
      template <typename vecT>
      inline vecT get_vector(const int i) const {
        return payload.template get_vector<vecT>( i );
      }
      inline T get_element(const int i) const {
        return payload.get_element( i );
      }

      template <typename vecT>
      inline void set_vector(const vecT & val, const int i) {
        return payload.set_vector( val, i );
      }
      inline void set_element(const T & val, const int i) {
        return payload.set_element( val, i );
      }
#endif

      /// Gather boundary elements for communication
      void gather_comm_elements(direction d, parity par, T * RESTRICT buffer, 
                                const lattice_struct::comm_node_struct & to_node) const {
#ifndef VECTORIZED
#ifdef SPECIAL_BOUNDARY_CONDITIONS
        // note: -d in is_on_edge, because we're about to send stuff to that direction
        // (fetching from direction +d)
        if (boundary_condition[d] == boundary_condition_t::ANTIPERIODIC &&
            lattice->special_boundaries[-d].is_on_edge) {
          payload.gather_comm_elements(buffer, to_node, par, lattice, true);
        } else {
          payload.gather_comm_elements(buffer, to_node, par, lattice, false);
        }
#else
        payload.gather_comm_elements(buffer, to_node, par, lattice, false);
#endif

#else
        // this is vectorized branch
        bool antiperiodic = false;
        if (boundary_condition[d] == boundary_condition_t::ANTIPERIODIC &&
            lattice->special_boundaries[-d].is_on_edge) antiperiodic = true;

        if constexpr (is_vectorizable_type<T>::value) {
          // now vectorized layout
          if (vector_lattice->is_boundary_permutation[abs(d)]) {
            // with boundary permutation need to fetch elems 1-by-1
            int n;
            const unsigned * index_list = to_node.get_sitelist(par,n);
            if (!antiperiodic)
              payload.gather_elements(buffer, index_list, n, lattice);
            else {
              payload.gather_elements_negated(buffer, index_list, n, lattice);
            }
          } else {
            // without it, can do the full block
            payload.gather_comm_vectors(buffer, to_node, par, vector_lattice, antiperiodic);
          }
        } else {
          // not vectoizable, standard methods
          int n;
          const unsigned * index_list = to_node.get_sitelist(par,n);
          if (!antiperiodic)
            payload.gather_elements(buffer, index_list, n, lattice);
          else {
            payload.gather_elements_negated(buffer, index_list, n, lattice);
          }
        }
#endif
      }


      /// Place boundary elements from neighbour
      void place_comm_elements(direction d, parity par, T * RESTRICT buffer, 
                               const lattice_struct::comm_node_struct & from_node){
// #ifdef USE_MPI
#ifdef VECTORIZED
        if constexpr (is_vectorizable_type<T>::value) {
          // now vectorized layout, act accordingly
          if (vector_lattice->is_boundary_permutation[abs(d)]) {
            payload.place_recv_elements(buffer, d, par, vector_lattice);
          } else {
            // nothing to do here, comms directly in place
          }
        } else {
          // non-vectorized, using vanilla method, again nothing to do
        }
#else
        // this one is only for CUDA
        payload.place_comm_elements(d, par, buffer, from_node, lattice);
#endif
// #endif
      }

      /// Place boundary elements from local lattice (used in vectorized version)
      void set_local_boundary_elements(direction dir, parity par){
        bool antiperiodic =
          (boundary_condition[dir] == boundary_condition_t::ANTIPERIODIC && lattice->special_boundaries[dir].is_on_edge);
        payload.set_local_boundary_elements(dir, par, lattice, antiperiodic);
      }

      /// Gather a list of elements to a single node
      void gather_elements(T * buffer, std::vector<coordinate_vector> coord_list, int root=0) const;
      void send_elements(T * buffer, std::vector<coordinate_vector> coord_list, int  root=0);

#if defined(USE_MPI)

      /// get the receive buffer pointer for the communication.    
      T * get_receive_buffer( direction d, parity par,
                              const lattice_struct::comm_node_struct & from_node ) {
#if defined(VANILLA)

        return (T *)payload.get_buffer() + from_node.offset(par);

#elif defined(CUDA)

        int offs = 0;
        if (par == ODD) offs = from_node.sites/2;
        if (receive_buffer[d] == nullptr) {
          receive_buffer[d] = payload.allocate_mpi_buffer( from_node.sites );
        }
        return receive_buffer[d] + offs;

#elif defined(VECTORIZED)

        if constexpr (!is_vectorizable_type<T>::value) {
          // use vanilla type, field laid out in std fashion
          return (T *)payload.get_buffer() + from_node.offset(par);
        } else {
          int offs = 0;
          if (par == ODD) offs = from_node.sites/2;

          if (vector_lattice->is_boundary_permutation[abs(d)]) {
            // extra copy operation needed
            if (receive_buffer[d] == nullptr) { 
              receive_buffer[d] = payload.allocate_mpi_buffer( from_node.sites );
            }
            return receive_buffer[d] + offs;
          } else {
            // directly to halo buffer
            constexpr int vector_size = vector_info<T>::vector_size;
            return ((T *)payload.get_buffer()
                    + (vector_lattice->halo_offset[d]*vector_size + offs) );
          }
        }
#endif
      } // end of get_receive_buffer
#endif  // USE_MPI

  };

  // static_assert( std::is_pod<T>::value, "Field expects only pod-type elements (plain data): default constructor, copy and delete");
  static_assert( std::is_trivial<T>::value && std::is_standard_layout<T>::value,
                "Field expects only pod-type elements (plain data): default constructor, copy and delete");
  
 public:

  field_struct * RESTRICT fs;
  
  field() {
    // std::cout << "In constructor 1\n";
    fs = nullptr;             // lazy allocation on 1st use
  }
  
  // Straightforward copy constructor seems to be necessary
  field(const field & other) {
    fs = nullptr;  // this is probably unnecessary
    if(other.fs != nullptr){
      (*this)[ALL] = other[X];
    }
  }
    
  // copy constructor - from fields which can be assigned
  template <typename A,
            std::enable_if_t<std::is_convertible<A,T>::value, int> = 0 >  
  field(const field<A>& other) {
    fs = nullptr;  // this is probably unnecessary
    if(other.fs != nullptr){
      (*this)[ALL] = other[X];
    }
  }

  // constructor with compatible scalar
  template <typename A,
            std::enable_if_t<std::is_convertible<A,T>::value, int> = 0 >  
  field(const A& val) {
    fs = nullptr;
    // static_assert(!std::is_same<A,int>::value, "in int constructor");
    (*this)[ALL] = val;
  }
  
  // move constructor - steal the content
  field(field && rhs) {
    // std::cout << "in move constructor\n";
    fs = rhs.fs;
    rhs.fs = nullptr;
  }

  ~field() {
    free();
  }
    
  void allocate() {
    assert(fs == nullptr);
    if (lattice == nullptr) {
      output0 << "Can not allocate field variables before lattice.setup()\n";
      hila::terminate(0); 
    }
    fs = new field_struct;
    fs->lattice = lattice;
    fs->allocate_payload();
    fs->initialize_communication();
    mark_changed(ALL);      // guarantees communications will be done
    fs->assigned_to = 0;    // and this means that it is not assigned

    for (direction d=(direction)0; d<NDIRS; ++d){
#ifndef CUDA
      fs->neighbours[d] = lattice->neighb[d];
#else
      fs->payload.neighbours[d] = lattice->backend_lattice->d_neighb[d];
#endif
    }
#ifdef SPECIAL_BOUNDARY_CONDITIONS
    foralldir(dir){
      fs->boundary_condition[dir] = boundary_condition_t::PERIODIC;
      fs->boundary_condition[opp_dir(dir)] = boundary_condition_t::PERIODIC;
    }
#endif

    #ifdef VECTORIZED
    fs->vector_lattice = lattice->backend_lattice->get_vectorized_lattice< vector_info<T>::vector_size >();
    #endif
  }

  void free() {
    // don't call destructors when exiting - either MPI or cuda can already
    // be off.  
    if (fs != nullptr && !hila::about_to_finish) {
      for (direction d=(direction)0; d<NDIRS; ++d) drop_comms(d,ALL);
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

    
  fetch_status move_status(parity p, int d) const { 
    assert(parity_bits(p) && d>=0 && d<NDIRS);
    return fs->move_status[(int)p - 1][d]; 
  }
  void set_move_status(parity p, int d, fetch_status stat) const { 
    assert(parity_bits(p) && d>=0 && d<NDIRS);
    fs->move_status[(int)p - 1][d] = stat;
  }

  // check that field is allocated, and if not do it (if not const)
  // call this BEFORE the var is actually written to
  void check_alloc() { 
    if (!is_allocated()) allocate();
  }
  // If field is const specified, we should not be able to write to it in the first
  // place
  void check_alloc() const { 
    assert(is_allocated());
  }

  // If ALL changes, both parities invalid; if p != ALL, then p and ALL.
  void mark_changed(const parity p) const {

    for (direction i=(direction)0; i<NDIRS; ++i) {
      // check if there's ongoing comms, invalidate it!
      drop_comms(i,opp_parity(p));

      set_move_status(opp_parity(p),i,fetch_status::NOT_DONE);
      if (p != ALL){
        set_move_status(ALL,i,fetch_status::NOT_DONE );
      } else {
        set_move_status(EVEN,i,fetch_status::NOT_DONE);
        set_move_status(ODD,i,fetch_status::NOT_DONE);
      }
    }
    fs->assigned_to |= parity_bits(p);
  }

  
  /// Mark the field parity fetched from direction
  // In case p=ALL we could mark everything fetched, but we'll be conservative here 
  // and mark only this parity, because there might be other parities on the fly and corresponding
  // waits should be done,  This should never happen in automatically generated loops.
  // In any case start_fetch, is_fetched, get_move_parity has intelligence to figure out the right thing to do
  //

  void mark_fetched( int dir, const parity p) const {
    set_move_status(p,dir,fetch_status::DONE);
  }

  // Check if the field has been fetched since the previous communication
  // par = ALL:   ALL or (EVEN+ODD) are OK
  // par != ALL:  ALL or par are OK
  bool is_fetched( int dir, parity par) const {
    if (par != ALL) {
      return move_status(par,dir) == fetch_status::DONE || move_status(ALL,dir) == fetch_status::DONE;
    } else {
      return move_status(ALL,dir) == fetch_status::DONE ||
           ( move_status(EVEN,dir) == fetch_status::DONE && move_status(ODD,dir) == fetch_status::DONE );
    }
  }
   
  // Mark communication started -- this must be just the one
  // going on with MPI
  void mark_move_started( int dir, parity p) const{
    set_move_status(p,dir,fetch_status::STARTED);
  }

  /// Check if communication has started.  This is strict, checks exactly this parity
  bool is_move_started( int dir, parity par) const {
    return move_status(par,dir) == fetch_status::STARTED;
  }
    
  bool move_not_done( int dir, parity par) const {
    return move_status(par,dir) == fetch_status::NOT_DONE;
  }
 
  void set_boundary_condition( direction dir, boundary_condition_t bc) {

    #ifdef SPECIAL_BOUNDARY_CONDITIONS
    // TODO: This works as intended only for periodic/antiperiodic b.c.
    check_alloc();
    fs->boundary_condition[dir] = bc;
    fs->boundary_condition[-dir] = bc;
    #ifndef CUDA
    fs->neighbours[dir]  = lattice->get_neighbour_array(dir,bc);
    fs->neighbours[-dir] = lattice->get_neighbour_array(-dir,bc);
    #else
    if(bc == boundary_condition_t::PERIODIC){
      fs->payload.neighbours[dir] = lattice->backend_lattice->d_neighb[dir];
      fs->payload.neighbours[-dir] = lattice->backend_lattice->d_neighb[-dir];
    } else {
      fs->payload.neighbours[dir] = lattice->backend_lattice->d_neighb_special[dir];
      fs->payload.neighbours[-dir] = lattice->backend_lattice->d_neighb_special[-dir];
    }
    #endif

    // Make sure boundaries get refreshed
    mark_changed(ALL);
    #endif
  }

  boundary_condition_t get_boundary_condition( direction dir ) const {
    #ifdef SPECIAL_BOUNDARY_CONDITIONS
    return fs->boundary_condition[dir];
    #else
    return boundary_condition_t::PERIODIC;
    #endif
  }

  void print_boundary_condition() {
    check_alloc();
    output0 << " ( ";
    for(int dir=0; dir<NDIRS; dir++){
      output0 << (int)fs->boundary_condition[dir] << " ";
    }
    output0 << ")\n";
  }

  template<typename A>
  void copy_boundary_condition(const field<A> & rhs) {
    foralldir(dir){
      set_boundary_condition(dir, rhs.get_boundary_condition(dir));
    }
  }

  // Overloading [] 
  // placemarker, should not be here
  // T& operator[] (const int i) { return data[i]; }

  // declarations -- WILL BE implemented by hilapp, not written here
  element<T>& operator[] (const parity p) const;             // f[EVEN]
  element<T>& operator[] (const X_index_type) const;         // f[X]
  element<T>& operator[] (const X_plus_direction p) const;   // f[X+dir]
  element<T>& operator[] (const X_plus_offset p) const;      // f[X+dir1+dir2] and others


  // TEMPORARY HACK: return ptr to bare array
  inline auto field_buffer() const { return fs->payload.get_buffer(); }


#ifndef VECTORIZED
  /// Get an individual element outside a loop. This is also used as a getter in the vanilla code.
  inline auto get_value_at(int i) const { return fs->get_element(i); }
#else
  inline auto get_value_at(int i) const { return fs->get_element(i); }
  template <typename vecT>
  inline auto get_vector_at(int i) const { return fs->template get_vector<vecT>(i); }
  inline auto get_value_at_nb_site(direction d, int i) const {
    return fs->get_element( fs->vector_lattice->site_neighbour(d,i) );
  }
#endif

#ifndef VECTORIZED
  /// Set an individual element outside a loop. This is also used as a setter in the vanilla code.
  template<typename A>
  inline void set_value_at(const A & value, int i) { fs->set_element( value, i); }

#else
  template<typename vecT>
  inline void set_vector_at(const vecT & value, int i) { fs->set_vector( value, i); }

  template<typename A>
  inline void set_value_at(const A & value, int i) { fs->set_element( value, i); }
#endif

  // fetch the element at this loc
  // T get(int i) const;
  
  // Basic copy constructor (cannot be a template)
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
  dir_mask_t start_fetch(direction d, parity p) const;
  dir_mask_t start_fetch(direction d) const { return start_fetch(d, ALL);}
  void wait_fetch(direction d, parity p) const;
  void fetch(direction d, parity p) const;
  void drop_comms(direction d, parity p) const;
  void cancel_comm(direction d, parity p) const; 

  // Declaration of shift methods
  field<T> shift(const coordinate_vector &v, parity par) const;
  field<T> shift(const coordinate_vector &v) const { return shift(v,ALL); }

  // General getters and setters
  void set_elements(T * elements, std::vector<coordinate_vector> coord_list);
  void set_element(T element, coordinate_vector coord);
  void get_elements(T * elements, std::vector<coordinate_vector> coord_list) const;
  T get_element(coordinate_vector coord) const;

  // Fourier transform declarations
  void FFT(fft_direction fdir = fft_direction::forward);

  // Writes the field to disk
  void write_to_stream(std::ofstream & outputfile);
  void write_to_file(std::string filename);
  void read_from_stream(std::ifstream & inputfile);
  void read_from_file(std::string filename);
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
 * For simplicity, these functions do not use field expressions and
 * can be ignored by the hilapp. Since the hilapp does not
 * have access to mpi.h, it cannot process this branch.
 */

/// start_fetch(): Communicate the field at parity par from direction
/// d. Uses accessors to prevent dependency on the layout.
/// return the direction mask bits where something is happening
template<typename T>
dir_mask_t field<T>::start_fetch(direction d, parity p) const {

  // get the mpi message tag right away, to ensure that we are always synchronized with the
  // mpi calls -- some nodes might not need comms, but the tags must be in sync

  int tag = get_next_msg_tag();


  lattice_struct::nn_comminfo_struct  & ci = lattice->nn_comminfo[d];
  lattice_struct::comm_node_struct & from_node = ci.from_node;
  lattice_struct::comm_node_struct & to_node = ci.to_node;

  // check if this is done - either fetched or no comm to be done in the 1st place

  if (is_fetched(d,p)) {
    lattice->n_gather_avoided++; 
    return 0;   // nothing to wait for
  }

  // No comms to do, nothing to wait for -- we'll use the is_fetched
  // status to keep track of vector boundary shuffle anyway

  if (from_node.rank == hila::myrank() && to_node.rank == hila::myrank()) {
    fs->set_local_boundary_elements(d,p);
    mark_fetched(d,p);
    return 0;
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

  T * receive_buffer;
  T * send_buffer;



  if (from_node.rank != hila::myrank()) {

    // HANDLE RECEIVES: get node which will send here
    post_receive_timer.start();

    // buffer can be separate or in field buffer
    receive_buffer = fs->get_receive_buffer(d,par,from_node);
  
    unsigned sites = from_node.n_sites(par);

    // c++ version does not return errors??
    MPI_Irecv( receive_buffer, sites*size, MPI_BYTE, from_node.rank,
	             tag, lattice->mpi_comm_lat, &fs->receive_request[par_i][d] );
    
    post_receive_timer.stop();
  }

  if (to_node.rank != hila::myrank()) {
    // HANDLE SENDS: Copy field elements on the boundary to a send buffer and send
    start_send_timer.start();

    unsigned sites = to_node.n_sites(par);

    if(fs->send_buffer[d] == nullptr)
      fs->send_buffer[d] = fs->payload.allocate_mpi_buffer(to_node.sites);
    send_buffer = fs->send_buffer[d] + to_node.offset(par);

    fs->gather_comm_elements(d,par,send_buffer,to_node);
 
    MPI_Isend( send_buffer, sites*size, MPI_BYTE, to_node.rank, 
               tag, lattice->mpi_comm_lat, &fs->send_request[par_i][d]);
    start_send_timer.stop();
  }

  // and do the boundary shuffle here, after MPI has started
  // NOTE: there should be no danger of MPI and shuffle overwriting, MPI writes
  // to halo buffers only if no permutation is needed.  With a permutation MPI
  // uses special receive buffer
  fs->set_local_boundary_elements(d,par);

  return get_dir_mask(d);

}

///  wait_fetch(): Wait for communication at parity par from
///  direction d completes the communication in the function.
///  If the communication has not started yet, also calls
///  start_fetch()
///
///  NOTE: This will be called even if the field is marked const.
///  Therefore this function is const, even though it does change
///  the internal content of the field, the halo. From the point
///  of view of the user, the value of the field does not change.
template<typename T>
void field<T>::wait_fetch(direction d, parity p) const {

  lattice_struct::nn_comminfo_struct  & ci = lattice->nn_comminfo[d];
  lattice_struct::comm_node_struct & from_node = ci.from_node;
  lattice_struct::comm_node_struct & to_node = ci.to_node;

  // check if this is done - either fetched or no comm to be done in the 1st place
  if (is_fetched(d,p)) return;

  // this is the branch if no comms -- shuffle was done in start_fetch
  if (from_node.rank == hila::myrank() && to_node.rank == hila::myrank()) return;

  // if (!is_move_started(d,p)) {
  //   output0 << "Wait move error - wait_fetch without corresponding start_fetch\n";
  //   exit(-1);
  // }

  // Note: the move can be parity p OR ALL -- need to wait for it in any case
  // set par to be the "sum" over both parities
  // There never should be ongoing ALL and other parity fetch -- start_fetch takes care

  // check here consistency, this should never happen
  if (p != ALL && is_move_started(d,p) && is_move_started(d,ALL)) {
    exit(-1);
  }

  parity par;
  int n_wait = 1;
  // what par to wait for?
  if (is_move_started(d,p)) par = p;            // standard match
  else if (p != ALL) {
    if (is_move_started(d,ALL)) par = ALL;      // if all is running wait for it
    else {
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
      exit(-1);
    }
  }

  if (n_wait == 2) par = EVEN; // we'll flip both

  for (int wait_i = 0; wait_i < n_wait; ++wait_i ) {

    int par_i = (int)par - 1;

    if (from_node.rank != hila::myrank()) {
      wait_receive_timer.start();

      MPI_Status status;
      MPI_Wait( &fs->receive_request[par_i][d], &status);

      wait_receive_timer.stop();

#ifndef VANILLA
      fs->place_comm_elements(d, par, fs->get_receive_buffer(d,par,from_node), from_node);
#endif
    }

    // then wait for the sends
    if (to_node.rank != hila::myrank()) {
      wait_send_timer.start();
      MPI_Status status;
      MPI_Wait( &fs->send_request[par_i][d], &status );
      wait_send_timer.stop();
    }

    // Mark the parity fetched from direction dir
    mark_fetched(d, par);
  
    // Keep count of communications
    lattice->n_gather_done += 1;

    par = opp_parity(par);  // flip if 2 loops
  }
}


///  drop_comms():  if field is changed or deleted,
///  cancel ongoing communications.  This should happen very seldom,
///  only if there are "by-hand" start_fetch operations and these are not needed
template<typename T>
void field<T>::drop_comms(direction d, parity p) const {

  if (is_comm_initialized()) {
    if (is_move_started(d, ALL)) cancel_comm(d, ALL);
    if (p != ALL) {
      if (is_move_started(d, p)) cancel_comm(d, p);
    } else {
      if (is_move_started(d, EVEN)) cancel_comm(d, EVEN);
      if (is_move_started(d, ODD)) cancel_comm(d, ODD);
    }
  }
}

/// cancel ongoing send and receive

template<typename T>
void field<T>::cancel_comm(direction d, parity p) const {
  if (lattice->nn_comminfo[d].from_node.rank != hila::myrank()) {
    cancel_receive_timer.start();
    MPI_Cancel( &fs->receive_request[(int)p-1][d] );
    cancel_receive_timer.stop();
  }
  if (lattice->nn_comminfo[d].to_node.rank != hila::myrank()) {
    cancel_send_timer.start();
    MPI_Cancel( &fs->send_request[(int)p-1][d] );
    cancel_send_timer.stop();
  }
}

#else  // No MPI now

///* Trivial implementation when no MPI is used

template<typename T>
dir_mask_t field<T>::start_fetch(direction d, parity p) const {
  // Update local elements in the halo (necessary for vectorized version)
  // We use here simpler tracking than in MPI, may lead to slight extra work
  if (!is_fetched(d,p)) {
    fs->set_local_boundary_elements(d, p);
    mark_fetched(d,p);
  }
  return 0;
}

template<typename T>
void field<T>::wait_fetch(direction d, parity p) const {}

template<typename T>
void field<T>::drop_comms(direction d, parity p) const {}

#endif  // MPI

/// And a convenience combi function
template<typename T>
void field<T>::fetch(direction d, parity p) const {
  start_fetch(d,p);
  wait_fetch(d,p);
}







#if defined(USE_MPI)

/// Gather a list of elements to a single node
template<typename T>
void field<T>::field_struct::gather_elements(T * buffer, std::vector<coordinate_vector> coord_list, int root) const {
  std::vector<unsigned> index_list;
  std::vector<unsigned> node_list(lattice->n_nodes());
  std::fill(node_list.begin(), node_list.end(),0);
  
  for(coordinate_vector c : coord_list){
    if( lattice->is_on_node(c) ){
      index_list.push_back(lattice->site_index(c));
    }

    node_list[lattice->node_rank(c)]++;
  }
  
  std::vector<T> send_buffer(index_list.size());
  payload.gather_elements((T*) send_buffer.data(), index_list.data(), send_buffer.size(), lattice);
  if(hila::myrank() != root && node_list[hila::myrank()] > 0){
    MPI_Send((char*) send_buffer.data(), node_list[hila::myrank()]*sizeof(T), MPI_BYTE, root, hila::myrank(), 
             lattice->mpi_comm_lat);
  }
  if(hila::myrank() == root) {
    for( int n=0; n<node_list.size(); n++ ) if(node_list[n] > 0) {
      if(n!=root) {
        MPI_Status status;
        MPI_Recv(buffer, node_list[n]*sizeof(T), MPI_BYTE, n, n, lattice->mpi_comm_lat, &status);
      } else {
        std::memcpy( buffer, (char *) send_buffer.data(), node_list[n]*sizeof(T) );
      }
      buffer += node_list[n];
    }
  }
}


/// Send elements from a single node to a list of coordinates
template<typename T>
void field<T>::field_struct::send_elements(T * buffer, std::vector<coordinate_vector> coord_list, int root) {
  std::vector<unsigned> index_list;
  std::vector<unsigned> node_list(lattice->n_nodes());
  std::fill(node_list.begin(), node_list.end(),0);

  for(coordinate_vector c : coord_list){
    if( lattice->is_on_node(c) ){
      index_list.push_back(lattice->site_index(c));
    }

    node_list[lattice->node_rank(c)]++;
  }

  std::vector<T> recv_buffer(index_list.size());
  payload.gather_elements((T*) recv_buffer.data(), index_list.data(), recv_buffer.size(), lattice);
  if(hila::myrank() != root && node_list[hila::myrank()] > 0){
    MPI_Status status;
    MPI_Recv((char*) recv_buffer.data(), node_list[hila::myrank()]*sizeof(T), MPI_BYTE, root, hila::myrank(), 
              lattice->mpi_comm_lat, &status);
  }
  if(hila::myrank() == root) {
    for( int n=0; n<node_list.size(); n++ ) if(node_list[n] > 0) {
      if(n!=root) {
        MPI_Send(buffer, node_list[n]*sizeof(T), MPI_BYTE, n, n, lattice->mpi_comm_lat);
      } else {
        std::memcpy( (char *) recv_buffer.data(), buffer, node_list[n]*sizeof(T) );
      }
      buffer += node_list[n];
    }
  }
  payload.place_elements((T*) recv_buffer.data(), index_list.data(), recv_buffer.size(), lattice);
}


#else  // Now not USE_MPI

/// Gather a list of elements to a single node
template<typename T>
void field<T>::field_struct::gather_elements(T * buffer, std::vector<coordinate_vector> coord_list, int root) const {
  std::vector<unsigned> index_list;
  for(coordinate_vector c : coord_list){
    index_list.push_back(lattice->site_index(c));
  }
  
  payload.gather_elements(buffer, index_list.data(), index_list.size(), lattice);
}


/// Send elements from a single node to a list of coordinates
template<typename T>
void field<T>::field_struct::send_elements(T * buffer, std::vector<coordinate_vector> coord_list, int root) {
  std::vector<unsigned> index_list;
  for(coordinate_vector c : coord_list){
    index_list.push_back(lattice->site_index(c));
  }
  
  payload.place_elements(buffer, index_list.data(), index_list.size(), lattice);
}

#endif




/// Functions for manipulating individual elements in an array

/// Set an element. Assuming that each node calls this with the same value, it is
/// sufficient to set the elements locally
template<typename T>
void field<T>::set_elements( T * elements, std::vector<coordinate_vector> coord_list) {
  std::vector<unsigned> my_indexes;
  std::vector<unsigned> my_elements;
  for(int i=0; i<coord_list.size(); i++){
    coordinate_vector c = coord_list[i];
    if( lattice->is_on_node(c) ){
      my_indexes.push_back(lattice->site_index(c));
      my_elements.push_back(elements[i]);
    }
  }
  fs->payload.place_elements(my_elements.data(), my_indexes.data(), my_indexes.size(), lattice);
  mark_changed(ALL);
}

// Set a single element. Assuming that each node calls this with the same value, it is
/// sufficient to set the element locally
template<typename T>
void field<T>::set_element( T element, coordinate_vector coord) {
  if( lattice->is_on_node(coord) ){
    set_value_at( element, lattice->site_index(coord));
  }
  mark_changed(ALL);
}



/// Get an element and return it on all nodes
#if defined(USE_MPI)
/// This is not local, the element needs to be communicated to all nodes
template<typename T>
T field<T>::get_element( coordinate_vector coord) const {
  T element;
  int owner = lattice->node_rank(coord);

  if( hila::myrank() == owner ){
    element = get_value_at( lattice->site_index(coord) );
  }

  MPI_Bcast( &element, sizeof(T), MPI_BYTE, owner, lattice->mpi_comm_lat);
  return element;
}


/// Get a list of elements and store them into an array on all nodes
template<typename T>
void field<T>::get_elements( T * elements, std::vector<coordinate_vector> coord_list) const {
  struct node_site_list_struct {
    std::vector<int> indexes;
    std::vector<coordinate_vector> coords;
  };
  std::vector<node_site_list_struct> nodelist(lattice->n_nodes());
  // Reorganize the list according to nodes
  for( int i=0; i<coord_list.size(); i++){
    coordinate_vector c = coord_list[i];
    int node = lattice->node_rank(c);
    nodelist[node].indexes.push_back(i);
    nodelist[node].coords.push_back(c);
  }

  // Fetch on each node found and communicate
  for(int n=0; n<nodelist.size(); n++){
    node_site_list_struct node = nodelist[n];
    if(node.indexes.size() > 0){
      T * element_buffer = (T *) malloc( sizeof(T)*node.indexes.size() );
      fs->payload.gather_elements( element_buffer, node.coords);

      MPI_Bcast( &element_buffer, sizeof(T), MPI_BYTE, n, lattice->mpi_comm_lat);

      // place in the array in original order
      for( int i=0; i<node.indexes.size(); i++){
        elements[i] = element_buffer[node.indexes[i]];
      }

      free(element_buffer);
    }
  }
}


#else
/// Without MPI, we just need to call get
template<typename T>
T field<T>::get_element( coordinate_vector coord) const {
  return get_value_at( lattice->site_index(coord) );
}

/// Without MPI, we just need to call get
template<typename T>
void field<T>::get_elements( T * elements, std::vector<coordinate_vector> coord_list) const {
  for( int i=0; i<coord_list.size(); i++){
    elements[i] = get_element(coord_list[i]);
  }
}
#endif




// Write the field to an file stream
template<typename T>
void field<T>::write_to_stream(std::ofstream& outputfile){
  constexpr size_t target_write_size = 1000000;
  constexpr size_t sites_per_write = target_write_size / sizeof(T);
  constexpr size_t write_size = sites_per_write * sizeof(T);

  std::vector<coordinate_vector> coord_list(sites_per_write);
  T * buffer = (T*) malloc(write_size);
  coordinate_vector size = lattice->size();

  int i=0;
  for(; i<lattice->volume(); i++){
    coordinate_vector site;
    int ii = i;
    foralldir(dir){
      site[dir] = ii%size[dir];
      ii = ii/size[dir];
    }

    coord_list[i%sites_per_write] = site;

    // Write the buffer when full
    if( (i+1)%sites_per_write == 0 ){
      fs->gather_elements(buffer, coord_list);
      if( hila::myrank()==0 )
        outputfile.write((char*)buffer,write_size);
    }
  }

  // Write the rest
  coord_list.resize(i%sites_per_write);
  fs->gather_elements(buffer, coord_list);
  double * v = (double*) buffer;
  if( hila::myrank() == 0 )
    outputfile.write((char*)buffer,sizeof(T)*(i%sites_per_write));

  std::free(buffer);
}


// Write the field to a file replacing the file
template<typename T>
void field<T>::write_to_file(std::string filename){
  std::ofstream outputfile;
  outputfile.open(filename, std::ios::out | std::ios::trunc | std::ios::binary);
  write_to_stream(outputfile);
  outputfile.close();
}



// Write a list of fields into an output stream
template<typename T>
static void write_fields(std::ofstream& outputfile, field<T>& last){
  last.write_to_stream(outputfile);
}

template<typename T, typename... fieldtypes>
static void write_fields(std::ofstream& outputfile, field<T>& next, fieldtypes&... fields){
  next.write_to_stream(outputfile);
  write_fields(outputfile, fields...);
}

// Write a list of fields to a file
template<typename... fieldtypes>
static void write_fields(std::string filename, fieldtypes&... fields){
  std::ofstream outputfile;
  outputfile.open(filename, std::ios::out | std::ios::trunc | std::ios::binary);
  write_fields(outputfile, fields...);
  outputfile.close();
}



// Read the field from a stream
template<typename T>
void field<T>::read_from_stream(std::ifstream& inputfile){
  constexpr size_t target_read_size = 1000000;
  constexpr size_t sites_per_read = target_read_size / sizeof(T);
  constexpr size_t read_size = sites_per_read * sizeof(T);

  mark_changed(ALL);

  std::vector<coordinate_vector> coord_list(sites_per_read);
  T * buffer = (T*) malloc(read_size);
  coordinate_vector size = lattice->size();

  int i=0;
  for(; i<lattice->volume(); i++){
    coordinate_vector site;
    int ii = i;
    foralldir(dir){
      site[dir] = ii%size[dir];
      ii = ii/size[dir];
    }

    coord_list[i%sites_per_read] = site;

    // Read the buffer when full
    if( (i+1)%sites_per_read == 0 ){
      if( hila::myrank()==0 )
        inputfile.read((char*)buffer,read_size);
      fs->send_elements(buffer, coord_list);
    }
  }

  // Read the rest
  coord_list.resize(i%sites_per_read);
  if( hila::myrank()==0 )
    inputfile.read((char*)buffer, sizeof(T)*(i%sites_per_read));
  double * v = (double*) buffer;
  fs->send_elements(buffer, coord_list);

  std::free(buffer);
}


// Read field contennts from the beginning of a file
template<typename T>
void field<T>::read_from_file(std::string filename){
  std::ifstream inputfile;
  inputfile.open(filename, std::ios::in | std::ios::binary);
  read_from_stream(inputfile);
  inputfile.close();
}


// Read a list of fields from an input stream
template<typename T>
static void read_fields(std::ifstream& inputfile, field<T>& last){
  last.read_from_stream(inputfile);
}

template<typename T, typename... fieldtypes>
static void read_fields(std::ifstream& inputfile, field<T>& next, fieldtypes&... fields){
  next.read_from_stream(inputfile);
  read_fields(inputfile, fields...);
}

// Read a list of fields from a file
template<typename... fieldtypes>
static void read_fields(std::string filename, fieldtypes&... fields){
  std::ifstream inputfile;
  inputfile.open(filename, std::ios::in | std::ios::binary);
  read_fields(inputfile, fields...);
  inputfile.close();
}

//HACK: force disable vectorization in a loop using
// if(disable_avx[X]==0){};
extern field<double> disable_avx;



// Include FFT implentations (mostly fftw, but for CUDA
// we need to specialize version)
#ifdef CUDA
#include "plumbing/backend_cuda/FFT.h"
#else
#include "plumbing/FFT.h"
#endif




#endif // FIELD_H


