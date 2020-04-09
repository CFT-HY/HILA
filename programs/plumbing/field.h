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
      unsigned is_fetched[NDIRS];     // is communication done
      unsigned move_started[NDIRS];   // is communication going on
      unsigned assigned_to;           // keeps track of first assignment to parities
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

      /// Gather a list of elements to a single node
#if defined(USE_MPI) && !defined(TRANSFORMER) 
      void gather_elements(char * buffer, std::vector<unsigned> index_list, std::vector<unsigned> node_list, int root, MPI_Comm Communicator) const;
      void gather_elements(char * buffer, std::vector<coordinate_vector> coord_list, int root=0, MPI_Comm Communicator=MPI_COMM_WORLD) const;
      void send_elements(char * buffer, std::vector<unsigned> index_list, std::vector<unsigned> node_list, int  root=0, MPI_Comm Communicator=MPI_COMM_WORLD);
      void send_elements(char * buffer, std::vector<coordinate_vector> coord_list, int  root=0, MPI_Comm Communicator=MPI_COMM_WORLD);
#else
      void gather_elements(char * buffer, std::vector<unsigned> index_list, int root=0) const;
      void gather_elements(char * buffer, std::vector<coordinate_vector> coord_list, int root=0) const;
      void send_elements(char * buffer, std::vector<unsigned> index_list, int root=0);
      void send_elements(char * buffer, std::vector<coordinate_vector> coord_list, int root=0);
#endif
  };

  static_assert( std::is_pod<T>::value, "Field expects only pod-type elements (plain data): default constructor, copy and delete");
  
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
    mark_changed(ALL);      // guarantees communications will be done
    fs->assigned_to = 0;    // and this means that it is not assigned
  }

  void free() {
    if (fs != nullptr) {
      fs->free_payload();
      delete fs;
      fs = nullptr;
    }
  }

  bool is_allocated() const { return (fs != nullptr); }

  bool is_initialized(parity p) const { 
    return fs != nullptr && ((fs->assigned_to & parity_bits(p)) != 0);
  }
  
  /// call this BEFORE the var is written to
  void mark_changed(const parity p) {
    if (fs == nullptr) allocate();
    else {
      // turn off bits corresponding to parity p
      assert( parity_bits(p) );
      for (int i=0; i<NDIRS; i++) fs->is_fetched[i]   &= parity_bits_inverse(p);
      for (int i=0; i<NDIRS; i++) fs->move_started[i] &= parity_bits_inverse(p);
    }
    fs->assigned_to |= parity_bits(p);
  }

  // Is const version of mark_changed needed?  Sounds strange
  void mark_changed(const parity p) const {
    assert(is_allocated());
    assert( parity_bits(p) );
    for (int i=0; i<NDIRS; i++) fs->is_fetched[i]   &= parity_bits_inverse(p);
    for (int i=0; i<NDIRS; i++) fs->move_started[i] &= parity_bits_inverse(p);
    fs->assigned_to |= parity_bits(p);
  }

  /// Mark the field parity fetched from direction
  void mark_fetched( int dir, const parity p) const {
    assert( parity_bits(p) );
    fs->is_fetched[dir] |= parity_bits(p);
  }

  /// Check if the field has been changed since the previous communication
  bool is_fetched( int dir, parity par) const {
    assert(dir < NDIRS);
    unsigned p = parity_bits(par);
    // true if all par-bits are on 
    return (fs->is_fetched[dir] & p) == p ;
  }

  /* Mark communication started */
  void mark_move_started( int dir, parity p) const{
    assert(dir < NDIRS);
    fs->move_started[dir] |= parity_bits(p);
  }

  /// Check if communication has started
  bool is_move_started( int dir, parity par) const{
    assert(dir < NDIRS);
    unsigned p = parity_bits(par);
    return (fs->move_started[dir] & p) == p ;
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
  void start_move(direction d) const { start_move(d, ALL);}
  void wait_move(direction d, parity p) const;

  // Declaration of shift methods
  field<T> shift(const coordinate_vector &v, parity par) const;
  field<T> shift(const coordinate_vector &v) const { return shift(v,ALL); }

  // General getters and setters
  void set_elements( T * elements, std::vector<coordinate_vector> coord_list) const;
  void set_elements( T element, coordinate_vector coord) const;

  // Fourier transform declarations
  void FFT();

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



/// Functions for manipulating lists of elements
template<typename T>
void field<T>::set_elements( T * elements, std::vector<coordinate_vector> coord_list) const {
  fs->send_elements( (char*) elements, coord_list);
}

template<typename T>
void field<T>::set_elements( T element, coordinate_vector coord) const {
  std::vector<coordinate_vector> coord_list;
  coord_list.push_back(coord);
  fs->send_elements( (char*) &element, coord_list);
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
    int index = static_cast<int>(d) + NDIRS*static_cast<int>(par);
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

      //printf("node %d, recv tag %d from %d\n", mynode(), tag, from_node.rank);

      MPI_Irecv( receive_buffer[n], sites*size, MPI_BYTE, from_node.rank, 
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
 
       //printf("node %d, send tag %d to %d\n", mynode(), tag, to_node.rank);

       /* And send */
       MPI_Isend( send_buffer[n], sites*size, MPI_BYTE, to_node.rank, 
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
    int index = static_cast<int>(d) + NDIRS*static_cast<int>(par);
    int tag =  fs->mpi_tag*3*NDIRS + index;

    if( is_fetched(d, par) ){
      // Not changed, return directly
      // Keep count of gathers optimized away
      lattice->n_gather_avoided += 1;
      return;
    }

    //printf("wait_move tag %d node %d\n",tag,mynode());

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





/// Gather a list of elements to a single node
#if defined(USE_MPI) && !defined(TRANSFORMER)

template<typename T>
void field<T>::field_struct::gather_elements(char * buffer, std::vector<unsigned> index_list, std::vector<unsigned> node_list, int root, MPI_Comm Communicator) const {
  std::vector<T> send_buffer(index_list.size());
  payload.gather_elements((char*) send_buffer.data(), index_list, lattice);
  if(mynode() != root && node_list[mynode()] > 0){
    MPI_Send((char*) send_buffer.data(), node_list[mynode()]*sizeof(T), MPI_BYTE, root, mynode(), MPI_COMM_WORLD);
  }
  if(mynode() == root) {
    for( int n=0; n<node_list.size(); n++ ) if(node_list[n] > 0) {
      if(n!=root) {
        MPI_Status status;
        MPI_Recv(buffer, node_list[n]*sizeof(T), MPI_BYTE, n, n, MPI_COMM_WORLD, &status);
      } else {
        std::memcpy( buffer, (char *) send_buffer.data(), node_list[n]*sizeof(T) );
      }
      buffer += node_list[n]*sizeof(T);
    }
  }
}

template<typename T>
void field<T>::field_struct::gather_elements(char * buffer, std::vector<coordinate_vector> coord_list, int root, MPI_Comm Communicator) const {
  std::vector<unsigned> index_list;
  std::vector<unsigned> node_list(lattice->n_nodes());
  std::fill(node_list.begin(), node_list.end(),0);
  
  for(coordinate_vector c : coord_list){
    if( lattice->is_on_node(c) ){
      index_list.push_back(lattice->site_index(c));
    }

    node_list[lattice->node_rank(c)]++;
  }
  
  gather_elements(buffer, index_list, node_list, root, Communicator);
}



template<typename T>
void field<T>::field_struct::send_elements(char * buffer, std::vector<unsigned> index_list, std::vector<unsigned> node_list, int root, MPI_Comm Communicator) {
  std::vector<T> recv_buffer(index_list.size());
  payload.gather_elements((char*) recv_buffer.data(), index_list, lattice);
  if(mynode() != root && node_list[mynode()] > 0){
    MPI_Status status;
    MPI_Recv((char*) recv_buffer.data(), node_list[mynode()]*sizeof(T), MPI_BYTE, root, mynode(), MPI_COMM_WORLD, &status);
  }
  if(mynode() == root) {
    for( int n=0; n<node_list.size(); n++ ) if(node_list[n] > 0) {
      if(n!=root) {
        MPI_Send(buffer, node_list[n]*sizeof(T), MPI_BYTE, n, n, MPI_COMM_WORLD);
      } else {
        std::memcpy( (char *) recv_buffer.data(), buffer, node_list[n]*sizeof(T) );
      }
      buffer += node_list[n]*sizeof(T);
    }
  }
  payload.place_elements((char*) recv_buffer.data(), index_list, lattice);
}

template<typename T>
void field<T>::field_struct::send_elements(char * buffer, std::vector<coordinate_vector> coord_list, int root, MPI_Comm Communicator) {
  std::vector<unsigned> index_list;
  std::vector<unsigned> node_list(lattice->n_nodes());
  std::fill(node_list.begin(), node_list.end(),0);

  for(coordinate_vector c : coord_list){
    if( lattice->is_on_node(c) ){
      index_list.push_back(lattice->site_index(c));
    }

    node_list[lattice->node_rank(c)]++;
  }
  
  send_elements(buffer, index_list, node_list, root, Communicator);
}


#else


template<typename T>
void field<T>::field_struct::gather_elements(char * buffer, std::vector<unsigned> index_list, int root) const {
  payload.gather_elements(buffer, index_list, lattice);
}

template<typename T>
void field<T>::field_struct::gather_elements(char * buffer, std::vector<coordinate_vector> coord_list, int root) const {
  std::vector<unsigned> index_list;
  for(coordinate_vector c : coord_list){
    index_list.push_back(lattice->site_index(c));
  }
  
  gather_elements(buffer, index_list);
}


template<typename T>
void field<T>::field_struct::send_elements(char * buffer, std::vector<unsigned> index_list, int root) {
  payload.place_elements(buffer, index_list, lattice);
}

template<typename T>
void field<T>::field_struct::send_elements(char * buffer, std::vector<coordinate_vector> coord_list, int root) {
  std::vector<unsigned> index_list;
  for(coordinate_vector c : coord_list){
    index_list.push_back(lattice->site_index(c));
  }
  
  send_elements(buffer, index_list, root);
}

#endif



// Write the field to an file stream
template<typename T>
void field<T>::write_to_stream(std::ofstream& outputfile){
  constexpr size_t target_write_size = 1000000;
  constexpr size_t sites_per_write = target_write_size / sizeof(T);
  constexpr size_t write_size = sites_per_write * sizeof(T);

  std::vector<coordinate_vector> coord_list(sites_per_write);
  char * buffer = (char*) malloc(write_size);
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
      if( mynode()==0 )
        outputfile.write(buffer,write_size);
    }
  }

  // Write the rest
  coord_list.resize(i%sites_per_write);
  fs->gather_elements(buffer, coord_list);
  double * v = (double*) buffer;
  if( mynode() == 0 )
    outputfile.write(buffer,sizeof(T)*(i%sites_per_write));

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
  char * buffer = (char*) malloc(read_size);
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
      if( mynode()==0 )
        inputfile.read(buffer,read_size);
      fs->send_elements(buffer, coord_list);
    }
  }

  // Read the rest
  coord_list.resize(i%sites_per_read);
  if( mynode()==0 )
    inputfile.read(buffer, sizeof(T)*(i%sites_per_read));
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


// Include Fourier transform
#include "../plumbing/FFT.h"


#endif // FIELD_H


