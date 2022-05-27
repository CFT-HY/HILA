#ifndef FIELD_H
#define FIELD_H
#include <sstream>
#include <iostream>
#include <string>
#include <cstring> //Memcpy is here...
#include <math.h>
#include <type_traits>

#include "plumbing/defs.h"
#include "plumbing/coordinates.h"
#include "plumbing/lattice.h"
#include "plumbing/field_storage.h"

#include "plumbing/backend_vector/vector_types.h"

#ifdef USE_MPI
#include "plumbing/com_mpi.h"
#endif

// This is a marker for hilapp -- will be removed by it
#define onsites(p) for (Parity par_dummy__(p); par_dummy__ == EVEN; par_dummy__ = ODD)

template <typename T>
class Field;

template <typename T>
void ensure_field_operators_exist(Field<T> &f);


/// Field class
/// This implements the standard methods for accessing fields
/// Hilapp replaces the parity access patterns, Field[par] with a loop over
/// the appropriate sites.
///
/// The Field class also contains member functions used by hilapp, as well
/// as members that may be useful for application developers.
///
/// The Field mainly implements the interface to the Field and not the
/// content.
///
/// The Field contains a pointer to Field::field_struct, which implements
/// MPI communication of the Field boundaries.
///
/// The Field::field_struct points to a field_storage, which is defined
/// by each backend. It implements storing and accessing the Field data,
/// including buffers for storing haloes returned from MPI communication.
///
/// Memory allocation (mainly automised by hilapp):
/// Field.allocate(): sets up memory for field content and communication.
/// Field.free(): destroys the data.
/// Field.is_allocated(): returns true if the Field data has been allocated
/// Field.is_initialized() returns true if the Field has been written
/// Field.check_alloc(): allocate if necessary
/// Field.check_alloc() const: assert that the Field is allocated
///
/// MPI related (automatically done by hilapp, but may be useful in apps):
/// Field.move_status(): returns current gather_status
/// Field.mark_changed(): make sure the Field gets communicated
/// Field.mark_gathered(): mark the Field already gathered, no need to
///        communicate.
///
/// Field.shift(): create a periodically shifted copy of the field
///
/// Others
/// Field.set_boundary_condition(): set the boundary conditions in a
///         given Direction (periodic or antiperiodic)
/// Field.get_boundary_condition(): get the boundary condition of the Field
/// Field.copy_boundary_condition(): copy the boundary condition to the
///        from another Field
/// Field.get_elements(): retrieve a list of elements to all nodes
/// Field.get_element(): retrieve an element to all nodes
/// Field.set_elements(): set elements in the Field
/// Field.set_element(): set an element in the Field
///
template <typename T>
class Field {

  public:
    enum class gather_status : unsigned { NOT_DONE, STARTED, DONE };

  private:
    /// The following struct holds the data + information about the field
    /// TODO: field-specific boundary conditions?
    class field_struct {
      public:
        field_storage<T>
            payload; // TODO: must be maximally aligned, modifiers - never null
        lattice_struct *lattice;
#ifdef VECTORIZED
        // get a direct ptr from here too, ease access
        vectorized_lattice_struct<hila::vector_info<T>::vector_size> *vector_lattice;
#endif
        unsigned assigned_to; // keeps track of first assignment to parities
        gather_status move_status[3][NDIRS]; // is communication done

        // neighbour pointers - because of boundary conditions, can be different for
        // diff. fields
        const unsigned *RESTRICT neighbours[NDIRS];
        BoundaryCondition boundary_condition[NDIRS];

#ifdef USE_MPI
        MPI_Request receive_request[3][NDIRS];
        MPI_Request send_request[3][NDIRS];
#ifndef VANILLA
        // vanilla needs no special receive buffers
        T *receive_buffer[NDIRS];
#endif
        T *send_buffer[NDIRS];

        void initialize_communication() {
            for (int d = 0; d < NDIRS; d++) {
                for (int p = 0; p < 3; p++)
                    move_status[p][d] = gather_status::NOT_DONE;
                send_buffer[d] = nullptr;
#ifndef VANILLA
                receive_buffer[d] = nullptr;
#endif
            }
        }

        void free_communication() {
            for (int d = 0; d < NDIRS; d++) {
                if (send_buffer[d] != nullptr)
                    payload.free_mpi_buffer(send_buffer[d]);
#ifndef VANILLA
                if (receive_buffer[d] != nullptr)
                    payload.free_mpi_buffer(receive_buffer[d]);
#endif
            }
        }

#else // Not MPI

        // empty stubs
        void initialize_communication() {}
        void free_communication() {}

#endif

        void allocate_payload() {
            payload.allocate_field(lattice);
        }
        void free_payload() {
            payload.free_field();
        }

#ifndef VECTORIZED
        /// Getter for an individual elements in a loop
        inline auto get(const unsigned i) const {
            return payload.get(i, lattice->field_alloc_size());
        }

        template <typename A>
        inline void set(const A &value, const unsigned i) {
            payload.set(value, i, lattice->field_alloc_size());
        }

        /// Getter for an element outside a loop. Used to manipulate the field directly
        /// outside loops.
        inline auto get_element(const unsigned i) const {
            return payload.get_element(i, lattice);
        }

        template <typename A>
        inline void set_element(const A &value, const unsigned i) {
            payload.set_element(value, i, lattice);
        }
#else
        template <typename vecT>
        inline vecT get_vector(const unsigned i) const {
            return payload.template get_vector<vecT>(i);
        }
        inline T get_element(const unsigned i) const {
            return payload.get_element(i);
        }

        template <typename vecT>
        inline void set_vector(const vecT &val, const unsigned i) {
            return payload.set_vector(val, i);
        }
        inline void set_element(const T &val, const unsigned i) {
            return payload.set_element(val, i);
        }
#endif

        /// Gather boundary elements for communication
        void
        gather_comm_elements(Direction d, Parity par, T *RESTRICT buffer,
                             const lattice_struct::comm_node_struct &to_node) const {
#ifndef VECTORIZED
#ifdef SPECIAL_BOUNDARY_CONDITIONS
            // note: -d in is_on_edge, because we're about to send stuff to that
            // Direction (gathering from Direction +d)
            if (boundary_condition[d] == BoundaryCondition::ANTIPERIODIC &&
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
#ifdef SPECIAL_BOUNDARY_CONDITIONS
            if (boundary_condition[d] == BoundaryCondition::ANTIPERIODIC &&
                lattice->special_boundaries[-d].is_on_edge) {
                antiperiodic = true;
            }
#endif

            if constexpr (hila::is_vectorizable_type<T>::value) {
                // now vectorized layout
                if (vector_lattice->is_boundary_permutation[abs(d)]) {
                    // with boundary permutation need to gather elems 1-by-1
                    int n;
                    const unsigned *index_list = to_node.get_sitelist(par, n);
                    if (!antiperiodic) {
                        payload.gather_elements(buffer, index_list, n, lattice);
                    } else {
                        payload.gather_elements_negated(buffer, index_list, n, lattice);
                    }
                } else {
                    // without it, can do the full block
                    payload.gather_comm_vectors(buffer, to_node, par, vector_lattice,
                                                antiperiodic);
                }
            } else {
                // not vectoizable, standard methods
                int n;
                const unsigned *index_list = to_node.get_sitelist(par, n);
                if (!antiperiodic)
                    payload.gather_elements(buffer, index_list, n, lattice);
                else {
                    payload.gather_elements_negated(buffer, index_list, n, lattice);
                }
            }
#endif
        }

        /// Place boundary elements from neighbour
        void place_comm_elements(Direction d, Parity par, T *RESTRICT buffer,
                                 const lattice_struct::comm_node_struct &from_node) {
// #ifdef USE_MPI
#ifdef VECTORIZED
            if constexpr (hila::is_vectorizable_type<T>::value) {
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
        void set_local_boundary_elements(Direction dir, Parity par) {

#ifdef SPECIAL_BOUNDARY_CONDITIONS
            bool antiperiodic =
                (boundary_condition[dir] == BoundaryCondition::ANTIPERIODIC &&
                 lattice->special_boundaries[dir].is_on_edge);
#else
            bool antiperiodic = false;
#endif
            payload.set_local_boundary_elements(dir, par, lattice, antiperiodic);
        }

        /// Gather a list of elements to a single node
        void gather_elements(T *buffer, const std::vector<CoordinateVector> &coord_list,
                             int root = 0) const;
        void send_elements(T *buffer, const std::vector<CoordinateVector> &coord_list,
                           int root = 0);

#if defined(USE_MPI)

        /// get the receive buffer pointer for the communication.
        T *get_receive_buffer(Direction d, Parity par,
                              const lattice_struct::comm_node_struct &from_node) {
#if defined(VANILLA)

            return (T *)payload.get_buffer() + from_node.offset(par);

#elif defined(CUDA) || defined(HIP)

            unsigned offs = 0;
            if (par == ODD)
                offs = from_node.sites / 2;
            if (receive_buffer[d] == nullptr) {
                receive_buffer[d] = payload.allocate_mpi_buffer(from_node.sites);
            }
            return receive_buffer[d] + offs;

#elif defined(VECTORIZED)

            if constexpr (!hila::is_vectorizable_type<T>::value) {
                // use vanilla type, field laid out in std fashion
                return (T *)payload.get_buffer() + from_node.offset(par);
            } else {
                unsigned offs = 0;
                if (par == ODD)
                    offs = from_node.sites / 2;

                if (vector_lattice->is_boundary_permutation[abs(d)]) {
                    // extra copy operation needed
                    if (receive_buffer[d] == nullptr) {
                        receive_buffer[d] =
                            payload.allocate_mpi_buffer(from_node.sites);
                    }
                    return receive_buffer[d] + offs;
                } else {
                    // directly to halo buffer
                    constexpr unsigned vector_size = hila::vector_info<T>::vector_size;
                    return ((T *)payload.get_buffer() +
                            (vector_lattice->halo_offset[d] * vector_size + offs));
                }
            }
#endif
        } // end of get_receive_buffer
#endif    // USE_MPI
    };

    // static_assert( std::is_pod<T>::value, "Field expects only pod-type elements
    // (plain data): default constructor, copy and delete");
    static_assert(std::is_trivial<T>::value && std::is_standard_layout<T>::value,
                  "Field expects only pod-type elements (plain data): default "
                  "constructor, copy and delete");

  public:
    ////////////////////////////////////////////////
    /// Field::fs keeps all of the field content
    ////////////////////////////////////////////////

    field_struct *RESTRICT fs;

    ////////////////////////////////////////////////
    /// Field constructors

    Field() {
        fs = nullptr; // lazy allocation on 1st use
    }

    // Straightforward copy constructor seems to be necessary
    Field(const Field &other) {
        fs = nullptr; // this is probably unnecessary
        if (other.fs != nullptr) {
            (*this)[ALL] = other[X];
        }
    }

    // copy constructor - from fields which can be assigned
    template <typename A, std::enable_if_t<std::is_convertible<A, T>::value, int> = 0>
    Field(const Field<A> &other) {
        fs = nullptr; // this is probably unnecessary
        if (other.fs != nullptr) {
            (*this)[ALL] = other[X];
        }
    }

    // constructor with compatible scalar
    template <typename A, std::enable_if_t<hila::is_assignable<T &, A>::value ||
                                               std::is_convertible<A, T>::value,
                                           int> = 0>
    Field(const A &val) {
        fs = nullptr;
        (*this)[ALL] = val;
    }

    // constructor from 0 - nullptr trick in use
    Field(const std::nullptr_t z) {
        fs = nullptr;
        (*this)[ALL] = 0;
    }

    // move constructor - steal the content
    Field(Field &&rhs) {
        fs = rhs.fs;
        rhs.fs = nullptr;
    }

    /////////////////////////////////////////////////
    /// Destructor

    ~Field() {
        free();

#ifdef HILAPP
        // Because destructor is instantiated for all fields,
        // use this to put in a hook for generating call to this function
        // in preprocessing stage
        ensure_field_operators_exist(*this);
#endif
    }

    void allocate() {
        assert(fs == nullptr);
        if (lattice == nullptr) {
            output0 << "Can not allocate Field variables before lattice.setup()\n";
            hila::terminate(0);
        }
        fs = (field_struct *)memalloc(sizeof(field_struct));
        fs->lattice = lattice;
        fs->allocate_payload();
        fs->initialize_communication();
        mark_changed(ALL);   // guarantees communications will be done
        fs->assigned_to = 0; // and this means that it is not assigned

        for (Direction d = (Direction)0; d < NDIRS; ++d) {

#if !defined(CUDA) && !defined(HIP)
            fs->neighbours[d] = lattice->neighb[d];
#else
            fs->payload.neighbours[d] = lattice->backend_lattice->d_neighb[d];
#endif
        }

#ifdef SPECIAL_BOUNDARY_CONDITIONS
        foralldir(dir) {
            fs->boundary_condition[dir] = BoundaryCondition::PERIODIC;
            fs->boundary_condition[-dir] = BoundaryCondition::PERIODIC;
        }
#endif

#ifdef VECTORIZED
        fs->vector_lattice =
            lattice->backend_lattice
                ->get_vectorized_lattice<hila::vector_info<T>::vector_size>();
#endif
    }

    void free() {
        // don't call destructors when exiting - either MPI or cuda can already
        // be off.
        if (fs != nullptr && !hila::about_to_finish) {
            for (Direction d = (Direction)0; d < NDIRS; ++d)
                drop_comms(d, ALL);
            fs->free_payload();
            fs->free_communication();
            std::free(fs);
            fs = nullptr;
        }
    }

    bool is_allocated() const {
        return (fs != nullptr);
    }

    bool is_initialized(Parity p) const {
        return fs != nullptr && ((fs->assigned_to & parity_bits(p)) != 0);
    }

    gather_status move_status(Parity p, int d) const {
        assert(parity_bits(p) && d >= 0 && d < NDIRS);
        return fs->move_status[(int)p - 1][d];
    }
    void set_move_status(Parity p, int d, gather_status stat) const {
        assert(parity_bits(p) && d >= 0 && d < NDIRS);
        fs->move_status[(int)p - 1][d] = stat;
    }

    /// check that Field is allocated, and if not do it (if not const)
    /// Must be called BEFORE the var is actually used
    /// "hilapp" will generate these calls as needed!
    void check_alloc() {
        if (!is_allocated())
            allocate();
    }

    /// If Field is const specified, we should not be able to write to it in the first
    /// place
    void check_alloc() const {
        assert(is_allocated());
    }

    // If ALL changes, both parities invalid; if p != ALL, then p and ALL.
    void mark_changed(const Parity p) const {

        for (Direction i = (Direction)0; i < NDIRS; ++i) {
            // check if there's ongoing comms, invalidate it!
            drop_comms(i, opp_parity(p));

            set_move_status(opp_parity(p), i, gather_status::NOT_DONE);
            if (p != ALL) {
                set_move_status(ALL, i, gather_status::NOT_DONE);
            } else {
                set_move_status(EVEN, i, gather_status::NOT_DONE);
                set_move_status(ODD, i, gather_status::NOT_DONE);
            }
        }
        fs->assigned_to |= parity_bits(p);
    }

    /// Mark the field parity gathered from Direction
    // In case p=ALL we could mark everything gathered, but we'll be conservative here
    // and mark only this parity, because there might be other parities on the fly and
    // corresponding waits should be done,  This should never happen in automatically
    // generated loops. In any case start_gather, is_gathered, get_move_parity has
    // intelligence to figure out the right thing to do
    //

    void mark_gathered(int dir, const Parity p) const {
        set_move_status(p, dir, gather_status::DONE);
    }

    // Check if the field has been gathered since the previous communication
    // par = ALL:   ALL or (EVEN+ODD) are OK
    // par != ALL:  ALL or par are OK
    bool is_gathered(int dir, Parity par) const {
        if (par != ALL) {
            return move_status(par, dir) == gather_status::DONE ||
                   move_status(ALL, dir) == gather_status::DONE;
        } else {
            return move_status(ALL, dir) == gather_status::DONE ||
                   (move_status(EVEN, dir) == gather_status::DONE &&
                    move_status(ODD, dir) == gather_status::DONE);
        }
    }

    // Mark communication started -- this must be just the one
    // going on with MPI
    void mark_move_started(int dir, Parity p) const {
        set_move_status(p, dir, gather_status::STARTED);
    }

    /// Check if communication has started.  This is strict, checks exactly this parity
    bool is_move_started(int dir, Parity par) const {
        return move_status(par, dir) == gather_status::STARTED;
    }

    bool move_not_done(int dir, Parity par) const {
        return move_status(par, dir) == gather_status::NOT_DONE;
    }

    void set_boundary_condition(Direction dir, BoundaryCondition bc) {

#ifdef SPECIAL_BOUNDARY_CONDITIONS
        // TODO: This works as intended only for periodic/antiperiodic b.c.
        check_alloc();
        fs->boundary_condition[dir] = bc;
        fs->boundary_condition[-dir] = bc;
#if !defined(CUDA) && !defined(HIP)
        fs->neighbours[dir] = lattice->get_neighbour_array(dir, bc);
        fs->neighbours[-dir] = lattice->get_neighbour_array(-dir, bc);
#else
        if (bc == BoundaryCondition::PERIODIC) {
            fs->payload.neighbours[dir] = lattice->backend_lattice->d_neighb[dir];
            fs->payload.neighbours[-dir] = lattice->backend_lattice->d_neighb[-dir];
        } else {
            fs->payload.neighbours[dir] =
                lattice->backend_lattice->d_neighb_special[dir];
            fs->payload.neighbours[-dir] =
                lattice->backend_lattice->d_neighb_special[-dir];
        }
#endif

        // Make sure boundaries get refreshed
        mark_changed(ALL);
#endif
    }

    BoundaryCondition get_boundary_condition(Direction dir) const {
#ifdef SPECIAL_BOUNDARY_CONDITIONS
        return fs->boundary_condition[dir];
#else
        return BoundaryCondition::PERIODIC;
#endif
    }

    void print_boundary_condition() {
        check_alloc();
        output0 << " ( ";
        for (int dir = 0; dir < NDIRS; dir++) {
            output0 << (int)fs->boundary_condition[dir] << " ";
        }
        output0 << ")\n";
    }

    template <typename A>
    void copy_boundary_condition(const Field<A> &rhs) {
        foralldir(dir) {
            set_boundary_condition(dir, rhs.get_boundary_condition(dir));
        }
    }

    // Overloading []
    // declarations -- WILL BE implemented by hilapp, not written here
    // let there be const and non-const protos
    T operator[](const Parity p) const;           // f[EVEN]
    T operator[](const X_index_type) const;       // f[X]
    T operator[](const X_plus_direction p) const; // f[X+dir]
    T operator[](const X_plus_offset p) const;    // f[X+dir1+dir2] and others

    T &operator[](const Parity p);     // f[EVEN]
    T &operator[](const X_index_type); // f[X]

    T &operator[](const CoordinateVector &v);       // f[CoordinateVector]
    T &operator[](const CoordinateVector &v) const; // f[CoordinateVector]

    // TEMPORARY HACK: return ptr to bare array
    inline auto field_buffer() const {
        return fs->payload.get_buffer();
    }

#ifndef VECTORIZED
    /// Get an individual element outside a loop. This is also used as a getter in the
    /// vanilla code.
    inline auto get_value_at(const unsigned i) const {
        return fs->get_element(i);
    }
#else
    inline auto get_value_at(const unsigned i) const {
        return fs->get_element(i);
    }
    template <typename vecT>
    inline auto get_vector_at(unsigned i) const {
        return fs->template get_vector<vecT>(i);
    }
    inline auto get_value_at_nb_site(Direction d, unsigned i) const {
        return fs->get_element(fs->vector_lattice->site_neighbour(d, i));
    }
#endif

#ifndef VECTORIZED
    /// Set an individual element outside a loop. This is also used as a setter in the
    /// vanilla code.
    template <typename A>
    inline void set_value_at(const A &value, unsigned i) {
        fs->set_element(value, i);
    }

#else
    template <typename vecT>
    inline void set_vector_at(const vecT &value, unsigned i) {
        fs->set_vector(value, i);
    }

    template <typename A>
    inline void set_value_at(const A &value, unsigned i) {
        fs->set_element(value, i);
    }
#endif

    /////////////////////////////////////////////////////////////////
    /// Standard arithmetic ops which fields should implement
    /// Not all are always callable, e.g. division may not be
    /// implemented by all field types
    /////////////////////////////////////////////////////////////////

    // Basic assignment operator
    Field<T> &operator=(const Field<T> &rhs) {
        (*this)[ALL] = rhs[X];
        return *this;
    }

    // More general = - possible only if T = A is OK
    template <typename A, std::enable_if_t<hila::is_assignable<T &, A>::value ||
                                               std::is_convertible<A, T>::value,
                                           int> = 0>
    Field<T> &operator=(const Field<A> &rhs) {
        (*this)[ALL] = rhs[X];
        return *this;
    }

    // Assign from element
    template <typename A, std::enable_if_t<hila::is_assignable<T &, A>::value ||
                                               std::is_convertible<A, T>::value,
                                           int> = 0>
    Field<T> &operator=(const A &d) {
        (*this)[ALL] = d;
        return *this;
    }

    // assignment of 0 - nullptr, zeroes field
    Field<T> &operator=(const std::nullptr_t &z) {
        (*this)[ALL] = 0;
        return *this;
    }

    // Do also move assignment
    Field<T> &operator=(Field<T> &&rhs) {
        if (this != &rhs) {
            free();
            fs = rhs.fs;
            rhs.fs = nullptr;
        }
        return *this;
    }

    // +=, -=  etc operators from compatible types
    template <
        typename A,
        std::enable_if_t<std::is_convertible<hila::type_plus<T, A>, T>::value, int> = 0>
    Field<T> &operator+=(const Field<A> &rhs) {
        (*this)[ALL] += rhs[X];
        return *this;
    }

    template <typename A,
              std::enable_if_t<std::is_convertible<hila::type_minus<T, A>, T>::value,
                               int> = 0>
    Field<T> &operator-=(const Field<A> &rhs) {
        (*this)[ALL] -= rhs[X];
        return *this;
    }

    template <
        typename A,
        std::enable_if_t<std::is_convertible<hila::type_mul<T, A>, T>::value, int> = 0>
    Field<T> &operator*=(const Field<A> &rhs) {
        (*this)[ALL] *= rhs[X];
        return *this;
    }

    template <
        typename A,
        std::enable_if_t<std::is_convertible<hila::type_div<T, A>, T>::value, int> = 0>
    Field<T> &operator/=(const Field<A> &rhs) {
        (*this)[ALL] /= rhs[X];
        return *this;
    }

    template <
        typename A,
        std::enable_if_t<std::is_convertible<hila::type_plus<T, A>, T>::value, int> = 0>
    Field<T> &operator+=(const A &rhs) {
        (*this)[ALL] += rhs;
        return *this;
    }

    template <typename A,
              std::enable_if_t<std::is_convertible<hila::type_minus<T, A>, T>::value,
                               int> = 0>
    Field<T> &operator-=(const A &rhs) {
        (*this)[ALL] -= rhs;
        return *this;
    }

    template <
        typename A,
        std::enable_if_t<std::is_convertible<hila::type_mul<T, A>, T>::value, int> = 0>
    Field<T> &operator*=(const A &rhs) {
        (*this)[ALL] *= rhs;
        return *this;
    }

    template <
        typename A,
        std::enable_if_t<std::is_convertible<hila::type_div<T, A>, T>::value, int> = 0>
    Field<T> &operator/=(const A &rhs) {
        (*this)[ALL] /= rhs;
        return *this;
    }

    // Unary + and -
    Field<T> operator+() const {
        return *this;
    }

    Field<T> operator-() const {
        Field<T> f;
        f[ALL] = -(*this)[X];
        return f;
    }

    hila::number_type<T> squarenorm() const {
        hila::number_type<T> n = 0;
        onsites(ALL) {
            n += ::squarenorm((*this)[X]);
        }
        return n;
    }

    ///////////////////////////////////////////////////////////////////////

    // Communication routines
    dir_mask_t start_gather(Direction d, Parity p = ALL) const;
    void wait_gather(Direction d, Parity p) const;
    void gather(Direction d, Parity p = ALL) const;
    void drop_comms(Direction d, Parity p) const;
    void cancel_comm(Direction d, Parity p) const;

    // Declaration of shift methods
    Field<T> &shift(const CoordinateVector &v, Field<T> &r, Parity par) const;
    Field<T> &shift(const CoordinateVector &v, Field<T> &r) const {
        return shift(v, r, ALL);
    }
    Field<T> shift(const CoordinateVector &v, Parity par) const;

    // General getters and setters
    void set_elements(T *elements, const std::vector<CoordinateVector> &coord_list);
    void set_element(const T &element, const CoordinateVector &coord);
    void get_elements(T *elements,
                      const std::vector<CoordinateVector> &coord_list) const;
    T get_element(const CoordinateVector &coord) const;

    template <typename A, std::enable_if_t<std::is_assignable<T &, A>::value, int> = 0>
    inline void set_element_at(const CoordinateVector coord, const A elem) {
        T e;
        e = elem;
        set_element(e, coord);
    }

    inline void set_element_at(const CoordinateVector coord, std::nullptr_t elem) {
        T e;
        e = 0;
        set_element(e, coord);
    }

    // Fourier transform declarations
    void FFT(fft_direction fdir = fft_direction::forward);

    // Writes the Field to disk
    void write_to_stream(std::ofstream &outputfile);
    void write_to_file(const std::string &filename);
    void read_from_stream(std::ifstream &inputfile);
    void read_from_file(const std::string &filename);

    void write_subvolume(std::ofstream &outputfile, const CoordinateVector &cmin,
                         const CoordinateVector &cmax, int precision = 6);
    void write_subvolume(const std::string &filenname, const CoordinateVector &cmin,
                         const CoordinateVector &cmax, int precision = 6);

    void write_slice(std::ofstream &outputfile, const CoordinateVector &slice,
                     int precision = 6);
    void write_slice(const std::string &outputfile, const CoordinateVector &slice,
                     int precision = 6);

    // and sum reduction
    T sum(Parity par = Parity::all, bool allreduce = true) const;

    T product(Parity par = Parity::all, bool allreduce = true) const;

    // Declare gpu_reduce here, defined only for GPU targets
    // For internal use only, preferably
    // T gpu_reduce_sum(bool allreduce = true, Parity par = Parity::all,
    //              bool do_mpi = true) const;
    /// Declare gpu_reduce here, defined only for GPU targets
    /// For internal use only, preferably

    T gpu_minmax(bool min_or_max, Parity par, CoordinateVector &loc) const;

    T min(Parity par = ALL) const;
    T min(CoordinateVector &loc) const;
    T min(Parity par, CoordinateVector &loc) const;
    T max(Parity par = ALL) const;
    T max(CoordinateVector &loc) const;
    T max(Parity par, CoordinateVector &loc) const;
    T minmax(bool is_min, Parity par, CoordinateVector &loc) const;

}; // End of class Field<>

// these operators rely on SFINAE, OK if field_hila::type_plus<A,B> exists i.e. A+B is
// OK
/// operator +
template <typename A, typename B>
auto operator+(Field<A> &lhs, Field<B> &rhs) -> Field<hila::type_plus<A, B>> {
    Field<hila::type_plus<A, B>> tmp;
    tmp[ALL] = lhs[X] + rhs[X];
    return tmp;
}

template <typename A, typename B>
auto operator+(const A &lhs, const Field<B> &rhs) -> Field<hila::type_plus<A, B>> {
    Field<hila::type_plus<A, B>> tmp;
    tmp[ALL] = lhs + rhs[X];
    return tmp;
}

template <typename A, typename B>
auto operator+(const Field<A> &lhs, const B &rhs) -> Field<hila::type_plus<A, B>> {
    Field<hila::type_plus<A, B>> tmp;
    tmp[ALL] = lhs[X] + rhs;
    return tmp;
}

/// operator -
template <typename A, typename B>
auto operator-(const Field<A> &lhs, const Field<B> &rhs)
    -> Field<hila::type_minus<A, B>> {
    Field<hila::type_minus<A, B>> tmp;
    tmp[ALL] = lhs[X] - rhs[X];
    return tmp;
}

template <typename A, typename B>
auto operator-(const A &lhs, const Field<B> &rhs) -> Field<hila::type_minus<A, B>> {
    Field<hila::type_minus<A, B>> tmp;
    tmp[ALL] = lhs - rhs[X];
    return tmp;
}

template <typename A, typename B>
auto operator-(const Field<A> &lhs, const B &rhs) -> Field<hila::type_minus<A, B>> {
    Field<hila::type_minus<A, B>> tmp;
    tmp[ALL] = lhs[X] - rhs;
    return tmp;
}

/// operator *
template <typename A, typename B>
auto operator*(const Field<A> &lhs, const Field<B> &rhs)
    -> Field<hila::type_mul<A, B>> {
    Field<hila::type_mul<A, B>> tmp;
    tmp[ALL] = lhs[X] * rhs[X];
    return tmp;
}

template <typename A, typename B>
auto operator*(const A &lhs, const Field<B> &rhs) -> Field<hila::type_mul<A, B>> {
    Field<hila::type_mul<A, B>> tmp;
    tmp[ALL] = lhs * rhs[X];
    return tmp;
}

template <typename A, typename B>
auto operator*(const Field<A> &lhs, const B &rhs) -> Field<hila::type_mul<A, B>> {
    Field<hila::type_mul<A, B>> tmp;
    tmp[ALL] = lhs[X] * rhs;
    return tmp;
}

/// operator /
template <typename A, typename B>
auto operator/(const Field<A> &lhs, const Field<B> &rhs)
    -> Field<hila::type_div<A, B>> {
    Field<hila::type_div<A, B>> tmp;
    tmp[ALL] = lhs[X] / rhs[X];
    return tmp;
}

template <typename A, typename B>
auto operator/(const A &lhs, const Field<B> &rhs) -> Field<hila::type_div<A, B>> {
    Field<hila::type_div<A, B>> tmp;
    tmp[ALL] = lhs / rhs[X];
    return tmp;
}

template <typename A, typename B>
auto operator/(const Field<A> &lhs, const B &rhs) -> Field<hila::type_div<A, B>> {
    Field<hila::type_div<A, B>> tmp;
    tmp[ALL] = lhs[X] / rhs;
    return tmp;
}


#define NAIVE_SHIFT
#if defined(NAIVE_SHIFT)

/// Definition of shift - this is currently OK only for short moves,
/// very inefficient for longer moves
/// TODO: make more advanced, switching to "global" move for long shifts
/// Returns a reference to parameter "res"

template <typename T>
Field<T> &Field<T>::shift(const CoordinateVector &v, Field<T> &res,
                          const Parity par) const {

    // use this to store remaining moves
    CoordinateVector rem = v;

    // check the parity of the move
    Parity par_s;

    int len = 0;
    foralldir(d) len += abs(rem[d]);

    // no move, just copy field
    if (len == 0) {
        res = *this;
        return res;
    }

    // opp_parity(ALL) == ALL
    if (len % 2 == 0)
        par_s = opp_parity(par);
    else
        par_s = par;

    // is this already gathered from one of the dirs in v?
    bool found_dir = false;
    Direction mdir;
    foralldir(d) {
        if (rem[d] > 0 && move_status(par_s, d) != gather_status::NOT_DONE) {
            mdir = d;
            found_dir = true;
            break;
        } else if (rem[d] < 0 && move_status(par_s, -d) != gather_status::NOT_DONE) {
            mdir = -d;
            found_dir = true;
            break;
        }
    }

    if (!found_dir) {
        // now did not find a 'ready' dir. Take the 1st available
        foralldir(d) {
            if (rem[d] > 0) {
                mdir = d;
                break;
            } else if (rem[d] < 0) {
                mdir = -d;
                break;
            }
        }
    }

    // Len 1, copy directly
    if (len == 1) {
        res[par_s] = (*this)[X + mdir];
        return res;
    }

    // now longer - need buffer
    Field<T> r1;
    Field<T> *from, *to;

    // this ensures that the final move lands on res
    if (len % 2 == 0) {
        from = &r1;
        to = &res;
    } else {
        from = &res;
        to = &r1;
    }
    // and copy initially to "from"
    (*from)[par_s] = (*this)[X + mdir];

    // and subtract remaining moves from rem
    rem = rem - mdir;
    par_s = opp_parity(par_s);

    foralldir(d) {
        if (rem[d] != 0) {
            mdir = (rem[d] > 0) ? d : -d;

            while (rem[d] != 0) {

                (*to)[par_s] = (*from)[X + mdir];

                par_s = opp_parity(par_s);
                rem = rem - mdir;
                std::swap(to, from);
            }
        }
    }

    return res;
}

template <typename T>
Field<T> Field<T>::shift(const CoordinateVector &v, const Parity par) const {
    Field<T> res;
    shift(v, res, par);
    return res;
}
#elif !defined(USE_MPI)

// this is junk at the moment
template <typename T>
Field<T> &Field<T>::shift(const CoordinateVector &v, Field<T> &res,
                          const Parity par) const {

    onsites(par) {
        if
    }
    r2 = *this;
    foralldir(d) {
        if (abs(v[d]) > 0) {
            Direction dir;
            if (v[d] > 0)
                dir = d;
            else
                dir = -d;

            for (int i = 0; i < abs(v[d]); i++) {
                r1[ALL] = r2[X + dir];
                r2 = r1;
            }
        }
    }
    return r2;
}

#endif

#if defined(USE_MPI)

/// start_gather(): Communicate the field at Parity par from Direction
/// d. Uses accessors to prevent dependency on the layout.
/// return the Direction mask bits where something is happening
template <typename T>
dir_mask_t Field<T>::start_gather(Direction d, Parity p) const {

    // get the mpi message tag right away, to ensure that we are always synchronized
    // with the mpi calls -- some nodes might not need comms, but the tags must be in
    // sync

    int tag = get_next_msg_tag();

    lattice_struct::nn_comminfo_struct &ci = lattice->nn_comminfo[d];
    lattice_struct::comm_node_struct &from_node = ci.from_node;
    lattice_struct::comm_node_struct &to_node = ci.to_node;

    // check if this is done - either gathered or no comm to be done in the 1st place

    if (is_gathered(d, p)) {
        lattice->n_gather_avoided++;
        return 0; // nothing to wait for
    }

    // No comms to do, nothing to wait for -- we'll use the is_gathered
    // status to keep track of vector boundary shuffle anyway

    if (from_node.rank == hila::myrank() && to_node.rank == hila::myrank()) {
        fs->set_local_boundary_elements(d, p);
        mark_gathered(d, p);
        return 0;
    }

    // if this parity or ALL-type gather is going on nothing to be done
    if (!move_not_done(d, p) || !move_not_done(d, ALL)) {
        lattice->n_gather_avoided++;
        return get_dir_mask(d); // nothing to do, but still need to wait
    }

    Parity par = p;
    // if p is ALL but ODD or EVEN is going on/done, turn off parity which is not needed
    // corresponding wait must do the same thing
    if (p == ALL) {
        if (!move_not_done(d, EVEN) && !move_not_done(d, ODD)) {
            // even and odd are going on or ready, nothing to be done
            lattice->n_gather_avoided++;
            return get_dir_mask(d);
        }
        if (!move_not_done(d, EVEN))
            par = ODD;
        else if (!move_not_done(d, ODD))
            par = EVEN;
        // if neither is the case par = ALL
    }

    mark_move_started(d, par);

    // Communication hasn't been started yet, do it now

    int par_i = static_cast<int>(par) - 1; // index to dim-3 arrays

    constexpr size_t size = sizeof(T);

    T *receive_buffer;
    T *send_buffer;

    int size_type;
    MPI_Datatype mpi_type = get_MPI_number_type<T>(size_type);

    if (from_node.rank != hila::myrank()) {

        // HANDLE RECEIVES: get node which will send here
        post_receive_timer.start();

        // buffer can be separate or in Field buffer
        receive_buffer = fs->get_receive_buffer(d, par, from_node);

        size_t n = from_node.n_sites(par) * size / size_type;

        if (n >= (1ULL << 31)) {
            hila::output << "Too large MPI message!  Size " << n << '\n';
            hila::terminate(1);
        }

        // c++ version does not return errors
        MPI_Irecv(receive_buffer, n, mpi_type, from_node.rank, tag,
                  lattice->mpi_comm_lat, &fs->receive_request[par_i][d]);

        post_receive_timer.stop();
    }

    if (to_node.rank != hila::myrank()) {
        // HANDLE SENDS: Copy Field elements on the boundary to a send buffer and send
        start_send_timer.start();

        unsigned sites = to_node.n_sites(par);

        if (fs->send_buffer[d] == nullptr)
            fs->send_buffer[d] = fs->payload.allocate_mpi_buffer(to_node.sites);
        send_buffer = fs->send_buffer[d] + to_node.offset(par);

        fs->gather_comm_elements(d, par, send_buffer, to_node);

        size_t n = sites * size / size_type;

        MPI_Isend(send_buffer, n, mpi_type, to_node.rank, tag, lattice->mpi_comm_lat,
                  &fs->send_request[par_i][d]);
        start_send_timer.stop();
    }

    // and do the boundary shuffle here, after MPI has started
    // NOTE: there should be no danger of MPI and shuffle overwriting, MPI writes
    // to halo buffers only if no permutation is needed.  With a permutation MPI
    // uses special receive buffer
    fs->set_local_boundary_elements(d, par);

    return get_dir_mask(d);
}

///  wait_gather(): Wait for communication at parity par from
///  Direction d completes the communication in the function.
///  If the communication has not started yet, also calls
///  start_gather()
///
///  NOTE: This will be called even if the field is marked const.
///  Therefore this function is const, even though it does change
///  the internal content of the field, the halo. From the point
///  of view of the user, the value of the field does not change.
template <typename T>
void Field<T>::wait_gather(Direction d, Parity p) const {

    lattice_struct::nn_comminfo_struct &ci = lattice->nn_comminfo[d];
    lattice_struct::comm_node_struct &from_node = ci.from_node;
    lattice_struct::comm_node_struct &to_node = ci.to_node;

    // check if this is done - either gathered or no comm to be done in the 1st place
    if (is_gathered(d, p))
        return;

    // this is the branch if no comms -- shuffle was done in start_gather
    if (from_node.rank == hila::myrank() && to_node.rank == hila::myrank())
        return;

    // if (!is_move_started(d,p)) {
    //   output0 << "Wait move error - wait_gather without corresponding
    //   start_gather\n"; exit(1);
    // }

    // Note: the move can be Parity p OR ALL -- need to wait for it in any case
    // set par to be the "sum" over both parities
    // There never should be ongoing ALL and other parity gather -- start_gather takes
    // care

    // check here consistency, this should never happen
    if (p != ALL && is_move_started(d, p) && is_move_started(d, ALL)) {
        exit(1);
    }

    Parity par;
    int n_wait = 1;
    // what par to wait for?
    if (is_move_started(d, p))
        par = p; // standard match
    else if (p != ALL) {
        if (is_move_started(d, ALL))
            par = ALL; // if all is running wait for it
        else {
            exit(1);
        }
    } else {
        // now p == ALL and ALL is not running
        if (is_gathered(d, EVEN) && is_move_started(d, ODD))
            par = ODD;
        else if (is_gathered(d, ODD) && is_move_started(d, EVEN))
            par = EVEN;
        else if (is_move_started(d, EVEN) && is_move_started(d, ODD)) {
            n_wait = 2; // need to wait for both!
            par = ALL;
        } else {
            exit(1);
        }
    }

    if (n_wait == 2)
        par = EVEN; // we'll flip both

    for (int wait_i = 0; wait_i < n_wait; ++wait_i) {

        int par_i = (int)par - 1;

        if (from_node.rank != hila::myrank()) {
            wait_receive_timer.start();

            MPI_Status status;
            MPI_Wait(&fs->receive_request[par_i][d], &status);

            wait_receive_timer.stop();

#ifndef VANILLA
            fs->place_comm_elements(d, par, fs->get_receive_buffer(d, par, from_node),
                                    from_node);
#endif
        }

        // then wait for the sends
        if (to_node.rank != hila::myrank()) {
            wait_send_timer.start();
            MPI_Status status;
            MPI_Wait(&fs->send_request[par_i][d], &status);
            wait_send_timer.stop();
        }

        // Mark the parity gathered from Direction dir
        mark_gathered(d, par);

        // Keep count of communications
        lattice->n_gather_done += 1;

        par = opp_parity(par); // flip if 2 loops
    }
}

///  drop_comms():  if field is changed or deleted,
///  cancel ongoing communications.  This should happen very seldom,
///  only if there are "by-hand" start_gather operations and these are not needed
template <typename T>
void Field<T>::drop_comms(Direction d, Parity p) const {

    if (is_comm_initialized()) {
        if (is_move_started(d, ALL))
            cancel_comm(d, ALL);
        if (p != ALL) {
            if (is_move_started(d, p))
                cancel_comm(d, p);
        } else {
            if (is_move_started(d, EVEN))
                cancel_comm(d, EVEN);
            if (is_move_started(d, ODD))
                cancel_comm(d, ODD);
        }
    }
}

/// cancel ongoing send and receive

template <typename T>
void Field<T>::cancel_comm(Direction d, Parity p) const {
    if (lattice->nn_comminfo[d].from_node.rank != hila::myrank()) {
        cancel_receive_timer.start();
        MPI_Cancel(&fs->receive_request[(int)p - 1][d]);
        cancel_receive_timer.stop();
    }
    if (lattice->nn_comminfo[d].to_node.rank != hila::myrank()) {
        cancel_send_timer.start();
        MPI_Cancel(&fs->send_request[(int)p - 1][d]);
        cancel_send_timer.stop();
    }
}

#else // No MPI now

///* Trivial implementation when no MPI is used

template <typename T>
dir_mask_t Field<T>::start_gather(Direction d, Parity p) const {
    // Update local elements in the halo (necessary for vectorized version)
    // We use here simpler tracking than in MPI, may lead to slight extra work
    if (!is_gathered(d, p)) {
        fs->set_local_boundary_elements(d, p);
        mark_gathered(d, p);
    }
    return 0;
}

template <typename T>
void Field<T>::wait_gather(Direction d, Parity p) const {}

template <typename T>
void Field<T>::drop_comms(Direction d, Parity p) const {}

#endif // MPI

/// And a convenience combi function
template <typename T>
void Field<T>::gather(Direction d, Parity p) const {
    start_gather(d, p);
    wait_gather(d, p);
}

#if defined(USE_MPI)

/// Gather a list of elements to a single node
/// coord_list must be same on all nodes, buffer is needed only on "root"
template <typename T>
void Field<T>::field_struct::gather_elements(
    T *RESTRICT buffer, const std::vector<CoordinateVector> &coord_list,
    int root) const {

    std::vector<unsigned> index_list;
    std::vector<int> sites_on_rank(lattice->n_nodes());
    std::vector<int> reshuffle_list(coord_list.size());

    std::fill(sites_on_rank.begin(), sites_on_rank.end(), 0);

    int nranks = 0;

    int i = 0;
    for (const CoordinateVector &c : coord_list) {
        int rank = lattice->node_rank(c);
        if (hila::myrank() == rank) {
            index_list.push_back(lattice->site_index(c));
        }

        if (sites_on_rank[rank] == 0 && rank != root)
            nranks++;
        sites_on_rank[rank]++;
        reshuffle_list[i++] = rank;
    }

    std::vector<T> send_buffer(index_list.size());
    payload.gather_elements((T *)send_buffer.data(), index_list.data(),
                            send_buffer.size(), lattice);
    if (hila::myrank() != root && sites_on_rank[hila::myrank()] > 0) {
        MPI_Send((char *)send_buffer.data(), sites_on_rank[hila::myrank()] * sizeof(T),
                 MPI_BYTE, root, hila::myrank(), lattice->mpi_comm_lat);
    }
    if (hila::myrank() == root) {

        // allocate buffer for receiving data
        T *b;
        std::vector<T> pb(coord_list.size() - sites_on_rank[root]);
        b = pb.data();
        // vector for node ptrs -- point to stuff from nodes
        std::vector<T *> nptr(lattice->n_nodes());

        std::vector<MPI_Request> mpi_req(nranks);
        int nreqs = 0;
        for (int n = 0; n < sites_on_rank.size(); n++) {
            if (sites_on_rank[n] > 0) {
                if (n != root) {
                    MPI_Status status;
                    MPI_Irecv(b, sites_on_rank[n] * sizeof(T), MPI_BYTE, n, n,
                              lattice->mpi_comm_lat, &mpi_req[nreqs++]);

                    nptr[n] = b;
                    b += sites_on_rank[n];

                } else {

                    nptr[n] = send_buffer.data();
                }
            }
        }

        if (nreqs > 0) {
            std::vector<MPI_Status> stat_arr(nreqs);
            MPI_Waitall(nreqs, mpi_req.data(), stat_arr.data());
        }

        // copy the data from bp to buffer, reordering
        for (int i = 0; i < coord_list.size(); i++) {
            buffer[i] = *nptr[reshuffle_list[i]];
            nptr[reshuffle_list[i]]++;
        }
    }
}

/// Send elements from a single node to a list of coordinates
/// coord_list must be the same on all nodes, but buffer is needed only on "root"!

template <typename T>
void Field<T>::field_struct::send_elements(
    T *RESTRICT buffer, const std::vector<CoordinateVector> &coord_list, int root) {

    std::vector<unsigned> index_list;
    std::vector<int> sites_on_rank(lattice->n_nodes());
    std::vector<int> reshuffle_list(coord_list.size());
    std::fill(sites_on_rank.begin(), sites_on_rank.end(), 0);

    int nranks = 0;
    int i = 0;
    for (CoordinateVector c : coord_list) {
        int rank = lattice->node_rank(c);
        if (hila::myrank() == rank) {
            index_list.push_back(lattice->site_index(c));
        }

        if (sites_on_rank[rank] == 0 && rank != root)
            nranks++;
        sites_on_rank[rank]++;
        reshuffle_list[i++] = rank;
    }

    // payload.gather_elements((T *)recv_buffer.data(), index_list.data(),
    //                         recv_buffer.size(), lattice);

    if (hila::myrank() != root && sites_on_rank[hila::myrank()] > 0) {
        std::vector<T> recv_buffer(index_list.size());
        MPI_Status status;

        MPI_Recv((char *)recv_buffer.data(), sites_on_rank[hila::myrank()] * sizeof(T),
                 MPI_BYTE, root, hila::myrank(), lattice->mpi_comm_lat, &status);

        payload.place_elements((T *)recv_buffer.data(), index_list.data(),
                               recv_buffer.size(), lattice);
    }
    if (hila::myrank() == root) {
        // reordering buffers
        std::vector<T> pb(coord_list.size());
        // vector for node counters -- point to stuff from nodes
        std::vector<unsigned> nloc(lattice->n_nodes());
        std::vector<unsigned> ncount(lattice->n_nodes());
        nloc[0] = ncount[0] = 0;

        for (int n = 1; n < lattice->n_nodes(); n++) {
            nloc[n] = nloc[n - 1] + sites_on_rank[n - 1];
            ncount[n] = 0;
        }
        for (int i = 0; i < coord_list.size(); i++) {
            int node = reshuffle_list[i];
            pb[nloc[node] + ncount[node]] = buffer[i];
            ncount[node]++;
        }

        std::vector<MPI_Request> mpi_req(nranks);
        int nreqs = 0;
        for (int n = 0; n < sites_on_rank.size(); n++) {
            if (sites_on_rank[n] > 0) {
                if (n != root) {
                    MPI_Isend(pb.data() + nloc[n], sites_on_rank[n] * sizeof(T),
                              MPI_BYTE, n, n, lattice->mpi_comm_lat, &mpi_req[nreqs++]);
                }
            }
        }

        payload.place_elements(pb.data() + nloc[root], index_list.data(),
                               index_list.size(), lattice);

        if (nreqs > 0) {
            std::vector<MPI_Status> stat_arr(nreqs);
            MPI_Waitall(nreqs, mpi_req.data(), stat_arr.data());
        }
    }
}

#else // Now not USE_MPI

/// Gather a list of elements to a single node
template <typename T>
void Field<T>::field_struct::gather_elements(
    T *buffer, const std::vector<CoordinateVector> &coord_list, int root) const {
    std::vector<unsigned> index_list;
    for (CoordinateVector c : coord_list) {
        index_list.push_back(lattice->site_index(c));
    }

    payload.gather_elements(buffer, index_list.data(), index_list.size(), lattice);
}

/// Send elements from a single node to a list of coordinates
template <typename T>
void Field<T>::field_struct::send_elements(
    T *buffer, const std::vector<CoordinateVector> &coord_list, int root) {
    std::vector<unsigned> index_list;
    for (CoordinateVector c : coord_list) {
        index_list.push_back(lattice->site_index(c));
    }

    payload.place_elements(buffer, index_list.data(), index_list.size(), lattice);
}

#endif

/// Functions for manipulating individual elements in an array

/// Set an element. Assuming that each node calls this with the same value, it is
/// sufficient to set the elements locally
template <typename T>
void Field<T>::set_elements(T *elements,
                            const std::vector<CoordinateVector> &coord_list) {
    std::vector<unsigned> my_indexes;
    std::vector<unsigned> my_elements;
    for (int i = 0; i < coord_list.size(); i++) {
        CoordinateVector c = coord_list[i];
        if (lattice->is_on_mynode(c)) {
            my_indexes.push_back(lattice->site_index(c));
            my_elements.push_back(elements[i]);
        }
    }
    fs->payload.place_elements(my_elements.data(), my_indexes.data(), my_indexes.size(),
                               lattice);
    mark_changed(ALL);
}

// Set a single element. Assuming that each node calls this with the same value, it
// is
/// sufficient to set the element locally
template <typename T>
void Field<T>::set_element(const T &element, const CoordinateVector &coord) {
    if (lattice->is_on_mynode(coord)) {
        set_value_at(element, lattice->site_index(coord));
    }
    mark_changed(ALL);
}

/// Get an element and return it on all nodes
#if defined(USE_MPI)
/// This is not local, the element needs to be communicated to all nodes
template <typename T>
T Field<T>::get_element(const CoordinateVector &coord) const {
    T element;

    int owner = lattice->node_rank(coord);

    if (hila::myrank() == owner) {
        element = get_value_at(lattice->site_index(coord));
    }

    MPI_Bcast(&element, sizeof(T), MPI_BYTE, owner, lattice->mpi_comm_lat);
    return element;
}

/// Get a list of elements and store them into an array on all nodes
template <typename T>
void Field<T>::get_elements(T *elements,
                            const std::vector<CoordinateVector> &coord_list) const {
    struct node_site_list_struct {
        std::vector<int> indexes;
        std::vector<CoordinateVector> coords;
    };

    std::vector<node_site_list_struct> nodelist(lattice->n_nodes());
    // Reorganize the list according to nodes
    for (int i = 0; i < coord_list.size(); i++) {
        CoordinateVector c = coord_list[i];
        int node = lattice->node_rank(c);
        nodelist[node].indexes.push_back(i);
        nodelist[node].coords.push_back(c);
    }

    // Fetch on each node found and communicate
    for (int n = 0; n < nodelist.size(); n++) {
        node_site_list_struct node = nodelist[n];
        if (node.indexes.size() > 0) {
            T *element_buffer = (T *)memalloc(sizeof(T) * node.indexes.size());
            fs->payload.gather_elements(element_buffer, node.coords);

            MPI_Bcast(&element_buffer, sizeof(T), MPI_BYTE, n, lattice->mpi_comm_lat);

            // place in the array in original order
            for (int i = 0; i < node.indexes.size(); i++) {
                elements[i] = element_buffer[node.indexes[i]];
            }

            std::free(element_buffer);
        }
    }
}

#else
/// Without MPI, we just need to call get
template <typename T>
T Field<T>::get_element(const CoordinateVector &coord) const {
    return get_value_at(lattice->site_index(coord));
}

/// Without MPI, we just need to call get
template <typename T>
void Field<T>::get_elements(T *elements,
                            const std::vector<CoordinateVector> &coord_list) const {
    for (int i = 0; i < coord_list.size(); i++) {
        elements[i] = (*this)[coord_list[i]];
    }
}
#endif


#ifdef HILAPP

////////////////////////////////////////////////////////////////////////////////
// A couple of placeholder functions, not included in produced code.
// These are here in order for hilapp to generate explicitly
// some Direction and CoordinateVector operations, which may not exist in
// original code as such.  It is easiest to let the general hilapp
// code generation to do it using this hack, instead of hard-coding these to hilapp.
//
// These are needed because hilapp changes X+d-d -> +d-d, which may involve an
// operator not met before

inline void dummy_X_f() {
    Direction d1 = e_x;
    CoordinateVector v1(0);
    onsites(ALL) {
        Direction d;
        d = +d1;
        d = -d1; // unaryops
        CoordinateVector vec;
        vec = +v1;
        vec = -v1;

        // Direction + vector combos
        vec = d + d1;
        vec = d - d1;
        vec = v1 + d1;
        vec = v1 - d1;
        vec = d1 + v1;
        vec = d1 - v1;
        vec = vec + v1;
        vec = vec - v1;

        // and Direction index func
        vec[e_x] = vec[0] + vec[e_y] + vec.e(e_y);
    }
}

/// Dummy function including Field<T> functions and methods which
/// need to be explicitly seen by hilapp during 1st pass in order to generater
/// necessary functions.  Add here ops as needed

template <typename T>
inline void ensure_field_operators_exist(Field<T> &f) {

    onsites(ALL) {
        f[X] = 0;     // set to zero
        f[X] = -f[X]; // unary -  -- needed for antiperiodic b.c.
    }
    // same for non-vectorized loop
    onsites(ALL) {
        if (X.coordinate(e_x) < X.coordinate(e_y)) {
            f[X] = 0;
            f[X] = -f[X];
        }
    }

    // make shift also explicit
    CoordinateVector v = 0;
    f = f.shift(v, ALL);
}

#endif

#endif // FIELD_H
