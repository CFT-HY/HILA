#ifndef HILA_BOUNDARY_CONDITION_H_
#define HILA_BOUNDARY_CONDITION_H_

///////////////////////////////////////////////////////////////////////////////
/// Boundary condition types for Field<T, BC>
///
/// The boundary condition (b.c.) of a Field is a compile-time property of its type:
///
///     Field<double> a;                                   // periodic
///     Field<double, hila::antiperiodic<e_t>> b;          // antiperiodic in t
///     Field<Complex<double>, hila::antiperiodic<e_x, e_y>> c;
///
/// If the BC argument is omitted, Field<T> uses the default b.c. of the element type T:
/// T::default_boundary_condition if T declares it, otherwise hila::periodic:
///
///     struct my_fermion {
///         ...
///         using default_boundary_condition = hila::antiperiodic<e_t>;
///     };
///     Field<my_fermion> psi;                             // antiperiodic in t
///
/// For types which cannot be modified, specialize the trait:
///
///     template <>
///     struct hila::default_boundary_condition<Wilson_vector<3, double>> {
///         using type = hila::antiperiodic<e_t>;
///     };
///
/// A b.c. type is an empty struct with static members:
///     static constexpr bool is_periodic;              // true only if periodic to all directions
///     static constexpr unsigned antiperiodic_mask;    // bit dim set if antiperiodic to dim
///     static constexpr hila::bc kind(Direction d);    // b.c. to direction d (up or down)
///     static std::string name();                      // for printing
///
/// The b.c. is applied when a neighbour value is read: the value of f[X + d] is negated if
/// the step from X to direction d crosses an antiperiodic lattice boundary, see
/// hila::bc_nb_value().  The field storage and the halo contain the plain (periodic) values,
/// thus the communications and the field memory layout do not depend on the b.c.
///////////////////////////////////////////////////////////////////////////////

#include <string>
#include <type_traits>
#include "plumbing/coordinates.h"

#ifndef HILA_HOST_DEVICE
#if defined(CUDA) || defined(HIP)
#define HILA_HOST_DEVICE __host__ __device__
#else
#define HILA_HOST_DEVICE
#endif
#endif

namespace hila {

/// list of boundary condition kinds per direction
enum class bc { PERIODIC, ANTIPERIODIC };

namespace bc_detail {
/// constexpr "abs" of Direction: e_x_down -> e_x etc.
constexpr unsigned up_dir(Direction d) {
    return (unsigned)d < NDIM ? (unsigned)d : NDIRS - 1 - (unsigned)d;
}
/// for static_assert in templates which must not be instantiated
template <typename...>
inline constexpr bool always_false = false;
} // namespace bc_detail

/// Periodic boundary conditions to all directions (the default)
struct periodic {
    static constexpr bool is_periodic = true;
    static constexpr unsigned antiperiodic_mask = 0;
    static constexpr hila::bc kind(Direction) {
        return hila::bc::PERIODIC;
    }
    static std::string name() {
        return "periodic";
    }
};

/// Antiperiodic boundary conditions to the listed directions, periodic to others.
/// Example: hila::antiperiodic<e_t>, hila::antiperiodic<e_x, e_y>
template <Direction... dirs>
struct antiperiodic {
    static_assert(sizeof...(dirs) > 0, "hila::antiperiodic<> needs at least one direction");
    static_assert(((dirs < NDIM) && ...),
                  "hila::antiperiodic<...>: give the up directions e_x, e_y, ...");

    static constexpr bool is_periodic = false;
    static constexpr unsigned antiperiodic_mask = ((1u << (unsigned)dirs) | ...);
    static constexpr hila::bc kind(Direction d) {
        return ((bc_detail::up_dir(d) == (unsigned)dirs) || ...) ? hila::bc::ANTIPERIODIC
                                                                 : hila::bc::PERIODIC;
    }
    static std::string name() {
        std::string s = "antiperiodic(";
        bool first = true;
        for (unsigned d : {(unsigned)dirs...}) {
            if (!first)
                s += ',';
            s += "xyzt"[d];
            first = false;
        }
        return s + ')';
    }
};

/// Default b.c. of element type T: T::default_boundary_condition if it exists, else periodic.
/// Can be specialized for types which cannot be modified.
template <typename T, typename = void>
struct default_boundary_condition {
    using type = hila::periodic;
};

template <typename T>
struct default_boundary_condition<T, std::void_t<typename T::default_boundary_condition>> {
    using type = typename T::default_boundary_condition;
};

template <typename T>
using default_boundary_condition_t = typename default_boundary_condition<T>::type;

/// dimension 0 .. NDIM-1 of direction d (up or down)
HILA_HOST_DEVICE constexpr int bc_dim(int d) {
    return d < NDIM ? d : NDIRS - 1 - d;
}

/// The value v of the neighbour of a site to direction d with the b.c. BC applied.
/// cross: the step crosses the lattice boundary.  Antiperiodic: v is negated if cross.
/// Used in the hilapp-generated neighbour accesses of non-periodic fields
template <typename BC, typename T>
HILA_HOST_DEVICE inline T bc_nb_value_cross(const T &v, int d, bool cross) {
    if constexpr (BC::is_periodic) {
        return v;
    } else {
        if (((BC::antiperiodic_mask >> bc_dim(d)) & 1u) && cross)
            return -v;
        return v;
    }
}

/// The same for the neighbour X + k*d (k >= 1): x is the global coordinate of X along the
/// dimension of d, L the lattice size along it
template <typename BC, typename T>
HILA_HOST_DEVICE inline T bc_nb_value(const T &v, int d, int k, int x, int L) {
    return bc_nb_value_cross<BC>(v, d, d < NDIM ? (x + k >= L) : (x < k));
}

/// Boundary condition of the result of a binary Field operation (e.g. Field + Field):
/// equal b.c.s are kept; if one of them is periodic, the other is used.
/// Two different non-periodic b.c.s are an error.
template <typename B1, typename B2>
struct combine_boundary_conditions {
    static_assert(B1::is_periodic || B2::is_periodic,
                  "Field operation with two different non-periodic boundary conditions");
    using type = std::conditional_t<B1::is_periodic, B2, B1>;
};

template <typename B>
struct combine_boundary_conditions<B, B> {
    using type = B;
};

template <typename B1, typename B2>
using combine_boundary_conditions_t = typename combine_boundary_conditions<B1, B2>::type;

} // namespace hila

#endif
