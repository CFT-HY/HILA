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
///     static constexpr hila::bc kind(Direction d);    // b.c. to direction d (up or down)
///     static std::string name();                      // for printing
///
/// Non-periodic b.c. requires that SPECIAL_BOUNDARY_CONDITIONS is defined.
///////////////////////////////////////////////////////////////////////////////

#include <string>
#include <type_traits>
#include "plumbing/coordinates.h"

namespace hila {

/// list of boundary condition kinds per direction
enum class bc { PERIODIC, ANTIPERIODIC, DIRICHLET };

/// False if we have b.c. which does not require communication
inline bool bc_need_communication(hila::bc bc) {
    if (bc == hila::bc::DIRICHLET) {
        return false;
    } else {
        return true;
    }
}

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
