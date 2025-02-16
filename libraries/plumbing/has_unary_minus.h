#ifndef HAS_UNARY_MINUS_H
#define HAS_UNARY_MINUS_H

#include <type_traits>

/**
 *@details Namespace `hila` contains most of `class templates` and `function templates`, which are
 *  necessary to run lattice filed simulations on distributed memory system as well as
 *  on GPU nodes concurrently.
 */
namespace hila {

/**
 *@brief Conditionally reture `bool` type `false` if type `T` does't have unary `-` operator
 *
 *@details If the type T does not implement `-T` (unary `-`) operator, i.e.
 *\code{.cpp}
 *    T T::operator-() const { ... }
 *\endcode,
 *`has_unary_minus::value` is `false`. This is needed for antiperiodic boundary conditions
 * @note `value` is false for `unsigned` type, whereas c++ allows `-unsigned`
 */
template <typename T, typename A = void>
class has_unary_minus {
  public:
    static constexpr bool value = false;
};

template <typename T>
class has_unary_minus<
    T, typename std::enable_if_t<!std::is_unsigned<hila::arithmetic_type<T>>::value &&
                                 hila::is_assignable<T &, decltype(-std::declval<T>())>::value>> {
  public:
    static constexpr bool value = true;
};

/**
 *@brief hila::has_assign_zero<T>::value  returns true if '= 0 is defined for T
 *
 */
template <typename T, typename A = void>
class has_assign_zero {
  public:
    static constexpr bool value = false;
};

template <typename T>
class has_assign_zero<T, typename std::enable_if_t<hila::is_arithmetic<T>::value ||
                             hila::is_assignable<T &, int>::value ||
                             hila::is_assignable<T &, std::nullptr_t>::value>> {
  public:
    static constexpr bool value = true;
};


} // namespace hila

#endif
