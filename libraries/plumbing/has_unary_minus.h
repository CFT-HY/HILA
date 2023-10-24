#ifndef HAS_UNARY_MINUS_H
#define HAS_UNARY_MINUS_H

#include <type_traits>

/////////////////////////////////////////////////////////////////////////////////////////////////
/// This file implements "hila::has_unary_minus<T>::value" -conditional, which is false
/// if the type T does not implement -T (unary -) operator, i.e.
///    T T::operator-() const { .. }
/// This is needed for antiperiodic b.c.
/// Note that this gives false for unsigned type, whereas c++ allows -unsigned


namespace hila {

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

} // namespace hila

#endif