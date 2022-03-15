#ifndef HAS_UNARY_MINUS_H
#define HAS_UNARY_MINUS_H

#include <type_traits>

/////////////////////////////////////////////////////////////////////////////////////////////////
/// This header implements "has_unary_minus<T>::value" -conditional, which is false
/// if the type T does not implement -T (unary -) operator, i.e.
///    T T::operator-() const { .. }
/// This is needed for antiperiodic b.c.
/// slightly modified from Valentin Milea's example in
/// https://stackoverflow.com/a/31539364

template <class C, class A = void> class has_unary_minus {
    template <class T> static std::true_type testSignature(T (T::*)(void) const);

    template <class T> static decltype(testSignature(&T::operator-)) test(std::nullptr_t);

    template <class T> static std::false_type test(...);

  public:
    using type = decltype(test<C>(nullptr));
    static const bool value = type::value;
};

/// specialize to elementary arithmetic types

template <typename T>
class has_unary_minus<T, typename std::enable_if_t<hila::is_arithmetic<T>::value>> {
  public:
    static const bool value = true;
};

#endif