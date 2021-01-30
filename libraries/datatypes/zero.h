#ifndef ZERO_H_
#define ZERO_H_

////////////////////////////////////////////////////////////////////
/// This class defines const "zero" and Zero() which can be assigned
/// to any built-in type and zeroes the value. 
///  A = 0 does not work e.g. for vectors or matrices, because A = 1 is ambiguous!
/// Here typecasts to basic types are implemented

#include <cstdint>

struct Zero {
  public:
    Zero() = default;
    operator int() const {return 0;}
    operator int64_t() const {return (int64_t)0;}
    operator float() const {return 0.0f; }
    operator double() const {return 0.0;}
    operator long double() const {return 0.0L;}
};

static constexpr Zero zero;

#endif