#ifndef ZERO_H_
#define ZERO_H_

////////////////////////////////////////////////////////////////////
/// This class defines const "zero" which can be assigned
/// to any built-in type and zeroes the value. 
///  A = 0 does not work e.g. for vectors or matrices, because A = 1 is ambiguous!
/// Here typecasts to basic types are implemented, 

#include <cstdint>

struct Zero {
  public:
    #pragma hila loop_function
    operator int() const {return 0;}
    #pragma hila loop_function
    operator int64_t() const {return (int64_t)0;}
    #pragma hila loop_function
    operator float() const {return 0.0f; }
    #pragma hila loop_function
    operator double() const {return 0.0;}
    #pragma hila loop_function
    operator long double() const {return 0.0L;}
};

static Zero zero;


#endif