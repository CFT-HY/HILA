#ifndef VECTORFIELD_H_
#define VECTORFIELD_H_

#include "hila.h"

template <typename T>
class VectorField {
private:
    std::array<Field<T>,NDIM> fdir;

public:
    // Default constructor
    VectorField() = default;

    // Straightforward copy constructor seems to be necessary
    VectorField(const VectorField &other) = default;

    // copy constructor - from fields which can be assigned
    template <typename A, std::enable_if_t<std::is_convertible<A, T>::value, int> = 0>
    VectorField(const VectorField<A> &other) {
        foralldir(d) fdir[d] = other[d];
    }

    // constructor with compatible scalar
    template <typename A, std::enable_if_t<hila::is_assignable<T &, A>::value, int> = 0>
    VectorField(const A &val) {
        foralldir(d) fdir[d] = val;
    }

    // constructor from 0 - nullptr trick in use
    VectorField(const std::nullptr_t z) {
        foralldir(d) fdir[d] = 0;
    }

    // move constructor - steal the content
    VectorField(VectorField &&rhs) = default;

    /////////////////////////////////////////////////
    /// Destructor

    ~VectorField() = default;

    /////////////////////////////////////////////////
    /// Access components with []

    inline Field<T> & operator[](Direction d) {
        return fdir[d];
    }

    inline const Field<T> & operator[](Direction d) const {
        return fdir[d];
    }

    //////////////////////////////////////////////////
    /// Assign from anything the field allows
    template <typename A>
    VectorField & operator=(const A &val) {
        foralldir(d) fdir[d] = val;
        return *this;
    }

    /// Separate 0 assignment
    VectorField & operator=(std::nullptr_t np) {
        foralldir(d) fdir[d] = 0;
        return *this;
    }

    

};


///////////////////////////////////////////////////////
/// Alias gaugefield to VectorField
template <typename T>
using GaugeField = VectorField<T>;


#endif