#ifndef GAUGEFIELD_H_
#define GAUGEFIELD_H_

#include "hila.h"

template <typename T>
class GaugeField {
private:
    std::array<Field<T>,NDIM> fdir;

public:
    // Default constructor
    GaugeField() = default;

    // Straightforward copy constructor seems to be necessary
    GaugeField(const GaugeField &other) = default;

    // copy constructor - from fields which can be assigned
    template <typename A, std::enable_if_t<std::is_convertible<A, T>::value, int> = 0>
    GaugeField(const GaugeField<A> &other) {
        foralldir(d) fdir[d] = other[d];
    }

    // constructor with compatible scalar
    template <typename A, std::enable_if_t<hila::is_assignable<T &, A>::value, int> = 0>
    GaugeField(const A &val) {
        foralldir(d) fdir[d] = val;
    }

    // constructor from 0 - nullptr trick in use
    GaugeField(const std::nullptr_t z) {
        foralldir(d) fdir[d] = 0;
    }

    // move constructor - steal the content
    GaugeField(GaugeField &&rhs) = default;

    /////////////////////////////////////////////////
    /// Destructor

    ~GaugeField() = default;

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
    GaugeField & operator=(const A &val) {
        foralldir(d) fdir[d] = val;
        return *this;
    }

    /// Separate 0 assignment
    GaugeField & operator=(std::nullptr_t np) {
        foralldir(d) fdir[d] = 0;
        return *this;
    }

    

};


///////////////////////////////////////////////////////
/// Alias VectorField to GaugeField
template <typename T>
using VectorField = GaugeField<T>;


#endif