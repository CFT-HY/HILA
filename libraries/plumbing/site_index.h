#ifndef SITEINDEX_H
#define SITEINDEX_H

#include "coordinates.h"
#include "lattice.h"


/**
 * @brief Running index for locating sites on the lattice
 *
 * The index is computed as x + nx*(y + ny*(z + nz*t) (in 4d)
 * Used sometimes instead of CoordinateVector, smaller storage.
 * Can be directly accessed by .value
 *
 */

class SiteIndex {
  public:
    size_t value;

    // std incantation for field types
    using base_type = size_t;
    using argument_type = size_t;

    SiteIndex() = default;
    SiteIndex(const SiteIndex &s) = default;
    SiteIndex(size_t v) : value(v) {}

    ~SiteIndex() = default;

    /**
     * @brief Construct from CoordinateVector
     *
     * Assumes cv is a valid lattice coordinate, undefined otherwise
     */
    SiteIndex(const CoordinateVector &cv) {
        value = 0;
        size_t m = 1;
        foralldir(d) {
            value += m * cv[d];
            m *= lattice.size(d);
        }
    }

    /**
     * @brief Convert to lattice coordinates
     */
    CoordinateVector coordinates() const {
        CoordinateVector res;
        size_t v = value;
        foralldir(d) {
            res.e(d) = v % lattice.size(d);
            v /= lattice.size(d);
        }
        return res;
    }
};

#endif
