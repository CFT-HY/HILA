#ifndef SITEINDEX_H
#define SITEINDEX_H


//////////////////////////////////////////////////////////////////////
/// SiteIndex type - indexes all sites on the lattice, so that the
/// first dimension runs fastest.  Equivalent to CoordinateVector
//////////////////////////////////////////////////////////////////////

class SiteIndex {
  public:
    uint64_t value;

    // std incantation for field types
    using base_type = uint64_t;
    using argument_type = uint64_t;

    SiteIndex() = default;
    SiteIndex(const SiteIndex &s) = default;
    SiteIndex(uint64_t v) : value(v) {}
    ~SiteIndex() = default;

    SiteIndex(const CoordinateVector &cv) {
        value = 0;
        uint64_t m = 1;
        foralldir (d) {
            value += m * cv[d];
            m *= lattice.size(d);
       }
    }

    CoordinateVector coordinates() const {
        CoordinateVector res;
        uint64_t v = value;
        foralldir (d) {
            res.e(d) = v % lattice.size(d);
            v /= lattice.size(d);
        }
        return res;
    }
};

#endif
