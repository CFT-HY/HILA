#ifndef SITEINDEX_H
#define SITEINDEX_H


//////////////////////////////////////////////////////////////////////
/// SiteIndex type - indexes all sites on the lattice, so that the
/// first dimension runs fastest.  Equivalent to CoordinateVector
//////////////////////////////////////////////////////////////////////

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

    SiteIndex(const CoordinateVector &cv) {
        value = 0;
        size_t m = 1;
        foralldir (d) {
            value += m * cv[d];
            m *= lattice.size(d);
       }
    }

    CoordinateVector coordinates() const {
        CoordinateVector res;
        size_t v = value;
        foralldir (d) {
            res.e(d) = v % lattice.size(d);
            v /= lattice.size(d);
        }
        return res;
    }
};

#endif
