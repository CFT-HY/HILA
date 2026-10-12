#ifndef HILA_SITE_LAYOUT_H_
#define HILA_SITE_LAYOUT_H_

/**
 * @file site_layout.h
 * @brief On-the-fly computation of site index <-> coordinates and neighbour indices
 *
 * @details
 * This reproduces the EVEN_SITES_FIRST layout (without GPU_OVERLAP_COMM) arithmetically, so
 * that the per-site arrays lattice.neighb[NDIRS][volume] and mynode.coordinates[volume]
 * (and their GPU copies) are not needed.
 *
 * Layout recap:
 *   - "lex" index of a site on the node: lc[0] + size[0]*(lc[1] + size[1]*(lc[2] + ...)),
 *     where lc = (global coordinates) - mynode.min
 *   - site index = lex/2                for even sites (parity of global coordinates)
 *                = lex/2 + evensites    for odd sites
 *     This works because in a lexicographically ordered box every aligned pair
 *     (lex 2m, 2m+1) contains one even and one odd site, for any box size and origin.
 *     Thus lex/2 is the rank of the site among the sites of the same parity.
 *   - halo for direction d (only if the node is divided to direction d, otherwise
 *     the neighbour wraps around on the node):
 *       halo_offset[d] + (site odd ? halo_evensites[d] : 0) + j/2
 *     where j is the lex index of the site on the boundary face perpendicular to d
 *     (x fastest, the direction d skipped), and the parity is that of the site whose
 *     neighbour is fetched.  Same pairing argument as above.
 *     halo_offset[d] are allocated in the order of directions e_x, e_y, ..., e_x_down.
 *
 * All functions marked HILA_HOST_DEVICE use only plain int arrays, so that they can be
 * called from GPU kernels and host code.  The struct is trivially copyable and small, so
 * it can be placed in __constant__ memory or passed as a kernel argument.
 *
 * NOT handled (yet): GPU_OVERLAP_COMM layout, SUBNODE_LAYOUT (vector backend),
 * SPECIAL_BOUNDARY_CONDITIONS which need extra halo elements.
 */

#include "plumbing/defs.h"
#include "plumbing/coordinates.h"

#if defined(CUDA) || defined(HIP)
#define HILA_HOST_DEVICE __host__ __device__
#else
#define HILA_HOST_DEVICE
#endif

namespace hila {

/// Coordinates and derived quantities of a site on this node, computed once per site
struct layout_site {
    int lc[NDIM];    // local coordinates, 0 <= lc[d] < size[d]
    unsigned lex;    // lexicographic index on the node
    unsigned parity; // 0 = even, 1 = odd (parity of the global coordinates)
};

struct site_layout_struct {

    int min[NDIM];  // node origin (global coordinates)
    int size[NDIM]; // node size
    unsigned size_factor[NDIM]; // 1, size[0], size[0]*size[1], ...
    unsigned volume;
    unsigned evensites;
    unsigned min_parity; // parity of the node origin
    bool wrap[NDIM];     // true if node is not divided to direction d - neighbours wrap on node

    // halo info for each direction
    unsigned halo_offset[NDIRS];
    unsigned halo_evensites[NDIRS];
    // face lex factors: face_factor[k][e] for face perpendicular to dim k, face_factor[k][k] = 0
    unsigned face_factor[NDIM][NDIM];

    unsigned halo_end;         // volume + all nn halos
    unsigned field_alloc_size; // halo_end rounded up to a multiple of 256, as in lattice setup
                               // (special boundary elements, if any, go after halo_end)

    //////////////////////////////////////////////////////////////////////////////////
    /// Host setup.  ndiv[d] is the number of node divisions to direction d
    //////////////////////////////////////////////////////////////////////////////////

    void setup(const CoordinateVector &nmin, const CoordinateVector &nsize,
               const CoordinateVector &ndiv) {

        unsigned v = 1;
        unsigned ps = 0;
        foralldir (d) {
            min[d] = nmin[d];
            size[d] = nsize[d];
            size_factor[d] = v;
            v *= size[d];
            ps += min[d];
            wrap[d] = (ndiv[d] == 1);
        }
        volume = v;
        min_parity = ps & 1;
        // pairs (2m,2m+1) contain one of each parity; if volume is odd, the last site
        // (lc = size-1) is unpaired
        evensites = volume / 2;
        if (volume % 2 == 1) {
            int last[NDIM];
            foralldir (d)
                last[d] = size[d] - 1;
            if (parity_of(last) == 0)
                evensites++;
        }

        foralldir (k) {
            unsigned f = 1;
            foralldir (e) {
                if (e == k) {
                    face_factor[k][e] = 0;
                } else {
                    face_factor[k][e] = f;
                    f *= size[e];
                }
            }
        }

        // halos in the order of directions, as in create_std_gathers()
        size_t offset = volume;
        for (int d = 0; d < NDIRS; d++) {
            int k = dim_of(d);
            if (wrap[k]) {
                halo_offset[d] = 0;
                halo_evensites[d] = 0;
            } else {
                // count even sites on the face
                unsigned face = volume / size[k];
                unsigned even = 0;
                int lc[NDIM];
                for (unsigned j = 0; j < face; j++) {
                    unsigned r = j;
                    foralldir (e) {
                        if (e == k) {
                            lc[e] = is_up(d) ? size[e] - 1 : 0;
                        } else {
                            lc[e] = r % size[e];
                            r /= size[e];
                        }
                    }
                    if (parity_of(lc) == 0)
                        even++;
                }
                halo_offset[d] = offset;
                halo_evensites[d] = even;
                offset += face;
            }
        }
        assert(offset < (1ULL << 32) && "site_layout: too large node for 32-bit indices");
        halo_end = offset;
        field_alloc_size = (halo_end + 255) / 256 * 256;
    }

    //////////////////////////////////////////////////////////////////////////////////
    /// Index arithmetic, usable on host and device
    //////////////////////////////////////////////////////////////////////////////////

    static HILA_HOST_DEVICE inline int dim_of(int d) {
        return d < NDIM ? d : NDIRS - 1 - d;
    }
    static HILA_HOST_DEVICE inline bool is_up(int d) {
        return d < NDIM;
    }

    HILA_HOST_DEVICE inline unsigned parity_of(const int (&lc)[NDIM]) const {
        unsigned s = min_parity;
        for (int d = 0; d < NDIM; d++)
            s += lc[d];
        return s & 1;
    }

    HILA_HOST_DEVICE inline unsigned lex_of(const int (&lc)[NDIM]) const {
        unsigned lex = 0;
        for (int d = 0; d < NDIM; d++)
            lex += lc[d] * size_factor[d];
        return lex;
    }

    /// site index from lex index and parity
    HILA_HOST_DEVICE inline unsigned index_of(unsigned lex, unsigned parity) const {
        return (lex >> 1) + (parity ? evensites : 0);
    }

    /// site index from local coordinates
    HILA_HOST_DEVICE inline unsigned index_of(const int (&lc)[NDIM]) const {
        return index_of(lex_of(lc), parity_of(lc));
    }

    /// local coordinates from lex index
    HILA_HOST_DEVICE inline void lex_to_local(unsigned lex, int (&lc)[NDIM]) const {
        for (int d = 0; d < NDIM - 1; d++) {
            unsigned q = lex / size[d];
            lc[d] = lex - q * size[d];
            lex = q;
        }
        lc[NDIM - 1] = lex;
    }

    /// Everything about site with index idx
    HILA_HOST_DEVICE inline layout_site site(unsigned idx) const {
        layout_site s;
        s.parity = (idx >= evensites) ? 1 : 0;
        unsigned lex = (idx - (s.parity ? evensites : 0)) << 1;
        lex_to_local(lex, s.lc);
        if (parity_of(s.lc) != s.parity) {
            // it is the other member of the pair, lex+1
            lex++;
            for (int d = 0; d < NDIM; d++) {
                if (++s.lc[d] < size[d])
                    break;
                s.lc[d] = 0;
            }
        }
        s.lex = lex;
        return s;
    }

    /// Index of the neighbour of site s to direction d in the field array (incl. halo)
    HILA_HOST_DEVICE inline unsigned neighbour(const layout_site &s, int d) const {
        int k = dim_of(d);
        if (is_up(d)) {
            if (s.lc[k] < size[k] - 1)
                return index_of(s.lex + size_factor[k], s.parity ^ 1);
            if (wrap[k])
                return index_of(s.lex - (size[k] - 1) * size_factor[k],
                                s.parity ^ ((size[k] - 1) & 1));
        } else {
            if (s.lc[k] > 0)
                return index_of(s.lex - size_factor[k], s.parity ^ 1);
            if (wrap[k])
                return index_of(s.lex + (size[k] - 1) * size_factor[k],
                                s.parity ^ ((size[k] - 1) & 1));
        }
        // now in halo
        return halo_index(s, d);
    }

    /// Halo index of the off-node neighbour of boundary site s to direction d
    HILA_HOST_DEVICE inline unsigned halo_index(const layout_site &s, int d) const {
        int k = dim_of(d);
        unsigned j = 0;
        for (int e = 0; e < NDIM; e++)
            j += s.lc[e] * face_factor[k][e];
        return halo_offset[d] + (s.parity ? halo_evensites[d] : 0) + (j >> 1);
    }

    /// true if the neighbour of s to direction d is off-node
    HILA_HOST_DEVICE inline bool is_off_node(const layout_site &s, int d) const {
        int k = dim_of(d);
        if (wrap[k])
            return false;
        return is_up(d) ? (s.lc[k] == size[k] - 1) : (s.lc[k] == 0);
    }

    HILA_HOST_DEVICE inline unsigned neighbour(unsigned idx, int d) const {
        return neighbour(site(idx), d);
    }

    //////////////////////////////////////////////////////////////////////////////////
    /// CoordinateVector interface (host, and device via hilapp loop functions)
    //////////////////////////////////////////////////////////////////////////////////

#pragma hila loop_function
    inline CoordinateVector coordinates(unsigned idx) const {
        layout_site s = site(idx);
        CoordinateVector c;
        foralldir (d)
            c.e(d) = s.lc[d] + min[d];
        return c;
    }

    /// site index of on-node global coordinate cv
#pragma hila loop_function
    inline unsigned site_index(const CoordinateVector &cv) const {
        int lc[NDIM];
        foralldir (d)
            lc[d] = cv[d] - min[d];
        return index_of(lc);
    }

    inline ::Parity site_parity(unsigned idx) const {
        return idx < evensites ? EVEN : ODD;
    }
};

static_assert(std::is_trivially_copyable<site_layout_struct>::value,
              "site_layout_struct must be trivially copyable (GPU __constant__ memory)");

} // namespace hila

#endif
