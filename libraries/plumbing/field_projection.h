#ifndef HILA_FIELD_PROJECTION_H_
#define HILA_FIELD_PROJECTION_H_

//////////////////////////////////////////////////////////////////////
/// Projections of Fields: sums over one or more directions,
/// Field<T>::get_projection(), write_projection() and write_projection_binary().
///
/// Method:
///  1. each MPI rank sums its own sites into an array over the kept directions,
///     using a node-local ReductionVector (works on CPU and GPU).  The array covers
///     only the rank's own range of the kept coordinates.
///  2. ranks with the same position in the node grid in the kept directions
///     form a sub-communicator, and an MPI_Reduce collects the sum to its lowest rank.
///  3. rank 0 receives the pieces one slab at a time (slab = one node division in the
///     highest kept direction) and hands the slab to the output function.
/// No rank holds more than its own piece, except rank 0, which also holds one slab.
//////////////////////////////////////////////////////////////////////

#include "plumbing/field.h"
#include "plumbing/field_io.h"
#include "plumbing/reductionvector.h"

#include <algorithm>

namespace hila {

/// @internal
/// Sum the field over the sites of this rank inside the box [lo, hi].
/// Site with coordinates c goes to index sum_d (c[d]-lo[d])*stride[d]; stride is 0 for the
/// projected directions.
template <typename T>
std::vector<T> projection_local_sum(const Field<T> &f, const CoordinateVector &lo,
                                    const CoordinateVector &hi, const CoordinateVector &stride,
                                    int n_local) {

    ReductionVector<T> rv(n_local);
    rv.node_local();

    onsites (ALL) {
        CoordinateVector c = X.coordinates();
        bool inside = true;
        int idx = 0;
        foralldir (d) {
            if (c[d] < lo[d] || c[d] > hi[d])
                inside = false;
            idx += (c[d] - lo[d]) * stride[d];
        }
        if (inside)
            rv[idx] += f[X];
    }

    return std::vector<T>(rv.begin(), rv.end());
}

} // namespace hila


template <typename T>
template <typename Func>
void Field<T>::projection_internal(const std::vector<Direction> &dirs,
                                   const CoordinateVector &cmin, const CoordinateVector &cmax,
                                   Func &&output) const {

    assert_all_ranks_present();

    assert(fs != nullptr && "projection of unallocated Field");
    assert(fs->mylattice.ptr() == lattice.ptr() &&
           "projection: the Field must belong to the current lattice");

    const lattice_struct &lat = *lattice.ptr();

    // which directions are projected
    CoordinateVector is_proj = 0;
    assert(dirs.size() > 0 && "projection: no directions given");
    for (auto d : dirs) {
        assert(is_up_dir(d) && "projection: directions must be e_x, e_y, ...");
        assert(is_proj[d] == 0 && "projection: direction given twice");
        is_proj[d] = 1;
    }
    foralldir (d) {
        assert(cmin[d] >= 0 && cmax[d] >= cmin[d] && cmax[d] < lat.l_size[d] &&
               "projection: range error");
    }

    // highest kept direction (-1 if all directions are projected)
    int d_slow = -1;
    foralldir (d)
        if (!is_proj[d])
            d_slow = d;

    MPI_Datatype dtype = get_MPI_number_type<T>();
    assert(dtype != MPI_BYTE && "projection: unknown number type");
    constexpr int n_numbers = sizeof(T) / sizeof(hila::arithmetic_type<T>);

    // Range of this rank's sites inside [cmin,cmax].  The kept part of the range is
    // the same for all ranks in the same sub-communicator
    CoordinateVector lo, hi, stride;
    int64_t n_local = 1;
    foralldir (d) {
        lo[d] = std::max(cmin[d], lat.mynode.min[d]);
        hi[d] = std::min(cmax[d], lat.mynode.min[d] + lat.mynode.size[d] - 1);
        if (is_proj[d]) {
            stride[d] = 0;
        } else {
            stride[d] = n_local;
            n_local *= std::max(hi[d] - lo[d] + 1, 0);
        }
    }

    // sub-communicator: ranks at the same node grid position in the kept directions
    int color = 0, mult = 1;
    foralldir (d) {
        if (!is_proj[d]) {
            const auto &div = lat.nodes.divisors[d];
            int nc = std::upper_bound(div.begin(), div.end(), (unsigned)lat.mynode.min[d]) -
                     div.begin() - 1;
            color += nc * mult;
            mult *= lat.nodes.n_divisions[d];
        }
    }
    const int myrank = hila::myrank();
    MPI_Comm subcomm;
    MPI_Comm_split(lat.mpi_comm_lat, color, myrank, &subcomm);
    int subrank;
    MPI_Comm_rank(subcomm, &subrank);

    // local sum and reduction within the sub-communicator.
    // Ranks with empty kept range have n_local == 0, and so do all the others in the same group
    std::vector<T> piece;
    if (n_local > 0) {
        piece = hila::projection_local_sum(*this, lo, hi, stride, n_local);

        if (subrank == 0)
            MPI_Reduce(MPI_IN_PLACE, piece.data(), n_local * n_numbers, dtype, MPI_SUM, 0,
                       subcomm);
        else
            MPI_Reduce(piece.data(), nullptr, n_local * n_numbers, dtype, MPI_SUM, 0, subcomm);
    }
    MPI_Comm_free(&subcomm);

    const bool is_root = (subrank == 0 && n_local > 0);

    // Rank 0 collects the list of pieces: (is_root, lo, hi) from each rank
    constexpr int hsize = 1 + 2 * NDIM;
    int header[hsize];
    header[0] = is_root;
    foralldir (d) {
        header[1 + d] = lo[d];
        header[1 + NDIM + d] = hi[d];
    }
    std::vector<int> headers;
    if (myrank == 0)
        headers.resize((size_t)hila::number_of_nodes() * hsize);

    MPI_Gather(header, hsize, MPI_INT, headers.data(), hsize, MPI_INT, 0, lat.mpi_comm_lat);

    constexpr int tag = 9137;

    if (myrank != 0) {
        if (is_root)
            MPI_Send(piece.data(), n_local * n_numbers, dtype, 0, tag, lat.mpi_comm_lat);
        return;
    }

    //////////////////////////////////////////////////////////////////////
    // rank 0 only from here on

    struct piece_info {
        int rank;
        CoordinateVector lo, hi;
    };
    std::vector<piece_info> pieces;
    for (int r = 0; r < hila::number_of_nodes(); r++) {
        const int *h = headers.data() + (size_t)r * hsize;
        if (h[0]) {
            piece_info p;
            p.rank = r;
            foralldir (d) {
                p.lo[d] = h[1 + d];
                p.hi[d] = h[1 + NDIM + d];
            }
            pieces.push_back(p);
        }
    }

    // order the pieces to slabs along the highest kept direction
    if (d_slow >= 0) {
        std::stable_sort(pieces.begin(), pieces.end(),
                         [d_slow](const piece_info &a, const piece_info &b) {
                             return a.lo[d_slow] < b.lo[d_slow];
                         });
    }

    // strides of the output array, highest kept direction separately
    CoordinateVector ostride = 0;
    size_t plane = 1;
    foralldir (d) {
        if (!is_proj[d] && d != d_slow) {
            ostride[d] = plane;
            plane *= cmax[d] - cmin[d] + 1;
        }
    }

    std::vector<T> slab, buf;
    size_t i = 0;
    while (i < pieces.size()) {
        int s_lo = 0, s_hi = 0;
        if (d_slow >= 0) {
            s_lo = pieces[i].lo[d_slow];
            s_hi = pieces[i].hi[d_slow];
        }
        size_t j = i;
        while (j < pieces.size() && (d_slow < 0 || pieces[j].lo[d_slow] == s_lo))
            j++;

        slab.resize(plane * (s_hi - s_lo + 1));

        for (size_t k = i; k < j; k++) {
            const piece_info &p = pieces[k];
            size_t np = 1;
            foralldir (d)
                if (!is_proj[d])
                    np *= p.hi[d] - p.lo[d] + 1;

            const T *src;
            if (p.rank == 0) {
                src = piece.data();
            } else {
                buf.resize(np);
                MPI_Recv(buf.data(), np * n_numbers, dtype, p.rank, tag, lat.mpi_comm_lat,
                         MPI_STATUS_IGNORE);
                src = buf.data();
            }

            // copy the piece to the slab, piece is in natural order over the kept dirs
            CoordinateVector c = p.lo;
            for (size_t e = 0; e < np; e++) {
                size_t idx = 0;
                foralldir (d) {
                    if (!is_proj[d]) {
                        if (d == d_slow)
                            idx += (c[d] - s_lo) * plane;
                        else
                            idx += (c[d] - cmin[d]) * ostride[d];
                    }
                }
                slab[idx] = src[e];

                foralldir (d) {
                    if (!is_proj[d]) {
                        if (++c[d] <= p.hi[d])
                            break;
                        c[d] = p.lo[d];
                    }
                }
            }
        }

        output(slab.data(), slab.size());
        i = j;
    }
}


template <typename T>
std::vector<T> Field<T>::get_projection(const std::vector<Direction> &dirs,
                                        const CoordinateVector &cmin, const CoordinateVector &cmax,
                                        bool broadcast) const {
    std::vector<T> res;
    projection_internal(dirs, cmin, cmax,
                        [&res](const T *data, size_t n) { res.insert(res.end(), data, data + n); });
    if (broadcast)
        hila::broadcast(res);
    return res;
}

template <typename T>
std::vector<T> Field<T>::get_projection(const std::vector<Direction> &dirs,
                                        bool broadcast) const {
    CoordinateVector cmin = 0, cmax = lattice.size();
    cmax.asArray() -= 1;
    return get_projection(dirs, cmin, cmax, broadcast);
}


template <typename T>
void Field<T>::write_projection_internal(std::ofstream &out, const std::vector<Direction> &dirs,
                                         const CoordinateVector &cmin,
                                         const CoordinateVector &cmax, int precision, bool binary,
                                         bool to_float) const {
    if (hila::myrank() == 0 && !binary)
        out.precision(precision);

    projection_internal(dirs, cmin, cmax, [&](const T *data, size_t n) {
        hila::write_elements_to_stream(out, data, n, binary, to_float);
    });
}

template <typename T>
void Field<T>::write_projection_internal(const std::string &filename,
                                         const std::vector<Direction> &dirs,
                                         const CoordinateVector &cmin,
                                         const CoordinateVector &cmax, int precision, bool binary,
                                         bool to_float) const {
    std::ofstream out;
    hila::open_output_file(filename, out, binary);
    write_projection_internal(out, dirs, cmin, cmax, precision, binary, to_float);
    hila::close_file(filename, out);
}

template <typename T>
template <typename Outf>
void Field<T>::write_projection(Outf &&output, const std::vector<Direction> &dirs,
                                const CoordinateVector &cmin, const CoordinateVector &cmax,
                                int precision) const {
    write_projection_internal(output, dirs, cmin, cmax, precision, false, false);
}

template <typename T>
template <typename Outf>
void Field<T>::write_projection(Outf &&output, const std::vector<Direction> &dirs,
                                int precision) const {
    CoordinateVector cmin = 0, cmax = lattice.size();
    cmax.asArray() -= 1;
    write_projection_internal(output, dirs, cmin, cmax, precision, false, false);
}

template <typename T>
template <typename Outf>
void Field<T>::write_projection_binary(Outf &&output, const std::vector<Direction> &dirs,
                                       const CoordinateVector &cmin, const CoordinateVector &cmax,
                                       bool to_float) const {
    write_projection_internal(output, dirs, cmin, cmax, 0, true, to_float);
}

template <typename T>
template <typename Outf>
void Field<T>::write_projection_binary(Outf &&output, const std::vector<Direction> &dirs,
                                       bool to_float) const {
    CoordinateVector cmin = 0, cmax = lattice.size();
    cmax.asArray() -= 1;
    write_projection_internal(output, dirs, cmin, cmax, 0, true, to_float);
}

#endif
