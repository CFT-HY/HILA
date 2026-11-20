#ifndef HILA_FIELD_OUTPUT_H_
#define HILA_FIELD_OUTPUT_H_

//////////////////////////////////////////////////////////////////////////////////
/// Field output functions with format flags (hila::io::output_flags):
///   f.write(filename, flags), f.write_subvolume(out, cmin, cmax, flags),
///   f.write_slice(out, slice, flags), f.write_projection(out, dirs, [cmin, cmax,] flags)
/// These dispatch to the ascii, binary (parallel MPI-IO if hila::parallel_io) and HDF5
/// (hila::hdf5_file) writers.
//////////////////////////////////////////////////////////////////////////////////

#include "plumbing/field_io.h"
#include "plumbing/field_projection.h"
#include "plumbing/hdf5_file.h"

namespace hila {
namespace io {

/// @internal Output format (ascii, binary or hdf5) from flags, with error checks
inline output_flags checked_format(output_flags flags, const std::string &filename,
                                   bool is_file) {
    output_flags f = output_format(flags, filename);
    if (f == none)
        hila::error("Field output: conflicting format flags for " +
                    (is_file ? filename : std::string("stream")) +
                    " (.h5/.hdf5 file names are always HDF5)");
    if (!is_file && f == hdf5)
        hila::error("Field output: HDF5 output needs a file name, not a stream");
    return f;
}

/// @internal to_float applies only to floating point numbers wider than float
template <typename T>
bool effective_to_float(output_flags flags) {
    using A = hila::arithmetic_type<T>;
    return (flags & to_float) && std::is_floating_point<A>::value && sizeof(A) > sizeof(float);
}

} // namespace io
} // namespace hila


template <typename T, typename BC>
void Field<T, BC>::write(const std::string &filename, hila::io::output_flags flags,
                         int precision) const {
    auto fmt = hila::io::checked_format(flags, filename, true);

    if (fmt == hila::io::hdf5) {
        hila::hdf5_file h5(filename);
        h5.write("data", *this, flags);
    } else if (fmt == hila::io::ascii) {
        write(filename, false, precision);
    } else if (!hila::io::effective_to_float<T>(flags)) {
        write(filename, true);
    } else {
        // single precision binary: subvolume covering the whole lattice
        CoordinateVector cmin, cmax;
        foralldir (d) {
            cmin[d] = 0;
            cmax[d] = fs->mylattice.size(d) - 1;
        }
        write_subvolume(filename, cmin, cmax, flags);
    }
}

template <typename T, typename BC>
void Field<T, BC>::write(std::ofstream &outputfile, hila::io::output_flags flags,
                         int precision) const {
    auto fmt = hila::io::checked_format(flags, "", false);

    if (fmt == hila::io::ascii) {
        write(outputfile, false, precision);
    } else if (!hila::io::effective_to_float<T>(flags)) {
        write(outputfile, true);
    } else {
        CoordinateVector cmin, cmax;
        foralldir (d) {
            cmin[d] = 0;
            cmax[d] = fs->mylattice.size(d) - 1;
        }
        write_subvolume_internal(outputfile, cmin, cmax, 0, true, true);
    }
}

template <typename T, typename BC>
template <typename Outf>
void Field<T, BC>::write_subvolume(Outf &&output, const CoordinateVector &cmin,
                                   const CoordinateVector &cmax, hila::io::output_flags flags,
                                   int precision) const {

    constexpr bool is_file = std::is_convertible<Outf, std::string>::value;
    const std::string fname = hila::io::output_name(output);
    auto fmt = hila::io::checked_format(flags, fname, is_file);
    const bool tof = hila::io::effective_to_float<T>(flags);

    if constexpr (is_file) {
        if (fmt == hila::io::hdf5) {
            hila::hdf5_file h5(fname);
            h5.write_subvolume("data", *this, cmin, cmax, flags);
            return;
        }
        if (fmt == hila::io::binary && hila::parallel_io) {
            MPI_File fh = hila::io::open_mpi_file(fname, true);
            hila::io::write_box_mpi(fh, 0, *this, cmin, cmax, tof, fname);
            hila::io::close_mpi_file(fh, fname);
            return;
        }
    }
    write_subvolume_internal(output, cmin, cmax, precision, fmt == hila::io::binary, tof);
}

template <typename T, typename BC>
template <typename Outf>
void Field<T, BC>::write_slice(Outf &&output, const CoordinateVector &slice,
                               hila::io::output_flags flags, int precision) const {

    constexpr bool is_file = std::is_convertible<Outf, std::string>::value;
    const std::string fname = hila::io::output_name(output);
    auto fmt = hila::io::checked_format(flags, fname, is_file);

    if constexpr (is_file) {
        if (fmt == hila::io::hdf5) {
            hila::hdf5_file h5(fname);
            h5.write_slice("data", *this, slice, flags);
            return;
        }
    }
    auto [cmin, cmax] = get_range_from_slice(slice);
    write_subvolume(output, cmin, cmax, flags, precision);
}

template <typename T, typename BC>
template <typename Outf>
void Field<T, BC>::write_projection(Outf &&output, const std::vector<Direction> &dirs,
                                    const CoordinateVector &cmin, const CoordinateVector &cmax,
                                    hila::io::output_flags flags, int precision) const {

    constexpr bool is_file = std::is_convertible<Outf, std::string>::value;
    const std::string fname = hila::io::output_name(output);
    auto fmt = hila::io::checked_format(flags, fname, is_file);

    if constexpr (is_file) {
        if (fmt == hila::io::hdf5) {
            hila::hdf5_file h5(fname);
            h5.write_projection("data", *this, dirs, cmin, cmax, flags);
            return;
        }
    }
    write_projection_internal(output, dirs, cmin, cmax, precision, fmt == hila::io::binary,
                              hila::io::effective_to_float<T>(flags));
}

template <typename T, typename BC>
template <typename Outf>
void Field<T, BC>::write_projection(Outf &&output, const std::vector<Direction> &dirs,
                                    hila::io::output_flags flags, int precision) const {
    CoordinateVector cmin, cmax;
    foralldir (d) {
        cmin[d] = 0;
        cmax[d] = lattice.size(d) - 1;
    }
    write_projection(output, dirs, cmin, cmax, flags, precision);
}

#endif
