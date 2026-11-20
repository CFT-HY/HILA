#ifndef HILA_IO_FLAGS_H_
#define HILA_IO_FLAGS_H_

#include <string>
#include <type_traits>

namespace hila {
namespace io {

/**
 * @brief Flags for Field output, combine with |, e.g. hila::io::hdf5 | hila::io::to_float
 * @details At most one format (ascii, binary, hdf5) may be given.  Files with extension .h5
 * or .hdf5 are always written in HDF5 format (also by the functions without flags, e.g.
 * f.write("phi.h5")); an ascii or binary flag with such a name is an error.  Otherwise the
 * default is binary.
 */
enum output_flags : unsigned {
    none = 0,     ///< no flags: format from file name
    ascii = 1,    ///< text output
    binary = 2,   ///< raw binary output
    hdf5 = 4,     ///< HDF5 file
    to_float = 8, ///< write floating point numbers in single precision (binary, hdf5)
};

constexpr output_flags operator|(output_flags a, output_flags b) {
    return output_flags((unsigned)a | (unsigned)b);
}

/// File names ending with .h5 or .hdf5 are written in HDF5 format, unless a format flag
/// says otherwise
inline bool is_hdf5_filename(const std::string &filename) {
    auto ends_with = [&](const std::string &ext) {
        return filename.size() >= ext.size() &&
               filename.compare(filename.size() - ext.size(), ext.size(), ext) == 0;
    };
    return ends_with(".h5") || ends_with(".hdf5");
}

/// @internal File name of an output target (std::string or char array), "" for streams
template <typename Outf>
std::string output_name(const Outf &output) {
    if constexpr (std::is_convertible<const Outf &, std::string>::value)
        return std::string(output);
    else
        return std::string();
}

/// @internal Output format (ascii, binary or hdf5) from flags and file name ("" for streams)
inline output_flags output_format(output_flags flags, const std::string &filename) {
    unsigned f = flags & (ascii | binary | hdf5);
    if (f == 0)
        return is_hdf5_filename(filename) ? hdf5 : binary;
    if (f != ascii && f != binary && f != hdf5)
        return none; // conflicting formats, caller reports
    if (f != hdf5 && is_hdf5_filename(filename))
        return none; // .h5 file name with ascii or binary flag
    return output_flags(f);
}

} // namespace io
} // namespace hila

#endif
