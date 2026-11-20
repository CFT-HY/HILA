#ifndef HILA_HDF5_FILE_H_
#define HILA_HDF5_FILE_H_

/**
 * @file hdf5_file.h
 * @brief Write HDF5 files without the HDF5 library
 *
 * hila::hdf5_file writes a small, self-contained subset of the HDF5 file format (format
 * specification 3.0, readable by HDF5 >= 1.8: h5py, h5dump, Julia, Matlab, ...):
 * - groups, with any number of datasets and subgroups
 * - datasets with contiguous storage of integer or floating point numbers (no chunking or
 *   compression)
 * - attributes (numbers, arrays of numbers, strings) on groups and datasets
 *
 * Field data is written with collective MPI-IO directly from all ranks.  The file is
 * write-once: the raw data of datasets is written as they come, and the metadata (object
 * headers and superblock) when the file is closed.  hila does not read HDF5 files.
 *
 * A Field<T> on a lattice of size (L_0, L_1, ..., L_{NDIM-1}) becomes a dataset of shape
 * [L_{NDIM-1}, ..., L_1, L_0, ncomp] (C order, last index fastest, i.e. x runs fastest as in
 * hila).  ncomp = number of real numbers in T (e.g. 2 for Complex<double>, 18 for SU<3,double>),
 * omitted if 1.  Complex numbers and matrices are stored in their memory order (re, im, ...).
 *
 * Example:
 * @code
 *   hila::hdf5_file h5("meas.h5");
 *   h5.attribute("beta", beta);
 *   h5.attribute("lattice size", lattice.size());
 *   h5.write("phi", phi).attribute("trajectory", traj);
 *   h5.write("gauge/U_x", U[e_x], hila::io::to_float);   // creates group "gauge"
 *   h5.close();                                           // or at destruction
 * @endcode
 * All ranks must call all methods.  Attribute values and std::vector data are taken from rank 0.
 */

#include "plumbing/field_io.h"
#include <cstdint>
#include <string>
#include <type_traits>
#include <vector>

namespace hila {

class hdf5_file {
  public:
    /// @internal number type stored in the file
    struct number_type {
        bool is_float;
        bool is_signed;
        uint8_t size; // in bytes
    };

    template <typename A>
    static number_type number_type_of() {
        static_assert(std::is_arithmetic<A>::value, "HDF5 output: not an arithmetic type");
        static_assert(!std::is_floating_point<A>::value || sizeof(A) == 4 || sizeof(A) == 8,
                      "HDF5 output: only float and double floating point types supported");
        return {std::is_floating_point<A>::value, std::is_signed<A>::value, (uint8_t)sizeof(A)};
    }

    /// Reference to a group or dataset in the file, used to add attributes
    class object_ref {
        hdf5_file *file;
        int obj;

      public:
        object_ref(hdf5_file *f, int i) : file(f), obj(i) {}

        /// Add attribute: a number, std::vector of numbers, string, or a hila type of numbers
        /// (Vector, CoordinateVector, Complex, ...).  Value from rank 0 is used.
        template <typename V>
        object_ref &attribute(const std::string &name, const V &value) {
            file->add_attribute(obj, name, value);
            return *this;
        }
    };

    /// Create (or replace) the file.  Collective.
    explicit hdf5_file(const std::string &filename);
    ~hdf5_file();

    hdf5_file(const hdf5_file &) = delete;
    hdf5_file &operator=(const hdf5_file &) = delete;

    /// Write the metadata and close the file.  Collective.  Called by the destructor if needed.
    void close();

    /// The root group "/"
    object_ref root() {
        return object_ref(this, 0);
    }

    /// Group with the given path ("a/b/c"), created with intermediate groups if needed
    object_ref group(const std::string &path) {
        return object_ref(this, find_or_create_group(path));
    }

    /// Attribute of the root group
    template <typename V>
    object_ref attribute(const std::string &name, const V &value) {
        return root().attribute(name, value);
    }

    /**
     * @brief Write a Field as a dataset at path ("name" or "group/.../name"). Collective.
     * @param flags hila::io::to_float: write floating point numbers in single precision
     */
    template <typename T, typename BC>
    object_ref write(const std::string &path, const Field<T, BC> &f, unsigned flags = 0);

    /**
     * @brief Write data held on rank 0 as a dataset.  Collective.
     * @param dims  array dimensions, dims[0] running fastest (as x in hila); the product of
     *              dims must equal data.size().  Must be the same on all ranks.
     * @param flags hila::io::to_float: write floating point numbers in single precision
     */
    template <typename T>
    object_ref write(const std::string &path, const std::vector<T> &data,
                     const std::vector<size_t> &dims, unsigned flags = 0);

  private:
    struct attribute_t {
        std::string name;
        std::vector<uint8_t> dtype, dspace, data;
    };

    struct object_t {
        bool is_group;
        std::string name; // link name in parent group
        std::vector<int> children;
        std::vector<uint8_t> dtype, dspace; // datasets
        uint64_t addr = 0, size = 0;        // datasets: raw data location
        std::vector<attribute_t> attrs;
    };

    std::string filename;
    MPI_File fh;
    bool is_open = false;
    uint64_t eof; // end of allocated file space
    std::vector<object_t> objects;

    uint64_t reserve(uint64_t bytes);
    int find_or_create_group(const std::string &path);
    int add_dataset(const std::string &path, number_type nt, const std::vector<uint64_t> &dims,
                    uint64_t addr, uint64_t bytes);
    void write_rank0(uint64_t addr, const void *buf, uint64_t bytes);

    void add_attribute_raw(int obj, const std::string &name, number_type nt, const void *data,
                           const std::vector<uint64_t> &dims, bool scalar);

    void add_attribute(int obj, const std::string &name, const std::string &value);
    void add_attribute(int obj, const std::string &name, const char *value) {
        add_attribute(obj, name, std::string(value));
    }

    template <typename V, std::enable_if_t<std::is_arithmetic<V>::value, int> = 0>
    void add_attribute(int obj, const std::string &name, const V &value) {
        add_attribute_raw(obj, name, number_type_of<V>(), &value, std::vector<uint64_t>(), true);
    }

    template <typename V, std::enable_if_t<std::is_arithmetic<V>::value, int> = 0>
    void add_attribute(int obj, const std::string &name, const std::vector<V> &value) {
        add_attribute_raw(obj, name, number_type_of<V>(), value.data(),
                          std::vector<uint64_t>(1, value.size()), false);
    }

    // hila types consisting of numbers: Vector, CoordinateVector, Complex, Matrix ...
    template <typename V, std::enable_if_t<std::is_class<V>::value, int> = 0>
    void add_attribute(int obj, const std::string &name, const V &value) {
        using A = hila::arithmetic_type<V>;
        static_assert(std::is_arithmetic<A>::value && sizeof(V) % sizeof(A) == 0,
                      "HDF5 attribute: unsupported type");
        add_attribute_raw(obj, name, number_type_of<A>(), &value,
                          std::vector<uint64_t>(1, sizeof(V) / sizeof(A)), false);
    }

    void serialize(int obj, std::vector<uint8_t> &meta, uint64_t meta_start,
                   std::vector<uint64_t> &header_addr) const;
};

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename BC>
hdf5_file::object_ref hdf5_file::write(const std::string &path, const Field<T, BC> &f,
                                       unsigned flags) {
    using A = hila::arithmetic_type<T>;
    static_assert(std::is_arithmetic<A>::value,
                  "HDF5 output: Field element must consist of arithmetic numbers");
    constexpr size_t ncomp = sizeof(T) / sizeof(A);
    static_assert(ncomp * sizeof(A) == sizeof(T), "HDF5 output: Field element has padding");

    if (!f.is_allocated())
        hila::error("hdf5_file::write: Field for dataset " + path + " is not allocated");

    Lattice flat = f.mylattice();
    const CoordinateVector gsize = flat.size();

    std::vector<uint64_t> dims;
    for (int d = NDIM - 1; d >= 0; d--)
        dims.push_back(gsize[d]);
    if (ncomp > 1)
        dims.push_back(ncomp);

    const bool tof =
        (flags & hila::io::to_float) && std::is_floating_point<A>::value && sizeof(A) > 4;
    const number_type nt = tof ? number_type_of<float>() : number_type_of<A>();

    const uint64_t bytes = (uint64_t)flat.volume() * ncomp * nt.size;
    const uint64_t addr = reserve(bytes);
    int obj = add_dataset(path, nt, dims, addr, bytes);

    if (!tof) {
        f.write_mpi(fh, addr, filename);
    } else {
        // copy_local_data() loops over the current lattice
        Lattice saved_lattice = lattice;
        lattice.switch_to(flat);

        std::vector<T> buf;
        f.copy_local_data(buf);
        const A *p = reinterpret_cast<const A *>(buf.data());
        std::vector<float> fbuf(buf.size() * ncomp);
        for (size_t i = 0; i < fbuf.size(); i++)
            fbuf[i] = (float)p[i];

        hila::io::write_block_mpi(fh, addr, fbuf.data(), ncomp * sizeof(float), gsize,
                                  lattice->mynode.size, lattice->mynode.min, filename);
        lattice.switch_to(saved_lattice);
    }
    return object_ref(this, obj);
}

template <typename T>
hdf5_file::object_ref hdf5_file::write(const std::string &path, const std::vector<T> &data,
                                       const std::vector<size_t> &dims, unsigned flags) {
    using A = std::conditional_t<std::is_arithmetic<T>::value, T, hila::arithmetic_type<T>>;
    static_assert(std::is_arithmetic<A>::value,
                  "HDF5 output: element must consist of arithmetic numbers");
    constexpr size_t ncomp = sizeof(T) / sizeof(A);
    static_assert(ncomp * sizeof(A) == sizeof(T), "HDF5 output: element has padding");

    size_t nelem = 1;
    std::vector<uint64_t> hdims;
    for (auto it = dims.rbegin(); it != dims.rend(); ++it) {
        hdims.push_back(*it);
        nelem *= *it;
    }
    if (ncomp > 1)
        hdims.push_back(ncomp);

    bool ok = true;
    if_rank0 () {
        ok = (data.size() == nelem);
    }
    if (!hila::broadcast(ok))
        hila::error("hdf5_file::write: data size does not match dims for dataset " + path);

    const bool tof =
        (flags & hila::io::to_float) && std::is_floating_point<A>::value && sizeof(A) > 4;
    const number_type nt = tof ? number_type_of<float>() : number_type_of<A>();

    const uint64_t bytes = (uint64_t)nelem * ncomp * nt.size;
    const uint64_t addr = reserve(bytes);
    int obj = add_dataset(path, nt, hdims, addr, bytes);

    std::vector<float> fbuf;
    const void *buf = data.data();
    if (tof && hila::myrank() == 0) {
        const A *p = reinterpret_cast<const A *>(data.data());
        fbuf.resize(nelem * ncomp);
        for (size_t i = 0; i < fbuf.size(); i++)
            fbuf[i] = (float)p[i];
        buf = fbuf.data();
    }
    write_rank0(addr, buf, bytes);

    return object_ref(this, obj);
}

} // namespace hila

#endif
