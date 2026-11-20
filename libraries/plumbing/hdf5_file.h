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
 *   h5.write_subvolume("phi_box", phi, cmin, cmax);       // attributes cmin, cmax
 *   h5.write_slice("phi_xy", phi, {-1, -1, 5});           // 2-dim dataset [y, x]
 *   h5.write_projection("phi_z", phi, {e_x, e_y});        // sum over x and y -> [z]
 *   h5.close();                                           // or at destruction
 * @endcode
 * All ranks must call all methods.  Attribute values and std::vector data are taken from rank 0.
 */

#include "plumbing/field_io.h"
#include <cstdint>
#include <cstring>
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

    /**
     * @brief Write the subvolume [cmin, cmax] (corners included) of a Field.  Collective,
     * each rank writes its own part with MPI-IO.
     * @details Dataset shape [cmax-cmin+1 for directions NDIM-1, ..., 0, ncomp].  Attributes
     * "cmin" and "cmax" (x first).
     */
    template <typename T, typename BC>
    object_ref write_subvolume(const std::string &path, const Field<T, BC> &f,
                               const CoordinateVector &cmin, const CoordinateVector &cmax,
                               unsigned flags = 0);

    /**
     * @brief Write a slice of a Field: slice[d] >= 0 fixes the coordinate, < 0 runs through
     * it (as in Field::write_slice()).  Collective.
     * @details The dataset has only the running directions, e.g. slice {-1, 3, -1} gives a
     * 2-dimensional dataset [L_z, L_x, ncomp].  Attribute "slice".
     */
    template <typename T, typename BC>
    object_ref write_slice(const std::string &path, const Field<T, BC> &f,
                           const CoordinateVector &slice, unsigned flags = 0);

    /**
     * @brief Write the projection (sum over directions dirs) of a Field, see
     * Field::get_projection().  Collective; rank 0 writes the result slab by slab.
     * @details The dataset has the kept directions (highest first, as for Fields), or one
     * element if all directions are projected.  Attributes "projected directions" (direction
     * indices, x = 0), "cmin" and "cmax".  The Field must belong to the current lattice.
     */
    template <typename T, typename BC>
    object_ref write_projection(const std::string &path, const Field<T, BC> &f,
                                const std::vector<Direction> &dirs, unsigned flags = 0);
    template <typename T, typename BC>
    object_ref write_projection(const std::string &path, const Field<T, BC> &f,
                                const std::vector<Direction> &dirs, const CoordinateVector &cmin,
                                const CoordinateVector &cmax, unsigned flags = 0);

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
    int write_at(uint64_t addr, const void *buf, uint64_t bytes);

    /// number type of T in the file, with to_float flag
    template <typename T>
    static number_type file_number_type(unsigned flags) {
        using A = std::conditional_t<std::is_arithmetic<T>::value, T, hila::arithmetic_type<T>>;
        static_assert(std::is_arithmetic<A>::value,
                      "HDF5 output: element must consist of arithmetic numbers");
        static_assert(sizeof(T) % sizeof(A) == 0, "HDF5 output: element has padding");
        if ((flags & hila::io::to_float) && std::is_floating_point<A>::value && sizeof(A) > 4)
            return number_type_of<float>();
        return number_type_of<A>();
    }

    /// number of arithmetic numbers in T
    template <typename T>
    static constexpr size_t n_components() {
        using A = std::conditional_t<std::is_arithmetic<T>::value, T, hila::arithmetic_type<T>>;
        return sizeof(T) / sizeof(A);
    }

    /// Append n elements of T as numbers of type nt (converted to float if needed) to out
    template <typename T>
    static void append_numbers(std::vector<uint8_t> &out, const T *src, size_t n,
                               number_type nt) {
        using A = std::conditional_t<std::is_arithmetic<T>::value, T, hila::arithmetic_type<T>>;
        const size_t m = n * n_components<T>();
        const size_t start = out.size();
        if (nt.size == sizeof(A)) {
            const uint8_t *p = reinterpret_cast<const uint8_t *>(src);
            out.insert(out.end(), p, p + m * sizeof(A));
        } else {
            const A *p = reinterpret_cast<const A *>(src);
            out.resize(start + m * sizeof(float));
            for (size_t i = 0; i < m; i++) {
                float v = (float)p[i];
                std::memcpy(out.data() + start + i * sizeof(float), &v, sizeof(float));
            }
        }
    }

    /// Dataset dimensions in file order from sizes x-first, adding ncomp
    static std::vector<uint64_t> file_dims(const std::vector<uint64_t> &sizes_x_first,
                                           size_t ncomp) {
        std::vector<uint64_t> dims(sizes_x_first.rbegin(), sizes_x_first.rend());
        if (ncomp > 1 || dims.empty())
            dims.push_back(ncomp);
        return dims;
    }

    template <typename T, typename BC>
    int write_box(const std::string &path, const Field<T, BC> &f, const CoordinateVector &cmin,
                  const CoordinateVector &cmax, const CoordinateVector &keep, unsigned flags);

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

/// Subvolume [cmin, cmax] of f; dataset has the directions with keep[d] != 0
template <typename T, typename BC>
int hdf5_file::write_box(const std::string &path, const Field<T, BC> &f,
                         const CoordinateVector &cmin, const CoordinateVector &cmax,
                         const CoordinateVector &keep, unsigned flags) {

    constexpr size_t ncomp = n_components<T>();
    const number_type nt = file_number_type<T>(flags);

    if (!f.is_allocated())
        hila::error("hdf5_file::write: Field for dataset " + path + " is not allocated");

    // copy_local_data() loops over the current lattice
    Lattice saved_lattice = lattice;
    Lattice flat = f.mylattice();
    lattice.switch_to(flat);
    const auto &node = lattice->mynode;

    // box size, and the part of it on this rank: [lo, hi]
    CoordinateVector bsize, lo, hi, isize, rel;
    std::vector<uint64_t> sizes;
    uint64_t nelem = 1;
    bool empty = false;
    foralldir (d) {
        if (cmin[d] < 0 || cmax[d] < cmin[d] || cmax[d] >= lattice.size(d))
            hila::error("hdf5_file " + filename + ": range error for dataset " + path);
        bsize[d] = cmax[d] - cmin[d] + 1;
        nelem *= bsize[d];
        if (keep[d])
            sizes.push_back(bsize[d]);
        lo[d] = std::max(cmin[d], node.min[d]);
        hi[d] = std::min(cmax[d], node.min[d] + node.size[d] - 1);
        if (hi[d] < lo[d])
            empty = true;
        isize[d] = hi[d] - lo[d] + 1;
        rel[d] = lo[d] - cmin[d];
    }

    const uint64_t bytes = nelem * ncomp * nt.size;
    const uint64_t addr = reserve(bytes);
    int obj = add_dataset(path, nt, file_dims(sizes, ncomp), addr, bytes);

    std::vector<T> buf;
    f.copy_local_data(buf);

    // collect this rank's part in lexicographic order, x-rows are contiguous in buf
    std::vector<uint8_t> out;
    if (!empty) {
        size_t rows = 1;
        foralldir (d)
            if (d > e_x)
                rows *= isize[d];
        CoordinateVector c = lo;
        for (size_t r = 0; r < rows; r++) {
            size_t idx = 0;
            foralldir (d)
                idx += (size_t)(c[d] - node.min[d]) * node.size_factor[d];
            append_numbers(out, buf.data() + idx, isize[e_x], nt);
            foralldir (d) {
                if (d == e_x)
                    continue;
                if (++c[d] <= hi[d])
                    break;
                c[d] = lo[d];
            }
        }
    }

    hila::io::write_block_mpi(fh, addr, out.data(), ncomp * nt.size, bsize, isize, rel, filename,
                              empty);

    lattice.switch_to(saved_lattice);
    return obj;
}

template <typename T, typename BC>
hdf5_file::object_ref hdf5_file::write_subvolume(const std::string &path, const Field<T, BC> &f,
                                                 const CoordinateVector &cmin,
                                                 const CoordinateVector &cmax, unsigned flags) {
    CoordinateVector keep;
    keep.fill(1);
    object_ref r(this, write_box(path, f, cmin, cmax, keep, flags));
    r.attribute("cmin", cmin).attribute("cmax", cmax);
    return r;
}

template <typename T, typename BC>
hdf5_file::object_ref hdf5_file::write_slice(const std::string &path, const Field<T, BC> &f,
                                             const CoordinateVector &slice, unsigned flags) {
    if (!f.is_allocated())
        hila::error("hdf5_file::write_slice: Field for dataset " + path + " is not allocated");
    CoordinateVector cmin, cmax, keep;
    foralldir (d) {
        if (slice[d] >= 0) {
            cmin[d] = cmax[d] = slice[d];
            keep[d] = 0;
        } else {
            cmin[d] = 0;
            cmax[d] = f.mylattice().size(d) - 1;
            keep[d] = 1;
        }
    }
    object_ref r(this, write_box(path, f, cmin, cmax, keep, flags));
    r.attribute("slice", slice);
    return r;
}

template <typename T, typename BC>
hdf5_file::object_ref hdf5_file::write_projection(const std::string &path, const Field<T, BC> &f,
                                                  const std::vector<Direction> &dirs,
                                                  const CoordinateVector &cmin,
                                                  const CoordinateVector &cmax, unsigned flags) {

    constexpr size_t ncomp = n_components<T>();
    const number_type nt = file_number_type<T>(flags);

    CoordinateVector is_proj;
    is_proj.fill(0);
    std::vector<int> pdirs;
    for (auto d : dirs) {
        if (!is_up_dir(d) || is_proj[d])
            hila::error("hdf5_file::write_projection: invalid directions for dataset " + path);
        is_proj[d] = 1;
        pdirs.push_back((int)d);
    }

    std::vector<uint64_t> sizes;
    uint64_t nelem = 1;
    foralldir (d) {
        if (cmin[d] < 0 || cmax[d] < cmin[d] || cmax[d] >= lattice.size(d))
            hila::error("hdf5_file " + filename + ": range error for dataset " + path);
        if (!is_proj[d]) {
            sizes.push_back(cmax[d] - cmin[d] + 1);
            nelem *= cmax[d] - cmin[d] + 1;
        }
    }

    const uint64_t bytes = nelem * ncomp * nt.size;
    const uint64_t addr = reserve(bytes);
    int obj = add_dataset(path, nt, file_dims(sizes, ncomp), addr, bytes);

    // projection_internal calls the output function on rank 0 only, slab by slab in order
    uint64_t pos = addr;
    int err = MPI_SUCCESS;
    std::vector<uint8_t> out;
    f.projection_internal(dirs, cmin, cmax, [&](const T *data, size_t n) {
        out.clear();
        append_numbers(out, data, n, nt);
        if (err == MPI_SUCCESS)
            err = write_at(pos, out.data(), out.size());
        pos += out.size();
    });
    hila::io::check_mpi_io(err, "writing", filename);

    object_ref r(this, obj);
    r.attribute("projected directions", pdirs).attribute("cmin", cmin).attribute("cmax", cmax);
    return r;
}

template <typename T, typename BC>
hdf5_file::object_ref hdf5_file::write_projection(const std::string &path, const Field<T, BC> &f,
                                                  const std::vector<Direction> &dirs,
                                                  unsigned flags) {
    CoordinateVector cmin, cmax;
    foralldir (d) {
        cmin[d] = 0;
        cmax[d] = lattice.size(d) - 1;
    }
    return write_projection(path, f, dirs, cmin, cmax, flags);
}

} // namespace hila

#endif
