#ifndef HILA_FIELD_IO_H_
#define HILA_FIELD_IO_H_

//////////////////////////////////////////////////////////////////////
// This file collects Field<> I/O routines

#include "plumbing/field.h"
#include <limits>

#if defined(USE_FMT)
#define FMT_HEADER_ONLY
#include <fmt/format.h>
// #include <print>
#endif

namespace hila {
// Define this as template, in order to avoid code generation if the function is not needed -
// a bit crazy
template <typename IS, std::enable_if_t<std::is_same<IS, std::ifstream>::value, int> = 0>
bool open_input_file(const std::string &filename, IS &inputfile) {

    bool ok = true;
    if_rank0 () {
        inputfile.open(filename, std::ios::in | std::ios::binary);
        if (inputfile.fail()) {
            hila::out0 << "ERROR in opening file " << filename << '\n';
            ok = false;
        }
    }
    hila::broadcast(ok);
    if (!ok)
        hila::terminate(5);

    return ok;
}

template <typename OS, std::enable_if_t<std::is_same<OS, std::ofstream>::value, int> = 0>
bool open_output_file(const std::string &filename, OS &outputfile, bool binary = true,
                      int precision = 8) {
    bool ok = true;
    if_rank0 () {
        std::ios::openmode mode = std::ios::out | std::ios::trunc;
        if (binary)
            mode |= std::ios::binary;

        outputfile.open(filename, mode);
        outputfile.precision(precision);
        if (outputfile.fail()) {
            hila::out0 << "ERROR in opening file " << filename << '\n';
            ok = false;
        }
    }
    hila::broadcast(ok);
    if (!ok)
        hila::terminate(4);

    return ok;
}


template <typename S, std::enable_if_t<std::is_same<S, std::ofstream>::value ||
                                           std::is_same<S, std::ifstream>::value,
                                       int> = 0>
bool close_file(const std::string &filename, S &fstream) {
    bool error = false;
    if_rank0 () {
        fstream.close();
        if (fstream.fail()) {
            hila::out0 << "ERROR in reading/writing file " << filename << '\n';
            error = true;
        }
    }
    hila::broadcast(error);

    if (error)
        hila::terminate(3);

    return !error;
}

//////////////////////////////////////////////////////////////////////////////////
// Parallel (MPI-IO) binary I/O

#ifndef HILA_PARALLEL_IO
#define HILA_PARALLEL_IO 1
#endif

/**
 * @brief Select parallel binary Field I/O.
 * @details If true, binary Field I/O to named files (Field::write(filename),
 * Field::read(filename), write_fields(), read_fields(), GaugeField::write/read/config_write/
 * config_read) uses collective MPI-IO, each rank reading/writing its own part of the file.
 * If false, all data is routed through rank 0.  The file format is identical in both cases.
 * Default is set by HILA_PARALLEL_IO (default 1).  ASCII I/O and I/O to std::ofstream /
 * std::ifstream streams always go through rank 0.
 */
inline bool parallel_io = (HILA_PARALLEL_IO != 0);

namespace io {

/**
 * @brief Flags for Field output, combine with |, e.g. hila::io::hdf5 | hila::io::to_float
 */
enum output_flags : unsigned {
    ascii = 1,    ///< text output
    binary = 2,   ///< raw binary output
    hdf5 = 4,     ///< HDF5 file
    to_float = 8, ///< convert floating point numbers to single precision
};

/// @internal Collective error check for MPI-IO calls: if any rank failed, report and terminate
inline void check_mpi_io(int err, const char *what, const std::string &filename,
                         int exitcode = 1) {
    int code = (err == MPI_SUCCESS) ? 0 : err;
    int anycode = 0;
    MPI_Allreduce(&code, &anycode, 1, MPI_INT, MPI_MAX, lattice->mpi_comm_lat);
    if (anycode != 0) {
        char msg[MPI_MAX_ERROR_STRING];
        int len = 0;
        MPI_Error_string(anycode, msg, &len);
        hila::out0 << "MPI-IO ERROR " << what << " file " << filename << ": "
                   << std::string(msg, len) << '\n';
        hila::terminate(exitcode);
    }
}

/**
 * @internal Open file for collective MPI-IO.
 * @details For writing, an existing file is replaced: rank 0 deletes it and the file is
 * created anew.  On Lustre this is much faster than truncating a large file, and the new file
 * gets the striping of the directory (or of MPI-IO hints, e.g. MPICH_MPIIO_HINTS), which is
 * fixed at file creation.  If the file cannot be deleted it is truncated instead.
 * Note: if filename is a symbolic link, the link is replaced, not its target.
 */
inline MPI_File open_mpi_file(const std::string &filename, bool for_writing) {
    MPI_Comm comm = lattice->mpi_comm_lat;
    MPI_File fh = MPI_FILE_NULL;
    int amode = for_writing ? (MPI_MODE_WRONLY | MPI_MODE_CREATE) : MPI_MODE_RDONLY;

    bool truncate = false;
    if (for_writing) {
        // delete old file; fails harmlessly if it does not exist.  If it fails for some
        // other reason (e.g. no write permission to the directory), truncate after opening.
        if_rank0 () {
            truncate = (MPI_File_delete(filename.c_str(), MPI_INFO_NULL) != MPI_SUCCESS);
        }
        // this also guarantees that preceding rank-0-only file operations are done
        hila::broadcast(truncate);
    }

    int err = MPI_File_open(comm, filename.c_str(), amode, MPI_INFO_NULL, &fh);
    check_mpi_io(err, "opening", filename, for_writing ? 4 : 5);
    if (truncate) {
        err = MPI_File_set_size(fh, 0);
        check_mpi_io(err, "truncating", filename, 4);
    }
    return fh;
}

/// @internal Close MPI-IO file
inline void close_mpi_file(MPI_File &fh, const std::string &filename) {
    int err = MPI_File_close(&fh);
    check_mpi_io(err, "closing", filename, 3);
}

/// @internal MPI datatypes describing the rectangular block [start, start+subsize) of a
/// lattice of size `size`, element = elem_size bytes, x-direction running fastest.
struct block_filetype {
    MPI_Datatype etype, ftype;

    block_filetype(size_t elem_size, const CoordinateVector &size,
                   const CoordinateVector &subsize, const CoordinateVector &start) {
        int gs[NDIM], ss[NDIM], st[NDIM];
        foralldir (d) {
            gs[d] = size[d];
            ss[d] = subsize[d];
            st[d] = start[d];
        }
        MPI_Type_contiguous((int)elem_size, MPI_BYTE, &etype);
        MPI_Type_commit(&etype);
        MPI_Type_create_subarray(NDIM, gs, ss, st, MPI_ORDER_FORTRAN, etype, &ftype);
        MPI_Type_commit(&ftype);
    }

    ~block_filetype() {
        MPI_Type_free(&ftype);
        MPI_Type_free(&etype);
    }
};

/// @internal Reset the file view to the default (bytes from offset 0), so that explicit
/// offsets of MPI_File_write_at / read_at are byte offsets.  Collective.  All hila functions
/// which set a file view restore the default view afterwards.
inline void reset_file_view(MPI_File fh, const std::string &filename) {
    int err = MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
    check_mpi_io(err, "resetting view of", filename);
}

/// @internal Collective write of this rank's block [bmin, bmin+bsize) of a lattice-shaped
/// array of size gsize, starting at file offset `offset`.  buf holds the block in
/// lexicographic order (x fastest), elements of elem_size bytes.
/// If empty == true, this rank writes nothing (bsize, bmin and buf are ignored), but takes
/// part in the collective write.
inline void write_block_mpi(MPI_File fh, MPI_Offset offset, const void *buf, size_t elem_size,
                            const CoordinateVector &gsize, const CoordinateVector &bsize,
                            const CoordinateVector &bmin, const std::string &filename,
                            bool empty = false) {
    size_t nelem = 0;
    if (!empty) {
        nelem = 1;
        foralldir (d)
            nelem *= bsize[d];
    }
    assert(nelem <= (size_t)std::numeric_limits<int>::max());

    // MPI subarrays must be non-empty: use a 1-element block and write 0 elements
    CoordinateVector one, zero;
    one.fill(1);
    zero.fill(0);

    block_filetype bt(elem_size, gsize, empty ? one : bsize, empty ? zero : bmin);
    int err = MPI_File_set_view(fh, offset, bt.etype, bt.ftype, "native", MPI_INFO_NULL);
    check_mpi_io(err, "setting view of", filename);
    err = MPI_File_write_all(fh, buf, (int)nelem, bt.etype, MPI_STATUS_IGNORE);
    check_mpi_io(err, "writing", filename);
    reset_file_view(fh, filename);
}

} // namespace io

} // namespace hila

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename BC>
MPI_Offset Field<T, BC>::write_mpi(MPI_File fh, MPI_Offset offset,
                                   const std::string &filename) const {

    assert_all_ranks_present();
    assert(is_allocated() && "Field::write_mpi: Field not allocated");

    // copy_local_data() loops over the current lattice
    Lattice saved_lattice = lattice;
    lattice.switch_to(fs->mylattice);

    const auto &node = lattice->mynode;

    std::vector<T> buffer;
    copy_local_data(buffer);

    hila::io::write_block_mpi(fh, offset, buffer.data(), sizeof(T), lattice.size(), node.size,
                              node.min, filename);

    MPI_Offset end = offset + (MPI_Offset)lattice.volume() * (MPI_Offset)sizeof(T);

    lattice.switch_to(saved_lattice);
    return end;
}

template <typename T, typename BC>
MPI_Offset Field<T, BC>::read_mpi(MPI_File fh, MPI_Offset offset, const std::string &filename) {
    if (!this->is_allocated())
        this->allocate();
    return read_mpi(fh, offset, fs->mylattice.size(), filename);
}

template <typename T, typename BC>
MPI_Offset Field<T, BC>::read_mpi(MPI_File fh, MPI_Offset offset, const CoordinateVector &insize,
                                  const std::string &filename) {

    assert_all_ranks_present();
    if (!this->is_allocated())
        this->allocate();

    // set_local_data() loops over the current lattice
    Lattice saved_lattice = lattice;
    lattice.switch_to(fs->mylattice);

    const auto &node = lattice->mynode;
    const CoordinateVector lsize = lattice.size();

    // Region of the input lattice this rank needs.  If insize < lattice size, the input is
    // replicated: local coordinate c maps to c mod insize.  Read the bounding box of the
    // needed input coordinates; in directions where the node block wraps around insize
    // read the full extent.
    CoordinateVector rmin, rsize;
    size_t invol = 1, rvol = 1;
    bool direct = true;
    foralldir (d) {
        assert(insize[d] > 0 && insize[d] <= lsize[d]);
        invol *= insize[d];
        if (insize[d] == lsize[d]) {
            rmin[d] = node.min[d];
            rsize[d] = node.size[d];
        } else {
            direct = false;
            int a = node.min[d] % insize[d];
            if (a + node.size[d] <= insize[d]) {
                rmin[d] = a;
                rsize[d] = node.size[d];
            } else {
                rmin[d] = 0;
                rsize[d] = insize[d];
            }
        }
        rvol *= rsize[d];
    }
    assert(rvol <= (size_t)std::numeric_limits<int>::max());

    MPI_Offset end = offset + (MPI_Offset)invol * (MPI_Offset)sizeof(T);

    MPI_Offset fsize;
    int err = MPI_File_get_size(fh, &fsize);
    hila::io::check_mpi_io(err, "getting size of", filename);
    if (fsize < end) {
        hila::out0 << "ERROR: file " << filename << " too short: size " << fsize
                   << " bytes, reading Field needs " << end << '\n';
        hila::terminate(3);
    }

    std::vector<T> rbuf(rvol);
    hila::io::block_filetype bt(sizeof(T), insize, rsize, rmin);
    err = MPI_File_set_view(fh, offset, bt.etype, bt.ftype, "native", MPI_INFO_NULL);
    hila::io::check_mpi_io(err, "setting view of", filename);
    err = MPI_File_read_all(fh, rbuf.data(), (int)rvol, bt.etype, MPI_STATUS_IGNORE);
    hila::io::check_mpi_io(err, "reading", filename);
    hila::io::reset_file_view(fh, filename);

    if (direct) {
        set_local_data(rbuf);
    } else {
        // replicate: fill local block (lexicographic, x fastest) from rbuf
        std::vector<T> buffer(node.volume);
        CoordinateVector lc;
        lc.fill(0);
        for (size_t i = 0; i < node.volume; i++) {
            size_t ri = 0, mul = 1;
            foralldir (d) {
                ri += (size_t)((node.min[d] + lc[d]) % insize[d] - rmin[d]) * mul;
                mul *= rsize[d];
            }
            buffer[i] = rbuf[ri];
            foralldir (d) {
                if (++lc[d] < node.size[d])
                    break;
                lc[d] = 0;
            }
        }
        set_local_data(buffer);
    }

    lattice.switch_to(saved_lattice);
    return end;
}

//////////////////////////////////////////////////////////////////////////////////

#if defined(USE_FMT)

/// Write the field to a file stream
template <typename T, typename BC>
void Field<T, BC>::write(std::FILE *outputfile, int precision) const {

    constexpr size_t sites_per_write = WRITE_BUFFER_SIZE / sizeof(T);
    constexpr size_t write_size = sites_per_write * sizeof(T);

    assert_all_ranks_present();

    std::vector<CoordinateVector> coord_list(sites_per_write);
    T *buffer = nullptr;
    if_rank0 ()
        buffer = (T *)memalloc(write_size);

    auto mylat = fs->mylattice;
    CoordinateVector size = mylat.size();

    for (size_t i = 0; i < mylat.volume(); i += sites_per_write) {
        size_t sites = std::min<size_t>(sites_per_write, mylat.volume() - i);
        for (size_t j = 0; j < sites; j++)
            coord_list[j] = mylat->global_coordinates(i + j);

        if (sites < sites_per_write)
            coord_list.resize(sites);

        fs->gather_elements(buffer, coord_list);
        if_rank0 () {
            for (size_t j = 0; j < sites; j++) {
                // outputfile << buffer[j] << '\n';
                fmt::print(outputfile, "{}\n", buffer[j]);
            }
        }
    }

    if_rank0 ()
        std::free(buffer);
}

#endif

/// Write the field to a file stream
template <typename T, typename BC>
void Field<T, BC>::write(std::ofstream &outputfile, bool binary, int precision) const {
    constexpr size_t sites_per_write = WRITE_BUFFER_SIZE / sizeof(T);
    constexpr size_t write_size = sites_per_write * sizeof(T);

    assert_all_ranks_present();

    if (!binary)
        outputfile.precision(precision);

    std::vector<CoordinateVector> coord_list(sites_per_write);
    T *buffer = nullptr;
    if (hila::myrank() == 0)
        buffer = (T *)memalloc(write_size);

    auto mylat = fs->mylattice;
    CoordinateVector size = mylat.size();

    for (size_t i = 0; i < mylat.volume(); i += sites_per_write) {
        size_t sites = std::min<size_t>(sites_per_write, mylat.volume() - i);
        for (size_t j = 0; j < sites; j++)
            coord_list[j] = mylat->global_coordinates(i + j);

        if (sites < sites_per_write)
            coord_list.resize(sites);

        fs->gather_elements(buffer, coord_list);
        if_rank0 () {
            if (binary) {
                outputfile.write((char *)buffer, sites * sizeof(T));
            } else {
                for (size_t j = 0; j < sites; j++) {
                    outputfile << buffer[j] << '\n';
                }
            }
        }
    }

    if (hila::myrank() == 0)
        std::free(buffer);
}


/// Write the Field to a named file replacing the file
template <typename T, typename BC>
void Field<T, BC>::write(const std::string &filename, bool binary, int precision) const {

    if (binary && hila::parallel_io) {
        MPI_File fh = hila::io::open_mpi_file(filename, true);
        write_mpi(fh, 0, filename);
        hila::io::close_mpi_file(fh, filename);
        return;
    }

#if defined(USE_FMT)
    if (!binary) {
        // ASCII write through libfmt
        std::FILE *outputfile = nullptr;
        if_rank0 () {
            outputfile = std::fopen(filename.c_str(), "w");
        }
        bool ok = (outputfile != nullptr);
        hila::broadcast(ok);
        if (!ok) {
            hila::error("Could not open file " + filename);
        }

        write(outputfile, precision);
        if_rank0 () {
            fclose(outputfile);
        }
        return;
    }
#endif

    std::ofstream outputfile;
    hila::open_output_file(filename, outputfile, binary);
    write(outputfile, binary, precision);
    hila::close_file(filename, outputfile);
}

/// Write a list of fields into an output stream
template <typename T, typename BC>
static void write_fields(std::ofstream &outputfile, Field<T, BC> &last) {
    last.write(outputfile);
}

template <typename T, typename BC, typename... fieldtypes>
static void write_fields(std::ofstream &outputfile, Field<T, BC> &next, fieldtypes &...fields) {
    next.write(outputfile);
    write_fields(outputfile, fields...);
}

/// Write a list of fields to a file
template <typename... fieldtypes>
static void write_fields(const std::string &filename, fieldtypes &...fields) {
    if (hila::parallel_io) {
        MPI_File fh = hila::io::open_mpi_file(filename, true);
        MPI_Offset offset = 0;
        ((offset = fields.write_mpi(fh, offset, filename)), ...);
        hila::io::close_mpi_file(fh, filename);
    } else {
        std::ofstream outputfile;
        hila::open_output_file(filename, outputfile);
        write_fields(outputfile, fields...);
        hila::close_file(filename, outputfile);
    }
}

/////////////////////////////////////////////////////////////////////////////////

/// Read the Field from a stream
template <typename T, typename BC>
void Field<T, BC>::read(std::ifstream &inputfile) {
    constexpr size_t sites_per_read = WRITE_BUFFER_SIZE / sizeof(T);
    constexpr size_t read_size = sites_per_read * sizeof(T);

    if (!this->is_allocated())
        this->allocate();

    assert_all_ranks_present();
    will_change();
    mark_changed(ALL);

    std::vector<CoordinateVector> coord_list(sites_per_read);
    T *buffer = nullptr;
    if (hila::myrank() == 0)
        buffer = (T *)memalloc(read_size);

    auto mylat = fs->mylattice;
    CoordinateVector size = mylat.size();

    for (size_t i = 0; i < mylat.volume(); i += sites_per_read) {
        size_t sites = std::min<size_t>(sites_per_read, mylat.volume() - i);
        for (size_t j = 0; j < sites; j++)
            coord_list[j] = mylat->global_coordinates(i + j);

        if (sites < sites_per_read)
            coord_list.resize(sites);

        if_rank0 ()
            inputfile.read((char *)buffer, sites * sizeof(T));

        fs->scatter_elements(buffer, coord_list);
    }

    if (hila::myrank() == 0)
        std::free(buffer);
}

/// Read the Field from a stream
template <typename T, typename BC>
void Field<T, BC>::read(std::ifstream &inputfile, const CoordinateVector &insize) {
    constexpr size_t sites_per_read = WRITE_BUFFER_SIZE / sizeof(T);
    constexpr size_t read_size = sites_per_read * sizeof(T);

    if (!this->is_allocated())
        this->allocate();

    assert_all_ranks_present();
    will_change();
    mark_changed(ALL);

    std::vector<CoordinateVector> coord_list(sites_per_read);
    T *buffer = nullptr;
    if_rank0 ()
        buffer = (T *)memalloc(read_size);

    CoordinateVector lsize = lattice.size();
    int scalef[NDIM];
    size_t tvol = 1;
    size_t scvol = 1;
    foralldir (d) {
        scalef[d] = lsize[d] / insize[d];
        tvol *= insize[d];
        scvol *= scalef[d];
    }

    for (size_t i = 0; i < tvol; i += sites_per_read) {
        size_t sites = std::min(sites_per_read, tvol - i);
        if_rank0 () {
            inputfile.read((char *)buffer, sites * sizeof(T));
        }

        for (size_t j = 0; j < sites; j++) {
            size_t ind = i + j;
            foralldir (dir) {
                coord_list[j][dir] = ind % insize[dir];
                ind /= insize[dir];
            }
        }

        if (sites < sites_per_read) {
            coord_list.resize(sites);
        }

        fs->scatter_elements(buffer, coord_list);

        // if input lattice fits multiple times on the current lattice size,
        // replicate input lattice to fill current lattice:
        for (size_t sci = 1; sci < scvol; ++sci) {
            std::vector<CoordinateVector> tcoord_list(coord_list.size());
            size_t ind = sci;
            foralldir (dir) {
                int tsf = ind % scalef[dir];
                ind /= scalef[dir];
                for (size_t j = 0; j < sites; j++) {
                    tcoord_list[j][dir] = (coord_list[j][dir] + tsf * insize[dir]) % lsize[dir];
                }
            }
            fs->scatter_elements(buffer, tcoord_list);
        }
    }

    if_rank0 ()
        std::free(buffer);
}

// Read Field contents from the beginning of a file
template <typename T, typename BC>
void Field<T, BC>::read(const std::string &filename) {
    if (hila::parallel_io) {
        MPI_File fh = hila::io::open_mpi_file(filename, false);
        read_mpi(fh, 0, filename);
        hila::io::close_mpi_file(fh, filename);
    } else {
        std::ifstream inputfile;
        hila::open_input_file(filename, inputfile);
        read(inputfile);
        hila::close_file(filename, inputfile);
    }
}

// Read a list of fields from an input stream
template <typename T, typename BC>
static void read_fields(std::ifstream &inputfile, Field<T, BC> &last) {
    last.read(inputfile);
}

template <typename T, typename BC, typename... fieldtypes>
static void read_fields(std::ifstream &inputfile, Field<T, BC> &next, fieldtypes &...fields) {
    next.read(inputfile);
    read_fields(inputfile, fields...);
}

// Read a list of fields from a file
template <typename... fieldtypes>
static void read_fields(const std::string &filename, fieldtypes &...fields) {
    if (hila::parallel_io) {
        MPI_File fh = hila::io::open_mpi_file(filename, false);
        MPI_Offset offset = 0;
        ((offset = fields.read_mpi(fh, offset, filename)), ...);
        hila::io::close_mpi_file(fh, filename);
    } else {
        std::ifstream inputfile;
        hila::open_input_file(filename, inputfile);
        read_fields(inputfile, fields...);
        hila::close_file(filename, inputfile);
    }
}


///////////////////////////////////////////////////////////////////////////////

namespace hila {
/// @internal
/// Write n elements of type T from buf to an open output stream (called on rank 0 only).
/// ASCII: each element on its own line, numbers separated by spaces.
/// Binary: raw element data, optionally converted to float.
/// Used by Field::write_subvolume and Field::write_projection.
template <typename T>
void write_elements_to_stream(std::ofstream &outputfile, const T *buffer, size_t n, bool binary,
                              bool to_float) {

    constexpr int n_numbers = sizeof(T) / sizeof(hila::arithmetic_type<T>);

    if (!binary) {
        for (size_t k = 0; k < n; k++) {
            for (int l = 0; l < n_numbers; l++) {
                outputfile << hila::get_number_in_var(buffer[k], l) << ' ';
            }
            outputfile << '\n';
        }
    } else if (!to_float || std::is_same<hila::arithmetic_type<T>, float>::value) {
        outputfile.write((const char *)buffer, n * sizeof(T));
    } else {
        std::vector<float> floatbuf(n * n_numbers);
        float *fp = floatbuf.data();
        for (size_t k = 0; k < n; k++) {
            for (int l = 0; l < n_numbers; l++) {
                *(fp++) = hila::get_number_in_var(buffer[k], l);
            }
        }
        outputfile.write((const char *)floatbuf.data(), floatbuf.size() * sizeof(float));
    }
}
} // namespace hila

/// Write a "subspace" of the original lattice
/// Each element is written on a single line
/// TODO: more formatting?


template <typename T, typename BC>
void Field<T, BC>::write_subvolume_internal(std::ofstream &outputfile, const CoordinateVector &cmin,
                                            const CoordinateVector &cmax, int precision,
                                            bool binary, bool to_float) const {

    constexpr size_t sites_per_write = WRITE_BUFFER_SIZE / sizeof(T);

    assert_all_ranks_present();

    size_t sites = 1;
    int line_len = 1; // number of elements on 1st non-trivial dimension
    foralldir (d) {
        assert(cmin[d] >= 0 && cmax[d] >= cmin[d] && cmax[d] < fs->mylattice.size(d) &&
               "subvolume size mismatch");
        sites *= cmax[d] - cmin[d] + 1;

        if (line_len == 1)
            line_len = sites;
    }

    size_t n_write = std::min(sites_per_write, sites);

    std::vector<CoordinateVector> coord_list(n_write);
    T *buffer = nullptr;
    if_rank0 ()
        buffer = (T *)memalloc(n_write * sizeof(T));

    CoordinateVector c;

    size_t i = 0, j = 0;

    if (hila::myrank() == 0 && !binary) {
        outputfile.precision(precision);
    }

    forcoordinaterange(c, cmin, cmax) {
        coord_list[i] = c;
        i++;
        j++;
        if (i == n_write || j == sites) {
            if (i < n_write)
                coord_list.resize(i);

            fs->gather_elements(buffer, coord_list);

            if_rank0 () {
                hila::write_elements_to_stream(outputfile, buffer, i, binary, to_float);
            }
            i = 0;
        }
    }
    if_rank0 () {
        std::free(buffer);
    }
}

template <typename T, typename BC>
void Field<T, BC>::write_subvolume_internal(const std::string &filename,
                                            const CoordinateVector &cmin,
                                            const CoordinateVector &cmax, int precision,
                                            bool binary, bool to_float) const {

    std::ofstream out;
    hila::open_output_file(filename, out, false);
    write_subvolume_internal(out, cmin, cmax, precision, binary, to_float);
    hila::close_file(filename, out);
}

template <typename T, typename BC>
template <typename Outf>
void Field<T, BC>::write_subvolume(Outf &&output, const CoordinateVector &cmin,
                                   const CoordinateVector &cmax, int precision) const {
    write_subvolume_internal(output, cmin, cmax, precision, false, false);
}


template <typename T, typename BC>
template <typename Outf>
void Field<T, BC>::write_subvolume_binary(Outf &&output, const CoordinateVector &cmin,
                                          const CoordinateVector &cmax, bool to_float) const {
    write_subvolume_internal(output, cmin, cmax, 0, true, to_float);
}


template <typename T, typename BC>
template <typename Outf>
void Field<T, BC>::write_slice(Outf &&outf, const CoordinateVector &slice, int precision) const {

    auto [cmin, cmax] = get_range_from_slice(slice);
    write_subvolume(outf, cmin, cmax, precision);
}

template <typename T, typename BC>
template <typename outf_type>
void Field<T, BC>::write_slice_binary(outf_type &&outf, const CoordinateVector &slice,
                                      bool to_float) const {

    auto [cmin, cmax] = get_range_from_slice(slice);
    write_subvolume_binary(outf, cmin, cmax, to_float);
}


#endif
