#ifndef HILA_FIELD_IO_H_
#define HILA_FIELD_IO_H_

//////////////////////////////////////////////////////////////////////
// This file collects Field<> I/O routines

#include "plumbing/field.h"

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

} // namespace hila

//////////////////////////////////////////////////////////////////////////////////

/// Write the field to a file stream
template <typename T>
void Field<T>::write(std::ofstream &outputfile, bool binary, int precision) const {
    constexpr size_t sites_per_write = WRITE_BUFFER_SIZE / sizeof(T);
    constexpr size_t write_size = sites_per_write * sizeof(T);

    assert_all_ranks_present();

    if (!binary)
        outputfile.precision(precision);

    std::vector<CoordinateVector> coord_list(sites_per_write);
    T *buffer = (T *)memalloc(write_size);
    auto mylat = fs->mylattice;
    CoordinateVector size = mylat.size();

    for (size_t i = 0; i < mylat.volume(); i += sites_per_write) {
        size_t sites = std::min(sites_per_write, mylat.volume() - i);
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

    std::free(buffer);
}

/// Write the Field to a named file replacing the file
template <typename T>
void Field<T>::write(const std::string &filename, bool binary, int precision) const {
    std::ofstream outputfile;
    hila::open_output_file(filename, outputfile, binary);
    write(outputfile, binary, precision);
    hila::close_file(filename, outputfile);
}

/// Write a list of fields into an output stream
template <typename T>
static void write_fields(std::ofstream &outputfile, Field<T> &last) {
    last.write(outputfile);
}

template <typename T, typename... fieldtypes>
static void write_fields(std::ofstream &outputfile, Field<T> &next, fieldtypes &...fields) {
    next.write(outputfile);
    write_fields(outputfile, fields...);
}

/// Write a list of fields to a file
template <typename... fieldtypes>
static void write_fields(const std::string &filename, fieldtypes &...fields) {
    std::ofstream outputfile;
    hila::open_output_file(filename, outputfile);
    write_fields(outputfile, fields...);
    hila::close_file(filename, outputfile);
}

/////////////////////////////////////////////////////////////////////////////////

/// Read the Field from a stream
template <typename T>
void Field<T>::read(std::ifstream &inputfile) {
    constexpr size_t sites_per_read = WRITE_BUFFER_SIZE / sizeof(T);
    constexpr size_t read_size = sites_per_read * sizeof(T);

    if (!this->is_allocated())
        this->allocate();

    assert_all_ranks_present();
    will_change();
    mark_changed(ALL);

    std::vector<CoordinateVector> coord_list(sites_per_read);
    T *buffer = (T *)memalloc(read_size);
    auto mylat = fs->mylattice;
    CoordinateVector size = mylat.size();

    for (size_t i = 0; i < mylat.volume(); i += sites_per_read) {
        size_t sites = std::min(sites_per_read, mylat.volume() - i);
        for (size_t j = 0; j < sites; j++)
            coord_list[j] = mylat->global_coordinates(i + j);

        if (sites < sites_per_read)
            coord_list.resize(sites);

        if_rank0 ()
            inputfile.read((char *)buffer, sites * sizeof(T));

        fs->scatter_elements(buffer, coord_list);
    }

    std::free(buffer);
}

/// Read the Field from a stream
template <typename T>
void Field<T>::read(std::ifstream &inputfile, const CoordinateVector &insize) {
    constexpr size_t sites_per_read = WRITE_BUFFER_SIZE / sizeof(T);
    constexpr size_t read_size = sites_per_read * sizeof(T);

    if (!this->is_allocated())
        this->allocate();

    assert_all_ranks_present();
    will_change();
    mark_changed(ALL);

    std::vector<CoordinateVector> coord_list(sites_per_read);
    T *buffer = (T *)memalloc(read_size);

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

    std::free(buffer);
}

// Read Field contents from the beginning of a file
template <typename T>
void Field<T>::read(const std::string &filename) {
    std::ifstream inputfile;
    hila::open_input_file(filename, inputfile);
    read(inputfile);
    hila::close_file(filename, inputfile);
}

// Read a list of fields from an input stream
template <typename T>
static void read_fields(std::ifstream &inputfile, Field<T> &last) {
    last.read_from_stream(inputfile);
}

template <typename T, typename... fieldtypes>
static void read_fields(std::ifstream &inputfile, Field<T> &next, fieldtypes &...fields) {
    next.read_from_stream(inputfile);
    read_fields(inputfile, fields...);
}

// Read a list of fields from a file
template <typename... fieldtypes>
static void read_fields(const std::string &filename, fieldtypes &...fields) {
    std::ifstream inputfile;
    hila::open_input_file(filename, inputfile);
    read_fields(inputfile, fields...);
    hila::close_file(filename, inputfile);
}


///////////////////////////////////////////////////////////////////////////////

/// Write a "subspace" of the original lattice
/// Each element is written on a single line
/// TODO: more formatting?

template <typename T>
void Field<T>::write_subvolume(std::ofstream &outputfile, const CoordinateVector &cmin,
                               const CoordinateVector &cmax, int precision) const {

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
    T *buffer = (T *)memalloc(n_write * sizeof(T));

    CoordinateVector c;

    size_t i = 0, j = 0;

    if_rank0 () {
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
                for (size_t k = 0; k < i; k++) {
                    for (int l = 0; l < sizeof(T) / sizeof(hila::arithmetic_type<T>); l++) {
                        outputfile << hila::get_number_in_var(buffer[k], l) << ' ';
                    }
                    outputfile << '\n';
                }
            }
            i = 0;
        }
    }
}


template <typename T>
void Field<T>::write_subvolume(const std::string &filename, const CoordinateVector &cmin,
                               const CoordinateVector &cmax, int precision) const {

    std::ofstream out;
    hila::open_output_file(filename, out, false);
    write_subvolume(out, cmin, cmax, precision);
    hila::close_file(filename, out);
}

// Write a subspace (slice) of the original lattice.  Here
// slice[d] = fixed coordinate value, or < 0: run the dimension through.
// For example, slice = {5,2,-1,0}; prints the line along z-axis at (x,y,t) = (5,2,0)
// Outf is either filename or ofstream

template <typename T>
template <typename outf_type>
void Field<T>::write_slice(outf_type &outf, const CoordinateVector &slice, int precision) const {

    static_assert(std::is_same<outf_type, std::string>::value ||
                      std::is_same<outf_type, std::ofstream>::value,
                  "file name / output stream argument in write_slice()?");

    CoordinateVector cmin, cmax;
    foralldir (d) {
        if (slice[d] < 0) {
            cmin[d] = 0;
            cmax[d] = fs->mylattice.size(d) - 1;
        } else {
            cmin[d] = cmax[d] = slice[d];
        }
    }
    write_subvolume(outf, cmin, cmax, precision);
}


#endif
