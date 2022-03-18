#ifndef FIELD_IO_H_
#define FIELD_IO_H_

//////////////////////////////////////////////////////////////////////
// This file collects Field<> I/O routines

#include "plumbing/field.h"

/// Write the field to a file stream
template <typename T>
void Field<T>::write_to_stream(std::ofstream &outputfile) {
    constexpr size_t sites_per_write = WRITE_BUFFER_SIZE / sizeof(T);
    constexpr size_t write_size = sites_per_write * sizeof(T);

    std::vector<CoordinateVector> coord_list(sites_per_write);
    T *buffer = (T *)memalloc(write_size);
    CoordinateVector size = lattice->size();

    for (size_t i = 0; i < lattice->volume(); i += sites_per_write) {
        size_t sites = std::min(sites_per_write, lattice->volume() - i);
        for (size_t j = 0; j < sites; j++)
            coord_list[j] = lattice->global_coordinates(i + j);

        if (sites < sites_per_write)
            coord_list.resize(sites);

        fs->gather_elements(buffer, coord_list);
        if (hila::myrank() == 0)
            outputfile.write((char *)buffer, sites * sizeof(T));
    }

    std::free(buffer);
}

/// Write the Field to a named file replacing the file
template <typename T>
void Field<T>::write_to_file(const std::string &filename) {
    std::ofstream outputfile;
    outputfile.open(filename, std::ios::out | std::ios::trunc | std::ios::binary);
    write_to_stream(outputfile);
    outputfile.close();
}

/// Write a list of fields into an output stream
template <typename T>
static void write_fields(std::ofstream &outputfile, Field<T> &last) {
    last.write_to_stream(outputfile);
}

template <typename T, typename... fieldtypes>
static void write_fields(std::ofstream &outputfile, Field<T> &next,
                         fieldtypes &... fields) {
    next.write_to_stream(outputfile);
    write_fields(outputfile, fields...);
}

/// Write a list of fields to a file
template <typename... fieldtypes>
static void write_fields(const std::string &filename, fieldtypes &... fields) {
    std::ofstream outputfile;
    outputfile.open(filename, std::ios::out | std::ios::trunc | std::ios::binary);
    write_fields(outputfile, fields...);
    outputfile.close();
}

/////////////////////////////////////////////////////////////////////////////////

/// Read the Field from a stream
template <typename T>
void Field<T>::read_from_stream(std::ifstream &inputfile) {
    constexpr size_t sites_per_read = WRITE_BUFFER_SIZE / sizeof(T);
    constexpr size_t read_size = sites_per_read * sizeof(T);

    mark_changed(ALL);

    std::vector<CoordinateVector> coord_list(sites_per_read);
    T *buffer = (T *)memalloc(read_size);
    CoordinateVector size = lattice->size();

    for (size_t i = 0; i < lattice->volume(); i += sites_per_read) {
        size_t sites = std::min(sites_per_read, lattice->volume() - i);
        for (size_t j = 0; j < sites; j++)
            coord_list[j] = lattice->global_coordinates(i + j);

        if (sites < sites_per_read)
            coord_list.resize(sites);

        if (hila::myrank() == 0)
            inputfile.read((char *)buffer, sites * sizeof(T));

        fs->send_elements(buffer, coord_list);
    }

    std::free(buffer);
}

// Read Field contennts from the beginning of a file
template <typename T>
void Field<T>::read_from_file(const std::string &filename) {
    std::ifstream inputfile;
    inputfile.open(filename, std::ios::in | std::ios::binary);
    read_from_stream(inputfile);
    inputfile.close();
}

// Read a list of fields from an input stream
template <typename T>
static void read_fields(std::ifstream &inputfile, Field<T> &last) {
    last.read_from_stream(inputfile);
}

template <typename T, typename... fieldtypes>
static void read_fields(std::ifstream &inputfile, Field<T> &next,
                        fieldtypes &... fields) {
    next.read_from_stream(inputfile);
    read_fields(inputfile, fields...);
}

// Read a list of fields from a file
template <typename... fieldtypes>
static void read_fields(const std::string &filename, fieldtypes &... fields) {
    std::ifstream inputfile;
    inputfile.open(filename, std::ios::in | std::ios::binary);
    read_fields(inputfile, fields...);
    inputfile.close();
}


///////////////////////////////////////////////////////////////////////////////

/// Write a "subspace" of the original lattice
/// With separator defined to other string than '\n' the lowest (non-unity) dimension is
/// written on a single line (TODO: more formatting?)

template <typename T>
void Field<T>::write_subvolume(std::ofstream &outputfile, const CoordinateVector &cmin,
                               const CoordinateVector &cmax,
                               const std::string &separator) {

    constexpr size_t sites_per_write = WRITE_BUFFER_SIZE / sizeof(T);

    size_t sites = 1;
    int line_len = 1; // number of elements on 1st non-trivial dimension
    foralldir (d) {
        assert(cmin[d] >= 0 && cmax[d] >= cmin[d] && cmax[d] < lattice->size(d) &&
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
    int line_pos = 0;

    forcoordinaterange(c, cmin, cmax) {
        coord_list[i] = c;
        i++;
        j++;
        if (i == n_write || j == sites) {
            if (i < n_write)
                coord_list.resize(i);

            fs->gather_elements(buffer, coord_list);

            if (hila::myrank() == 0) {
                for (size_t k = 0; k < i; k++) {
                    outputfile << buffer[k];

                    // write separator if not at line end
                    line_pos++;
                    if (line_pos < line_len)
                        outputfile << separator;
                    else {
                        outputfile << '\n';
                        line_pos = 0;
                    }
                }
            }
            i = 0;
        }
    }
}


template <typename T>
void Field<T>::write_subvolume(const std::string &filename,
                               const CoordinateVector &cmin,
                               const CoordinateVector &cmax,
                               const std::string &separator) {

    std::ofstream out;
    if (hila::myrank() == 0) {
        out.open(filename, std::ios::out | std::ios::trunc | std::ios::binary);
    }

    write_subvolume(out, cmin, cmax, separator);
    int fail = 0;
    if (hila::myrank() == 0) {
        out.close();
        if (out.fail()) {
            hila::output << "Error in writing file " << filename << '\n';
            fail = 1;
        }
    }

    hila::broadcast(fail);
    if (fail)
        hila::terminate(1);
}

// Write a subspace (slice) of the original lattice.  Here
// slice[d] = fixed coordinate value, or < 0: run the dimension through.
// For example, slice = {5,2,-1,0}; prints the line along z-axis at (x,y,t) = (5,2,0)
// Outf is either filename or ofstream

template <typename T>
template <typename F>
void Field<T>::write_slice(F &outf, const CoordinateVector &slice,
                           const std::string &separator) {
    static_assert(std::is_same<F, std::ofstream>::value ||
                      std::is_same<F, std::string>::value,
                  "First argument must be std::ofstream or std::string");

    CoordinateVector cmin, cmax;
    foralldir (d) {
        if (slice[d] < 0) {
            cmin[d] = 0;
            cmax[d] = lattice->size(d) - 1;
        } else {
            cmin[d] = cmax[d] = slice[d];
        }
    }
    write_subvolume(outf, cmin, cmax, separator);
}


#endif
