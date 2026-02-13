/**
 * @file PlaquetteField.h
 * @brief Definition of Gauge Field
 * @details This file contains the definition for the PlaquetteField class
 */
#ifndef PlaquetteField_H_
#define PlaquetteField_H_

#include "hila.h"

/**
 * @brief Plaquette field class
 * @details Stores and defines plaquette Lattice Field elements. Number of plaquettes is
 * `lattice.size()*NDIM*NDIM` since at each point any pair of directions (dir1,dir2) can span a
 * plaquette. Degenerate plaquettes (dir1=dir2) are included for simplicity to be able to define
 * PlaquetteField in terms of VectorField.
 *
 * @param fdir std::array<VectorField<T>,NDIM> type element which stores PlaquetteField.
 *
 * @tparam T Group that PlaquetteField consists of
 */
template <typename T>
class PlaquetteField {
  private:
    std::array<VectorField<T>, NDIM> fdir;

    // somewhat arbitrary fingerprint flag for configuration files
    static constexpr int64_t config_flag = 394824242;

  public:
    // Default constructor
    PlaquetteField() = default;

    // Straightforward copy constructor seems to be necessary
    PlaquetteField(const PlaquetteField &other) = default;

    // copy constructor - from fields which can be assigned
    template <typename A, std::enable_if_t<std::is_convertible<A, T>::value, int> = 0>
    PlaquetteField(const PlaquetteField<A> &other) {
        foralldir(d1) fdir[d1] = other[d1];
    }

    // constructor with compatible Vector field
    template <typename A, std::enable_if_t<std::is_convertible<A, T>::value, int> = 0>
    PlaquetteField(const VectorField<A> &other) {
        foralldir(d1) fdir[d1] = other;
    }

    // constructor with compatible scalar
    template <typename A, std::enable_if_t<hila::is_assignable<T &, A>::value, int> = 0>
    PlaquetteField(const A &val) {
        foralldir(d1) fdir[d1] = val;
    }

    // constructor from 0 - nullptr trick in use
    PlaquetteField(const std::nullptr_t z) {
        foralldir(d1) fdir[d1] = 0;
    }


    /////////////////////////////////////////////////
    /// Destructor

    ~PlaquetteField() = default;

    /////////////////////////////////////////////////
    /// Access components with []

    inline VectorField<T> &operator[](Direction d1) {
        return fdir[d1];
    }

    inline const VectorField<T> &operator[](Direction d1) const {
        return fdir[d1];
    }

    //////////////////////////////////////////////////
    /// Assign from anything the field allows
    template <typename A>
    PlaquetteField &operator=(const A &val) {
        foralldir(d1) fdir[d1] = val;
        return *this;
    }

    /// Separate 0 assignment
    PlaquetteField &operator=(std::nullptr_t np) {
        foralldir(d1) fdir[d1] = 0;
        return *this;
    }

    template <typename A>
    PlaquetteField &operator=(const PlaquetteField<A> &rhs) {
        foralldir(d1) fdir[d1] = rhs[d1];
        return *this;
    }

    ////////////////////////////////////////////////////////
    // I/O operations for plaquette fields (here only binary)

    void write(std::ofstream &outputfile) const {
        foralldir(d1) {
            fdir[d1].write(outputfile);
        }
    }

    void write(const std::string &filename) const {
        std::ofstream outputfile;
        hila::open_output_file(filename, outputfile);
        write(outputfile);
        hila::close_file(filename, outputfile);
    }

    void read(std::ifstream &inputfile) {
        foralldir(d1) {
            fdir[d1].read(inputfile);
        }
    }

    void read(std::ifstream &inputfile, CoordinateVector &insize) {
        foralldir(d1) {
            fdir[d1].read(inputfile, insize);
        }
    }

    void read(const std::string &filename) {
        std::ifstream inputfile;
        hila::open_input_file(filename, inputfile);
        read(inputfile);
        hila::close_file(filename, inputfile);
    }

    /// config_write writes the gauge field to file, with additional "verifying" header

    void config_write(const std::string &filename) const {
        std::ofstream outputfile;
        hila::open_output_file(filename, outputfile);

        // write header
        if_rank0() {
            int64_t f = config_flag;
            outputfile.write(reinterpret_cast<char *>(&f), sizeof(int64_t));
            f = NDIM;
            outputfile.write(reinterpret_cast<char *>(&f), sizeof(int64_t));
            f = sizeof(T);
            outputfile.write(reinterpret_cast<char *>(&f), sizeof(int64_t));

            foralldir(d) {
                f = lattice.size(d);
                outputfile.write(reinterpret_cast<char *>(&f), sizeof(int64_t));
            }
        }

        write(outputfile);
        hila::close_file(filename, outputfile);
    }

    void config_read(const std::string &filename) {
        std::ifstream inputfile;
        hila::open_input_file(filename, inputfile);
        std::string conferr("CONFIG ERROR in file " + filename + ": ");

        // read header
        bool ok = true;
        int64_t f;
        if_rank0() {
            inputfile.read(reinterpret_cast<char *>(&f), sizeof(int64_t));
            ok = (f == config_flag);
            if (!ok)
                hila::out0 << conferr << "wrong id, should be " << config_flag << " is " << f
                           << '\n';
        }

        if (ok && hila::myrank() == 0) {
            inputfile.read(reinterpret_cast<char *>(&f), sizeof(int64_t));
            ok = (f == NDIM);
            if (!ok)
                hila::out0 << conferr << "wrong dimensionality, should be " << NDIM << " is " << f
                           << '\n';
        }

        if (ok && hila::myrank() == 0) {
            inputfile.read(reinterpret_cast<char *>(&f), sizeof(int64_t));
            ok = (f == sizeof(T));
            if (!ok)
                hila::out0 << conferr << "wrong size of field element, should be " << sizeof(T)
                           << " is " << f << '\n';
        }

        CoordinateVector insize = lattice.size();
        if (ok && hila::myrank() == 0) {

            foralldir(d) {
                inputfile.read(reinterpret_cast<char *>(&f), sizeof(int64_t));
                insize[d] = f;
                ok = ok && (lattice.size(d) % insize[d] == 0);
                if (!ok)
                    hila::out0 << conferr << "incorrect lattice dimension " << hila::prettyprint(d)
                               << " is " << f << " should be (or divide) " << lattice.size(d) << '\n';
            }
        }

        if (!hila::broadcast(ok)) {
            hila::terminate(1);
        } else {
            hila::broadcast_array(insize.c, NDIM);
        }

        read(inputfile, insize);
        hila::close_file(filename, inputfile);
    }
};


/// Implement hila::swap for gauge fields
namespace hila {
template <typename T>
void swap(PlaquetteField<T> &A, PlaquetteField<T> &B) {
    foralldir(d1) hila::swap(A[d1], B[d1]);
}
} // namespace std



#endif