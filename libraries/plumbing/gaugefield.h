#ifndef GAUGEFIELD_H_
#define GAUGEFIELD_H_

#include "hila.h"

template <typename T>
class GaugeField {
  private:
    std::array<Field<T>, NDIM> fdir;

    // somewhat arbitrary fingerprint flag for configuration files
    static constexpr int64_t config_flag = 394824242;

  public:
    // Default constructor
    GaugeField() = default;

    // Straightforward copy constructor seems to be necessary
    GaugeField(const GaugeField &other) = default;

    // copy constructor - from fields which can be assigned
    template <typename A, std::enable_if_t<std::is_convertible<A, T>::value, int> = 0>
    GaugeField(const GaugeField<A> &other) {
        foralldir(d) fdir[d] = other[d];
    }

    // constructor with compatible scalar
    template <typename A, std::enable_if_t<hila::is_assignable<T &, A>::value, int> = 0>
    GaugeField(const A &val) {
        foralldir(d) fdir[d] = val;
    }

    // constructor from 0 - nullptr trick in use
    GaugeField(const std::nullptr_t z) {
        foralldir(d) fdir[d] = 0;
    }


    /////////////////////////////////////////////////
    /// Destructor

    ~GaugeField() = default;

    /////////////////////////////////////////////////
    /// Access components with []

    inline Field<T> &operator[](Direction d) {
        return fdir[d];
    }

    inline const Field<T> &operator[](Direction d) const {
        return fdir[d];
    }

    //////////////////////////////////////////////////
    /// Assign from anything the field allows
    template <typename A>
    GaugeField &operator=(const A &val) {
        foralldir(d) fdir[d] = val;
        return *this;
    }

    /// Separate 0 assignment
    GaugeField &operator=(std::nullptr_t np) {
        foralldir(d) fdir[d] = 0;
        return *this;
    }

    template <typename A>
    GaugeField &operator=(const GaugeField<A> &rhs) {
        foralldir(d) fdir[d] = rhs[d];
        return *this;
    }

    ////////////////////////////////////////////////////////
    /// I/O operations for gauge fields (here only binary)

    void write(std::ofstream &outputfile) const {
        foralldir(d) {
            fdir[d].write(outputfile);
        }
    }

    void write(const std::string &filename) const {
        std::ofstream outputfile;
        hila::open_output_file(filename, outputfile);
        write(outputfile);
        hila::close_file(filename, outputfile);
    }

    void read(std::ifstream &inputfile) {
        foralldir(d) {
            fdir[d].read(inputfile);
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
        if (hila::myrank() == 0) {
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
        if (hila::myrank() == 0) {
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

        if (ok && hila::myrank() == 0) {

            foralldir(d) {
                inputfile.read(reinterpret_cast<char *>(&f), sizeof(int64_t));
                ok = ok && (f == lattice.size(d));
                if (!ok)
                    hila::out0 << conferr << "incorrect lattice dimension " << hila::prettyprint(d)
                               << " is " << f << " should be " << lattice.size(d) << '\n';
            }
        }

        if (!hila::broadcast(ok)) {
            hila::terminate(1);
        }

        read(inputfile);
        hila::close_file(filename, inputfile);
    }
};


/// Implement std::swap for gauge fields
namespace std {
template <typename T>
void swap(GaugeField<T> &A, GaugeField<T> &B) {
    foralldir(d) std::swap(A[d], B[d]);
}
} // namespace std


///////////////////////////////////////////////////////
/// Alias VectorField to GaugeField
template <typename T>
using VectorField = GaugeField<T>;


#endif