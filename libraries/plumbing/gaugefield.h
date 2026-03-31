#ifndef HILA_GAUGEFIELD_H_
#define HILA_GAUGEFIELD_H_

/**
 * @file gaugefield.h
 * @brief Definition of Gauge Field
 * @details This file contains the definition for the GaugeField class
 */

#include "hila.h"

/**
 * @brief Gauge field class
 * @details Stores and defines links between Lattice Field elements. Number of links is
 * `lattice.size()*NDIM`, since for each point there is a link in all directions.
 *
 * @param fdir std::array<Field<T>,NDIM> type element which stores GaugeField links in back to back
 * direction wise ordering.
 *
 * @tparam T Group that GaugeField consists of
 */
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
        foralldir (d)
            fdir[d] = other[d];
    }

    // constructor with compatible scalar
    template <typename A, std::enable_if_t<hila::is_assignable<T &, A>::value, int> = 0>
    GaugeField(const A &val) {
        foralldir (d)
            fdir[d] = val;
    }

    // constructor from 0 - nullptr trick in use
    GaugeField(const std::nullptr_t z) {
        foralldir (d)
            fdir[d] = 0;
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
        foralldir (d)
            fdir[d] = val;
        return *this;
    }

    /// Separate 0 assignment
    GaugeField &operator=(std::nullptr_t np) {
        foralldir (d)
            fdir[d] = 0;
        return *this;
    }

    template <typename A>
    GaugeField &operator=(const GaugeField<A> &rhs) {
        foralldir (d)
            fdir[d] = rhs[d];
        return *this;
    }

    // This is now explicitly needed since it is implicitly deleted when move constructor is defined
    GaugeField &operator=(const GaugeField &rhs) {
        foralldir (d)
            fdir[d] = rhs[d];
        return *this;
    }

    GaugeField(GaugeField &&rhs) noexcept {
        foralldir (d)
            (*this)[d] = std::move(rhs[d]);
    }

    GaugeField &operator=(GaugeField &&rhs) noexcept {
        if (this != &rhs) {
            foralldir (d)
                (*this)[d] = std::move(rhs[d]);
        }
        return *this;
    }

    void clear() {
        foralldir (d)
            fdir[d].clear();
    }

    /**
     * @brief Reunitarize Gauge Field consisting of \f$ SU(N)\f$ matrices
     * @details Only defined for \f$ SU(N) \f$ matrices and Fields
     */
    void reunitarize_gauge() {
        foralldir (d) {
            onsites (ALL)
                (*this)[d][X].reunitarize();
        }
    }

    /**
     * @brief Computes Wilson action
     * @details \f{align}{ S &=  \beta\sum_{\textbf{dir}_1 < \textbf{dir}_2}\sum_{X} \frac{1}{N}
     * \Re\mathrm{Tr}\left[ 1- U_{\textbf{dir}_1 \textbf{dir}_2}(X) \right] \f} Where \f$\beta =
     * 2N/g^2\f$
     *
     * @return double
     */
    double measure_plaq() const {
        Reduction<double> plaq;
        plaq.allreduce(false);

        foralldir (dir1)
            foralldir (dir2)
                if (dir1 < dir2) {

                    onsites (ALL) {
                        plaq += 1.0 - real(trace((*this)[dir1][X] * (*this)[dir2][X + dir1] *
                                                 (*this)[dir1][X + dir2].dagger() *
                                                 (*this)[dir2][X].dagger())) /
                                          T::size();
                    }
                }

        return plaq.value();
    }
    ////////////////////////////////////////////////////////
    // I/O operations for gauge fields (here only binary)

    void write(std::ofstream &outputfile) const {
        foralldir (d) {
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
        foralldir (d) {
            fdir[d].read(inputfile);
        }
    }

    void read(std::ifstream &inputfile, CoordinateVector &insize) {
        foralldir (d) {
            fdir[d].read(inputfile, insize);
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

            foralldir (d) {
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
            foralldir (d) {
                inputfile.read(reinterpret_cast<char *>(&f), sizeof(int64_t));
                insize[d] = f;
                ok = ok && (lattice.size(d) % insize[d] == 0);
                if (!ok)
                    hila::out0 << conferr << "incorrect lattice dimension " << hila::prettyprint(d)
                               << " is " << f << " should be (or divide)" << lattice.size(d)
                               << '\n';
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

    /**
     * @brief Block the gauge field from parent gauge - form "long links" to connect blocked sites
     */

    void block_gauge(const GaugeField<T> &parent) {
        foralldir (d) {
            assert(parent[d].is_initialized(ALL));
            (*this)[d].check_alloc();
        }
        lattice_struct *blocklat = (*this)[e_x].fs->mylattice.ptr();
        lattice_struct *parentlat = parent[e_x].fs->mylattice.ptr();
        lattice_struct *currentlat = lattice.ptr();

        assert(blocklat->parent == parentlat && "blocking must happen from parent lattice Field");

        // alloc temp array, size of the blocked lattice
        size_t bufsize = blocklat->mynode.volume;
        T *buf = (T *)d_malloc(bufsize * sizeof(T));

        CoordinateVector blockfactor = parentlat->l_size.element_div(blocklat->l_size);
        CoordinateVector cvmin = blocklat->mynode.min;
        auto size_factor = blocklat->mynode.size_factor;

        foralldir (d) {
            assert(blockfactor[d] <= 2 &&
                   "block_gauge() can be used only with blocking factors 1 or 2");
        }

        foralldir (d) {
            // switch to parent
            lattice.switch_to(parentlat);

            // if no blocking to d, just copy the field
            if (blockfactor[d] == 1) {
                #pragma hila direct_access(buf)
                onsites (ALL) {
                    if (X.coordinates().is_divisible(blockfactor)) {
                        // get blocked coords logically on this
                        Vector<NDIM, unsigned> cv =
                            X.coordinates().element_div(blockfactor) - cvmin;
                        buf[cv.dot(size_factor)] = parent[d][X];
                    }
                }
            } else {
                // Now there is blocking by factor of 2, multiply long link
                #pragma hila direct_access(buf)
                onsites (ALL) {
                    if (X.coordinates().is_divisible(blockfactor)) {
                        // get blocked coords logically on this
                        Vector<NDIM, unsigned> cv =
                            X.coordinates().element_div(blockfactor) - cvmin;
                        buf[cv.dot(size_factor)] = parent[d][X] * parent[d][X + d];
                    }
                }
            }

            lattice.switch_to(blocklat);

            #pragma hila direct_access(buf)
            onsites (ALL) {
                // get blocked coords logically on this node
                Vector<NDIM, unsigned> cv = X.coordinates() - cvmin;
                (*this)[d][X] = buf[cv.dot(size_factor)];
            }
        } // directions

        lattice.switch_to(currentlat);

        d_free(buf);
    }

    /**
     * @brief Block the gauge field to the currently active lattice - form "long links" to connect
     * blocked sites
     * @details This Gaugefield has to be defined on the current active lattice's parent lattice and
     * the blocking factor has to be 1 or 2. Useful when blocking in a loop.
     */
    void block_gauge_to_current_lattice() {
        foralldir (d) {
            (*this)[d].check_alloc();
        }
        lattice_struct *thislat = (*this)[e_x].fs->mylattice.ptr();
        lattice_struct *currentlat = lattice.ptr();
        if (thislat == currentlat)
            return; // nothing to do

        // TODO: is this a necesary assumption?
        assert(currentlat->parent == thislat &&
               "blocking must happen to the next lattice blocked from thislat");

        // TODO: we can save this one alloc and one onsites(ALL) loop by constructing the
        // blocked field straight to the new field. This requires some hilapp magic to get the
        // indices right on each different platform. WIP.

        // alloc temp array, size of the blocked lattice
        size_t bufsize = currentlat->mynode.volume;
        T *buf = (T *)d_malloc(bufsize * sizeof(T));

        CoordinateVector blockfactor = thislat->l_size.element_div(currentlat->l_size);
        CoordinateVector cvmin = currentlat->mynode.min;
        auto size_factor = currentlat->mynode.size_factor;

        foralldir (d) {
            assert(blockfactor[d] <= 2 &&
                   "block_gauge() can be used only with blocking factors 1 or 2");
        }

        foralldir (d) {
            // switch to this fields lattice
            lattice.switch_to(thislat);

            // if no blocking to d, just copy the field
            if (blockfactor[d] == 1) {
                #pragma hila direct_access(buf)
                onsites (ALL) {
                    if (X.coordinates().is_divisible(blockfactor)) {
                        // get blocked coords logically on this
                        Vector<NDIM, unsigned> cv =
                            X.coordinates().element_div(blockfactor) - cvmin;
                        buf[cv.dot(size_factor)] = (*this)[d][X];
                    }
                }
            } else {
                // Now there is blocking by factor of 2, multiply long link
                #pragma hila direct_access(buf)
                onsites (ALL) {
                    if (X.coordinates().is_divisible(blockfactor)) {
                        // get blocked coords logically on this
                        Vector<NDIM, unsigned> cv =
                            X.coordinates().element_div(blockfactor) - cvmin;
                        buf[cv.dot(size_factor)] = (*this)[d][X] * (*this)[d][X + d];
                    }
                }
            }

            lattice.switch_to(currentlat);
            (*this)[d].clear();

            #pragma hila direct_access(buf)
            onsites (ALL) {
                // get blocked coords logically on this node
                Vector<NDIM, unsigned> cv = X.coordinates() - cvmin;
                (*this)[d][X] = buf[cv.dot(size_factor)];
            }
        } // directions

        d_free(buf);
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    void unblock_gauge(GaugeField<T> &target) const {
        foralldir (d) {
            (*this)[d].unblock_to(target[d]);
        }
    }
};


/// Implement hila::swap for gauge fields
namespace hila {
template <typename T>
void swap(GaugeField<T> &A, GaugeField<T> &B) {
    foralldir (d)
        hila::swap(A[d], B[d]);
}
} // namespace hila


///////////////////////////////////////////////////////
/// Alias VectorField to GaugeField
template <typename T>
using VectorField = GaugeField<T>;


#endif