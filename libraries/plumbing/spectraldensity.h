#ifndef SPECTRALDENSITY_H
#define SPECTRALDENSITY_H

#include "hila.h"
#include "fft.h"


// mod the k-space around origin, -L/2 -> L/2
inline Vector<NDIM, double> k_vector(const CoordinateVector &loc) {

    Vector<NDIM, double> k;
    foralldir (d) {
        int n = loc[d];
        if (n > lattice->size(d) / 2)
            n -= lattice->size(d);

        k[d] = n * 2.0 * M_PI / lattice->size(d);
    }
    return k;
}


/// Generic k-space field binner routine - returns field element type vector

template <typename T>
std::vector<T> bin_k_field(const Field<T> &f, int bins, double max_k = M_PI) {

    double mult = bins / max_k;

    VectorReduction<T> s(bins);
    s.allreduce(false);
    s = 0;

    onsites (ALL) {

        double kr = k_vector(X.coordinates()).norm();

        int b = kr * mult;
        if (b >= 0 && b < bins) {
            s[b] += f[X];
        }
    }

    std::vector<T> res(bins);

    for (int i = 0; i < bins; i++)
        res[i] = s[i];
    return res;
}

template <typename T, std::enable_if_t<hila::contains_complex<T>::value, int> = 0>
std::vector<double> bin_k_field_squarenorm(const Field<T> &f, int bins,
                                           double max_k = M_PI) {

    // This is for complex type
    using cmplx_t = Complex<hila::number_type<T>>;
    constexpr int n_cmplx = sizeof(T) / sizeof(cmplx_t);

    double mult = bins / max_k;

    VectorReduction<double> s(bins);
    s.allreduce(false);
    s = 0;

    onsites (ALL) {

        Vector<NDIM, double> k;
        k = k_vector(X.coordinates());
        double kr = k.norm();

        int b = kr * mult;
        if (b >= 0 && b < bins) {
            double ps = 0;
            for (int i = 0; i < n_cmplx; i++) {
                ps += hila::get_complex_element(f[X], i).squarenorm();
                // ps += c.squarenorm();
            }
            s[b] += ps;
        }
    }

    std::vector<double> res(bins);

    for (int i = 0; i < bins; i++)
        res[i] = s[i];
    return res;
}


//////////////////////////////////////////////////////////////////////////////////
/// Spectral density
/// This version takes in complex field

template <typename T, std::enable_if_t<hila::contains_complex<T>::value, int> = 0>
std::vector<double> spectraldensity(const Field<T> &f, int bins, double max_k = M_PI) {

    Field<T> ftrans;
    FFT_field(f, ftrans);

    return bin_k_field_squarenorm(ftrans, bins, max_k);
}

/// interface for real fields - an extra copy which could be avoided
template <typename T, std::enable_if_t<!hila::contains_complex<T>::value, int> = 0>
std::vector<double> spectraldensity(const Field<T> &f, int bins, double max_k = M_PI) {

    using cmplx_t = Complex<hila::number_type<T>>;

    if constexpr (sizeof(T) % sizeof(Complex<hila::number_type<T>>) == 0) {
        // real field, size is even -- cast the field to pseudo-complex
        // This works because layouts are compatible in all archs - if this changes then
        // need copy
        constexpr int nc = sizeof(T) / sizeof(cmplx_t);

        return spectraldensity(*reinterpret_cast<const Field<Vector<nc, cmplx_t>> *>(&f),
                             bins, max_k);
    } else {
        // now the size of input is not evenly divisible by sizeof complex.
        // new complex field
        constexpr int nc = sizeof(T) / sizeof(cmplx_t) + 1;

        Field<Vector<nc, cmplx_t>> cfield;
        onsites (ALL) {
            union {
                T tval;
                Vector<nc, cmplx_t> cvec;
            } u;
            u.tval = f[X];
            u.cvec.e(nc - 1).im = 0;
            cfield[X] = u.cvec;
        }

        return spectraldensity(cfield, bins, max_k);
    }
}


/// Data structure to hold the binning info - two vectors,
/// holding average k value in a bin and count of lattice points
struct binning_info {
    std::vector<double> k;
    std::vector<long> count;

    binning_info(int bins) : k(bins), count(bins) {}
};


/// Return the information about the binning  in binning_info 

inline binning_info bin_k_info(int bins, double max_k = M_PI) {

    double mult = bins / max_k;

    binning_info res(bins);

    VectorReduction<double> s(bins);
    VectorReduction<long> count(bins);
    s.allreduce(false);
    count.allreduce(false);
    s = 0;
    count = 0;

    onsites (ALL) {

        double kr = k_vector(X.coordinates()).norm();

        int b = kr * mult;
        if (b >= 0 && b < bins) {
            s[b] += kr;
            count[b] += 1;
        }
    }

    for (int i = 0; i < bins; i++) {
        res.count[i] = count[i];
        res.k[i] = s[i]/count[i];
    }

    return res;
}


#endif