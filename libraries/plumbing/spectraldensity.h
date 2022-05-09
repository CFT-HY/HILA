#ifndef POWERSPECTRUM_H
#define POWERSPECTRUM_H

#include "hila.h"
#include "fft.h"

//////////////////////////////////////////////////////////////////////////////////
/// Spectral density
/// This version takes in complex field

template <typename T, std::enable_if_t<hila::contains_complex<T>::value, int> = 0>
std::vector<double> powerspectrum(const Field<T> &f, int bins, double max_k = M_PI) {

    // This is for complex type
    using cmplx_t = Complex<hila::number_type<T>>;
    constexpr int n_cmplx = sizeof(T) / sizeof(cmplx_t);

    Field<T> ftrans;
    FFT_field(f, ftrans);

    double mult = bins / max_k;

    VectorReduction<double> s(bins);
    s.allreduce(false);
    s = 0;

    onsites (ALL) {

        // mod the k-space around origin, -L/2 -> L/2
        CoordinateVector n;
        Vector<NDIM, double> k;

        foralldir (d) {
            int n = X.coordinate(d);
            if (n > lattice->size(d) / 2)
                n -= lattice->size(d);

            k[d] = n * 2.0 * M_PI / lattice->size(d);
        }

        double kr = k.norm();

        int b = kr * mult;
        if (b >= 0 && b < bins) {
            double ps = 0;
            for (int i = 0; i < n_cmplx; i++) {
                ps += hila::get_complex_element(ftrans[X], i).squarenorm();
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


/// interface for real fields - an extra copy which could be avoided
template <typename T, std::enable_if_t<!hila::contains_complex<T>::value, int> = 0>
std::vector<double> powerspectrum(const Field<T> &f, int bins, double max_k = M_PI) {

    using cmplx_t = Complex<hila::number_type<T>>;

    if constexpr (sizeof(T) % sizeof(Complex<hila::number_type<T>>) == 0) {
        // real field, size is even -- cast the field to pseudo-complex
        // This works because layouts are compatible in all archs - if this changes then
        // need copy
        constexpr int nc = sizeof(T) / sizeof(cmplx_t);

        return powerspectrum(*reinterpret_cast<const Field<Vector<nc, cmplx_t>> *>(&f),
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

        return powerspectrum(cfield, bins, max_k);
    }
}


inline std::vector<double>
powerspectrum_k_info(int bins, double max_k,
                    std::vector<int> *countp = nullptr) {
    // This is for complex type

    double mult = bins / max_k;

    VectorReduction<double> s(bins);
    VectorReduction<int> count(bins);
    s.allreduce(false);
    count.allreduce(false);
    s = 0;
    count = 0;

    onsites (ALL) {

        // mod the k-space around origin, -L/2 -> L/2
        CoordinateVector n;
        Vector<NDIM, double> k;

        foralldir (d) {
            int n = X.coordinate(d);
            if (n > lattice->size(d) / 2)
                n -= lattice->size(d);

            k[d] = n * 2.0 * M_PI / lattice->size(d);
        }

        double kr = k.norm();

        int b = kr * mult;
        if (b >= 0 && b < bins) {
            s[b] += kr;
            count[b] += 1;
        }
    }

    if (countp != nullptr) {
        countp->resize(bins);
        for (int i=0; i<bins; i++) (*countp)[i] = count[i];
    }

    std::vector<double>res(bins);
    for (int i = 0; i < bins; i++) {
        if (count[i] > 0)
            s[i] /= count[i];
        res[i] = s[i];        
    }

    return res;
}


#endif