#ifndef SPECTRALDENSITY_H
#define SPECTRALDENSITY_H

#include "hila.h"


/// sd_k_bin_parameters holds the parameters to define binning.

struct sd_k_bin_parameters {
    double max;   // min and max of k (default: 0 and pi)
    double power; // binned function k^p, default p=1  : note, min and max are "powered"
                  // too
    int bins;     // how many bins in output (default=max of lattice.size(d)/2)
    int exact;    // do this many low k-value "bins" exactly
};


/// get bin
/// TODO: this should be insider k_binning class but hilapp is not yet able to handle
/// it!

inline int sd_get_k_bin(const CoordinateVector &cv, const sd_k_bin_parameters &p) {

    double kr = convert_to_k(cv).norm();

    kr = kr / p.max;
    int b = pow(kr, p.power) * p.bins;
    return b;
}


namespace hila {

/// class hila::k_binning is used to bin "k-space" fields (typically fourier transformed
/// from real space).  Vector k is normalized so that -pi < k_i <= pi, i.e.
//  k_i = 2*pi*x_i/L_i, where L_i = lattice.size(i) and x_i is the L_i/2 modded
/// coordinate.
/// Bin is determined by formula
///   b = (int) ( pow(k/k_max)^p * n_bins )
///
///   hila::k_binning kb;
///   kb.bins(80).k_max(M_PI);
///   auto sd = kb.spectraldensity(f);  // get the spectral density of field f
///
/// methods:  (use as kb.xxx)
///   hila::k_binning & bins(int n)       set number of bins
///   int bins()                          get number of bins
///   hila::k_binning & k_max(double m)   set max value of k in binning
///   double k_max()                      get max value of k
///   hila::k_binning & power(double p)   set the "power" in binning
///   double power()                      get power
///
///   std::vector<T>      bin_k_field(const Field<T> &f)       bin the k-space field f
///   std::vector<double> bin_k_field_squarenorm(const Field<T & f)
///                            bin the square norm of k-space field f
///   std::vector<double> spectraldensity(const Field<T> &f)
///                            FFT real-space field f and bin the result in squarenorm
///
///   double k(int b)                  return the average k within bin b
///   long count(int b)                return the number of points within bin b
///   double bin_min(int b)            return the minimum k of bin b
///   double bin_max(int b)            maximum k in bin b
class k_binning {
  private:
    sd_k_bin_parameters par;
    std::vector<double> k_avg;
    std::vector<size_t> bin_count;

  public:
    k_binning() {
        par.max = M_PI;
        par.power = 1;
        par.exact = 0;
        par.bins = 1;
        foralldir(d) {
            if (lattice.size(d) > 2 * par.bins)
                par.bins = lattice.size(d) / 2;
        }
    }

    k_binning(int b) {
        par.max = M_PI;
        par.power = 1;
        par.exact = 0;
        par.bins = b;
    }

    /// Set number of bins in histogram
    k_binning &bins(int n) {
        assert(n > 0);
        par.bins = n;
        return *this;
    }

    int bins() {
        return par.bins;
    }

    /// Max value of k
    k_binning &k_max(double km) {
        assert(km > 0);
        par.max = km;
        return *this;
    }

    double k_max() {
        return par.max;
    }

    /// Bin quantity k^p
    k_binning &power(double p) {
        assert(p > 0);
        par.power = p;
        return *this;
    }

    double power() {
        return par.power;
    }

    /// Bin exactly this many bins
    k_binning &exact_bins(int e) {
        assert(e >= 0);
        par.exact = e;
        return *this;
    }

    int exact_bins() {
        return par.exact;
    }

    /// Generic k-space field binner routine - returns field element type vector

    template <typename T>
    std::vector<T> bin_k_field(const Field<T> &f) {

        if (k_avg.size() != par.bins)
            sd_calculate_bin_info();

        ReductionVector<T> s(par.bins);
        s.allreduce(false);
        s = 0;

        onsites(ALL) {

            int b = sd_get_k_bin(X.coordinates(), par);
            if (b >= 0 && b < par.bins) {
                s[b] += f[X];
            }
        }

        return s.vector();
    }

    /// sum up the square norm of all elements

    template <typename T>
    std::vector<double> bin_k_field_squarenorm(const Field<T> &f) {

        using float_t = hila::number_type<T>;
        constexpr int n_float = sizeof(T) / sizeof(float_t);

        if (k_avg.size() != par.bins)
            sd_calculate_bin_info();

        ReductionVector<double> s(par.bins);
        s.allreduce(false);
        s = 0;

        onsites(ALL) {

            int b = sd_get_k_bin(X.coordinates(), par);
            if (b >= 0 && b < par.bins) {
                double ps = 0;
                for (int i = 0; i < n_float; i++) {
                    auto a = hila::get_number_in_var(f[X], i);
                    ps += a * a;
                }
                s[b] += ps;
            }
        }

        return s.vector();
    }


    //////////////////////////////////////////////////////////////////////////////////
    /// Spectral density
    /// This version takes in complex field

    template <typename T, std::enable_if_t<hila::contains_complex<T>::value, int> = 0>
    std::vector<double> spectraldensity(const Field<T> &f) {

        Field<T> ftrans;
        FFT_field(f, ftrans);

        return bin_k_field_squarenorm(ftrans);
    }

    //////////////////////////////////////////////////////////////////////////////////
    /// interface for real fields - an extra copy which could be avoided
    template <typename T, std::enable_if_t<!hila::contains_complex<T>::value, int> = 0>
    std::vector<double> spectraldensity(const Field<T> &f) {

        using cmplx_t = Complex<hila::number_type<T>>;

        if constexpr (sizeof(T) % sizeof(Complex<hila::number_type<T>>) == 0) {
            // real field, size is even -- cast the field to pseudo-complex
            // This works because layouts are compatible in all archs - if this changes
            // then need copy
            constexpr int nc = sizeof(T) / sizeof(cmplx_t);

            return spectraldensity(*reinterpret_cast<const Field<Vector<nc, cmplx_t>> *>(&f));
        } else {
            // now the size of input is not evenly divisible by sizeof complex.
            // new complex field
            constexpr int nc = sizeof(T) / sizeof(cmplx_t) + 1;

            Field<Vector<nc, cmplx_t>> cfield;

            onsites(ALL) {

                cfield[X] = 0;
                for (int i = 0; i < sizeof(T) / sizeof(hila::number_type<T>); i++) {
                    auto a = hila::get_number_in_var(f[X], i);
                    hila::set_number_in_var(cfield[X], i, a);
                }
            }

            return spectraldensity(cfield);
        }
    }


    /// Data structure to hold the binning info - two vectors,
    /// holding average k value in a bin and count of lattice points

    void sd_calculate_bin_info() {

        k_avg.resize(par.bins);
        bin_count.resize(par.bins);

        ReductionVector<double> s(par.bins);
        ReductionVector<long> count(par.bins);
        s = 0;
        count = 0;

        onsites(ALL) {

            double kr = convert_to_k(X.coordinates()).norm();
            int b = sd_get_k_bin(X.coordinates(), par);

            if (b >= 0 && b < par.bins) {
                s[b] += kr;
                count[b] += 1;
            }
        }

        for (int i = 0; i < par.bins; i++) {
            bin_count[i] = count[i];
            k_avg[i] = s[i] / count[i];
        }
    }

    /// Get the average k-value of bin i

    double k(int i) {

        if (k_avg.size() == 0) {
            sd_calculate_bin_info();
        }

        if (i >= 0 && i < par.bins)
            return k_avg[i];
        else
            return 0.0;
    }

    /// get the count of points within bin i

    long count(int i) {

        if (bin_count.size() == 0) {
            sd_calculate_bin_info();
        }

        if (i >= 0 && i < par.bins)
            return bin_count[i];
        else
            return 0;
    }

    /// Bin limits
    /// bin is  b = floor((k / k_max)^p * n_bins)
    /// thus, k_min(b) = (b/n_bins)^(1/p) * k_max

    double bin_min(int i) {
        return pow(((double)i) / par.bins, 1.0 / par.power) * par.max;
    }

    double bin_max(int i) {
        return bin_min(i + 1);
    }
};

} // namespace hila

#endif