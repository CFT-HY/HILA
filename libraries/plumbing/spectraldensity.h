#ifndef HILA_SPECTRALDENSITY_H
#define HILA_SPECTRALDENSITY_H

#include "hila.h"


extern hila::timer binning_timer;

/// sd_k_bin_parameters holds the parameters to define binning.

struct sd_k_bin_parameters {
    double max;      // min and max of k (default: 0 and pi)
    double power;    // binned function k^p, default p=1  : note, min and max are "powered"
                     // too
    double binwidth; 
    int bins;        // how many bins in output (default=max of lattice.size(d)/2)
    
    int exact;       // do this many low k-value "bins" exactly NOT IMPLEMENTED YET
    bool bins_set;   // bool to keep track if number of bins have been set

    Vector<NDIM, double> k_scale; // k_i multiplied by this: 1/blocking factor of lattice
};


/// get bin
/// TODO: this should be insider k_binning class but hilapp is not yet able to handle
/// it!

inline double sd_k_norm(const CoordinateVector &cv, const sd_k_bin_parameters &p) {
    Vector<NDIM, double> k = cv.convert_to_k();
    foralldir (d)
        k[d] *= p.k_scale[d];
    return k.norm();
}

inline int sd_get_k_bin(const CoordinateVector &cv, const sd_k_bin_parameters &p) {

    double kr = sd_k_norm(cv, p);

    kr = kr / p.max;
    int b = pow(kr, p.power) * p.bins;
    return b;
}


namespace hila {

/// class hila::k_binning is used to bin "k-space" fields (typically fourier transformed
/// from real space).  Vector k is normalized so that -pi < k_i <= pi, i.e.
//  k_i = 2*pi*x_i/L_i, where L_i = lattice.size(i) and x_i is the L_i/2 modded
/// coordinate.
/// On a blocked lattice (lattice.block()) k is given in units of the base lattice:
///   k_i = 2*pi*x_i/(b_i*L_i), where b_i = (base lattice size)/L_i is the blocking factor.
/// Thus the binning of a blocked field matches the binning of the unblocked one, also
/// for anisotropic blocking.  For a smooth field FFT of the unblocked field
/// ~ (b_0*b_1*...) * FFT of blocked field; spectraldensity() and cross_spectrum() multiply
/// the result by (b_0*b_1*...)^2 so that the normalization matches the unblocked one.
/// The bin_k_field*() -functions bin the given k-space field as is, without rescaling.
/// Bin is determined by formula
///   b = (int) ( (k/k_max)^p * n_bins ), where p is the power in binning - usually 1
///
///   hila::k_binning kb;
///   kb.bins(80).k_max(M_PI);
///   auto sd = kb.spectraldensity(f);  // get the spectral density of field f
///
/// methods:  (use as kb.xxx)
///   hila::k_binning & bins(int n)       set number of bins (k_max const.)
///   int bins()                          get number of bins
///   hila::k_binning & k_max(double m)   set max value of k in binning 
///                  if bins have been modified it is kept constant, otherwise
///                  number of bins is changed so that binwidth is const.
///   double k_max()                      get max value of k
///   hila::k_binning & power(double p)   set the "power" in binning
///   double power()                      get power
///   hila::k_binning & binwidth(double w)  set binwidth (k_max const, so bins changes)
///   double binwidth()                   get binwidth
///
///   std::vector<T>      bin_k_field(const Field<T> &f)       bin the k-space field f
///   std::vector<double> bin_k_field_squarenorm(const Field<T & f)
///                            bin the square norm of k-space field f
///   std::vector<double> spectraldensity(const Field<T> &f)
///                            FFT real-space field f and bin the result in squarenorm
///   std::vector<double> bin_k_field_dot(const Field<T> &a, const Field<T> &b)
///                            bin sum_i Re(a_i^* b_i) of k-space fields a and b
///   std::vector<double> cross_spectrum(const Field<T> &f, const Field<T> &g)
///                            FFT real-space fields f and g and bin sum_i Re(F_i^* G_i)
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
    const lattice_struct *bin_info_lattice = nullptr; // lattice where k_avg was computed

    /// blocking factor of the current lattice relative to the base lattice
    static CoordinateVector blocking_factor() {
        CoordinateVector b;
        foralldir (d)
            b[d] = hila::defined_lattices.front()->l_size[d] / lattice.size(d);
        return b;
    }

    /// rescale FFT-based spectra on blocked lattice by (blocking volume)^2
    static std::vector<double> scale_to_base(std::vector<double> v) {
        auto bf = blocking_factor();
        double bvol = 1;
        foralldir (d)
            bvol *= bf[d];
        if (bvol != 1) {
            for (auto &e : v)
                e *= bvol * bvol;
        }
        return v;
    }

    /// (re)compute bin info if lattice or binning parameters have changed
    void check_bin_info() {
        if (bin_info_lattice != lattice.ptr() || k_avg.size() != par.bins)
            sd_calculate_bin_info();
    }

  public:
    /// Default: binwidth 2*pi/L, where L is the largest base lattice dimension, and
    /// k_max = pi/(minimum blocking factor), i.e. pi on unblocked lattice.
    k_binning() {
        par.power = 1;
        par.exact = 0;
        par.bins = 1;
        par.bins_set = false;
        auto bf = blocking_factor();
        int bmin = bf[0];
        foralldir (d) {
            int L = lattice.size(d) * bf[d];
            if (L > 2 * par.bins)
                par.bins = L / 2;
            bmin = std::min(bmin, bf[d]);
        }
        par.binwidth = M_PI / par.bins;
        par.bins = (par.bins + bmin - 1) / bmin;
        par.max = par.binwidth * par.bins;
    }

    k_binning(int b) {
        par.max = M_PI;
        par.power = 1;
        par.exact = 0;
        par.bins = b;
        par.bins_set = true;
        par.binwidth = par.max / par.bins;
    }

    /// Set number of bins in histogram
    k_binning &bins(int n) {
        assert(n > 0);
        par.bins = n;
        par.bins_set = true;
        par.binwidth = par.max / par.bins;
        k_avg.clear();
        return *this;
    }

    int bins(void) {
        return par.bins;
    }

    /// Max value of k
    k_binning &k_max(double km) {
        assert(km > 0);
        par.max = km;
        if (!par.bins_set) {
            par.bins = ceil(par.max / par.binwidth);
            par.max = par.binwidth * par.bins;
        }
        k_avg.clear();
        return *this;
    }

    double k_max(void) {
        return par.max;
    }

    k_binning &binwidth(double w) {
        assert(w > 0);
        par.binwidth = w;
        par.bins = ceil(par.max / w);
        par.max = par.binwidth * par.bins;
        par.bins_set = false;
        k_avg.clear();
        return *this;
    }

    double binwidth() {
        return par.binwidth;
    }

    /// Bin quantity k^p
    k_binning &power(double p) {
        assert(p > 0);
        par.power = p;
        k_avg.clear();
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

        binning_timer.start();

        check_bin_info();

        ReductionVector<T> s(par.bins);
        s.allreduce(false);
        s = 0;

        onsites(ALL) {

            int b = sd_get_k_bin(X.coordinates(), par);
            if (b >= 0 && b < par.bins) {
                s[b] += f[X];
            }
        }

        binning_timer.stop();

        return s.vector();
    }

    /// sum up the square norm of all elements

    template <typename T>
    std::vector<double> bin_k_field_squarenorm(const Field<T> &f) {

        using float_t = hila::arithmetic_type<T>;
        constexpr int n_float = sizeof(T) / sizeof(float_t);

        binning_timer.start();

        check_bin_info();

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

        binning_timer.stop();

        return s.vector();
    }


    //////////////////////////////////////////////////////////////////////////////////
    /// Spectral density
    /// This version takes in complex field

    template <typename T, std::enable_if_t<hila::contains_complex<T>::value, int> = 0>
    std::vector<double> spectraldensity(const Field<T> &f) {

        Field<T> ftrans;
        FFT_field(f, ftrans);

        return scale_to_base(bin_k_field_squarenorm(ftrans));
    }

    //////////////////////////////////////////////////////////////////////////////////
    /// interface for real fields: use the real Fourier transform, which gives
    /// the same binned result as the complex FFT (see FFT_real_field())
    template <typename T, std::enable_if_t<!hila::contains_complex<T>::value, int> = 0>
    std::vector<double> spectraldensity(const Field<T> &f) {

        Field<T> ftrans;
        FFT_real_field(f, ftrans);

        return scale_to_base(bin_k_field_squarenorm(ftrans));
    }


    /// sum up Re(a_i^* b_i) over all elements of a and b, i.e. the dot product of
    /// the real numbers in a and b (for complex types Re a Re b + Im a Im b).
    /// bin_k_field_dot(f,f) == bin_k_field_squarenorm(f)

    template <typename T>
    std::vector<double> bin_k_field_dot(const Field<T> &fa, const Field<T> &fb) {

        using float_t = hila::arithmetic_type<T>;
        constexpr int n_float = sizeof(T) / sizeof(float_t);

        binning_timer.start();

        check_bin_info();

        ReductionVector<double> s(par.bins);
        s.allreduce(false);
        s = 0;

        onsites(ALL) {

            int b = sd_get_k_bin(X.coordinates(), par);
            if (b >= 0 && b < par.bins) {
                double ps = 0;
                for (int i = 0; i < n_float; i++) {
                    ps += hila::get_number_in_var(fa[X], i) * hila::get_number_in_var(fb[X], i);
                }
                s[b] += ps;
            }
        }

        binning_timer.stop();

        return s.vector();
    }

    //////////////////////////////////////////////////////////////////////////////////
    /// Cross spectrum of fields f and g:  sum over k in a bin of sum_i Re(F_i^*(k) G_i(k)),
    /// where F, G are the Fourier transforms and i runs over the elements of the field
    /// type.  Normalized as spectraldensity() (also on blocked lattice), and cross_spectrum(f,f) == spectraldensity(f).
    /// Real fields use the real Fourier transform (see FFT_real_field()), which gives the
    /// same binned result as the complex FFT.

    template <typename T>
    std::vector<double> cross_spectrum(const Field<T> &f, const Field<T> &g) {

        Field<T> ft, gt;
        if constexpr (hila::contains_complex<T>::value) {
            FFT_field(f, ft);
            FFT_field(g, gt);
        } else {
            FFT_real_field(f, ft);
            FFT_real_field(g, gt);
        }

        return scale_to_base(bin_k_field_dot(ft, gt));
    }


    /// Data structure to hold the binning info - two vectors,
    /// holding average k value in a bin and count of lattice points.
    /// Sets also the k scaling for the current (possibly blocked) lattice.

    void sd_calculate_bin_info() {

        auto bf = blocking_factor();
        foralldir (d)
            par.k_scale[d] = 1.0 / bf[d];
        bin_info_lattice = lattice.ptr();

        k_avg.resize(par.bins);
        bin_count.resize(par.bins);

        ReductionVector<double> s(par.bins);
        ReductionVector<long> count(par.bins);
        s = 0;
        count = 0;

        onsites(ALL) {

            double kr = sd_k_norm(X.coordinates(), par);
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

        check_bin_info();

        if (i >= 0 && i < par.bins)
            return k_avg[i];
        else
            return 0.0;
    }

    /// get the count of points within bin i

    long count(int i) {

        check_bin_info();

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