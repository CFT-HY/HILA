/**
 * @file hila_healthcheck.cpp
 * @author Kari Rummukainen
 * @brief HILA healthcheck application.
 * @details Performs a comprehensive health check on HILA libraries.
 */
#include "hila.h"

#include "clusters.h"

// unistd.h needed for isatty()
#include <unistd.h>


// report the result of the test -- TODO: nicer formatting?
bool report_pass(std::string message, double eps, double limit) {

    static int is_terminal = -1;
    static std::string ok_string, fail_string;

    if (is_terminal < 0) {
        is_terminal = isatty(fileno(stdout));

        if (is_terminal) {
            ok_string = "\x1B[32m --- \033[0m ";
            fail_string = "\x1B[31m *** \033[0m ";
        } else {
            ok_string = " ---  ";
            fail_string = " ***  ";
        }
    }

    if (abs(eps) < limit) {
        hila::out0 << ok_string << message << " passed" << std::endl;
        return true;
    } else {
        hila::out0 << fail_string << message << " FAILED: eps " << eps << " limit " << limit
                   << std::endl;
        return false;
    }
}

//////////////////////////////////////////////////////////////////////////////////////
// hila::global variables in namespaces, used in onsites() and loop functions

namespace hc_globals {
hila::global<double> a;

// unqualified reference inside the namespace
#pragma hila loop_function
inline double times_a(double x) {
    return x * a();
}

namespace inner {
hila::global<double> b;
}

// partially qualified reference
#pragma hila loop_function
inline double get_b() {
    return inner::b();
}
} // namespace hc_globals

namespace {
hila::global<double> hc_anon_global;
}

hila::global<double> hc_top_global;

/**
 * @brief Test hila::global variables in namespaces
 * @details On GPUs globals are __constant__ device variables, whose names hilapp
 * generates from the namespaces of the declaration
 */
void test_globals() {

    // function-local using-directive, hilapp copies it to GPU kernels
    using namespace hc_globals;

    hc_globals::a = 1;
    hc_globals::inner::b = 10;
    hc_anon_global = 100;
    hc_top_global = 1000;

    Field<double> f;
    onsites (ALL) {
        f[X] = hc_globals::a() + times_a(1.0) + hc_globals::inner::b() + get_b() +
               hc_anon_global() + a() + hc_top_global();
    }
    // 1 + 1 + 10 + 10 + 100 + 1 + 1000
    report_pass("hila::global in namespaces", f.sum() / lattice.volume() - 1123, 1e-12);
}

//////////////////////////////////////////////////////////////////////////////////////
// Field of a hila type (template class with base_type and argument_type) declared in a
// namespace: hilapp generates specializations of Field methods for it, which must be
// placed outside the namespace

namespace hc_types {
template <typename T>
struct pt {
    using base_type = hila::arithmetic_type<T>;
    using argument_type = T;
    T x;
};

// Field<pt<double>> used inside the namespace, after the type
inline double neighbour_product_sum() {
    Field<pt<double>> p;
    onsites (ALL) p[X].x = X.x();
    double s = 0;
    onsites (ALL) s += p[X + e_x].x * p[X].x;
    return s;
}
} // namespace hc_types

void test_namespace_field_type() {
    // sum over x of x * ((x+1) mod L), times the number of x-rows
    int L = lattice.size(e_x);
    double expect = 0;
    for (int x = 0; x < L; x++)
        expect += x * ((x + 1) % L);
    expect *= lattice.volume() / L;
    report_pass("Field of a hila type in a namespace",
                (hc_types::neighbour_product_sum() - expect) / expect, 1e-12);
}

/**
 * @brief Test various functions on fields.
 * @details Test functions exp (on both real and complex field) and sin.
 *
 */
void test_functions() {

    Field<double> df = 0;
    report_pass("Field functions: real exp", exp(df).sum() - lattice.volume(), 1e-8);

    Field<Complex<double>> cf = 0;
    report_pass("Field functions: complex exp", abs(exp(cf).sum() - lattice.volume()), 1e-8);

    df[ALL] = sin(X.x() * 2 * M_PI / lattice.size(e_x));
    report_pass("Field functions: sin", df.sum(), 1e-8);
}


/////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test reduction operations on fields
 * @details Test normal reduction `+=` and ReductionVector on fields with different types such as
 * real, complex, matrix and SU matrix.
 *
 */
void test_reductions() {


    // test reductions

    {
        Field<Complex<double>> f;
        f[ALL] = expi(2 * M_PI * X.coordinate(e_x) / lattice.size(e_x));

        Complex<double> sum = 0;
        onsites (ALL)
            sum += f[X];

        sum /= lattice.volume();
        report_pass("Complex reduction value " + hila::prettyprint(sum), abs(sum), 1e-4);

        ReductionVector<Complex<double>> rv(lattice.size(e_x));

        onsites (ALL) {
            rv[X.coordinate(e_x)] += f[X];
        }

        sum = 0;
        for (int i = 0; i < lattice.size(e_x); i++) {
            sum += expi(2 * M_PI * i / lattice.size(e_x)) -
                   rv[i] / (lattice.volume() / lattice.size(e_x));
        }
        report_pass("Complex ReductionVector, sum " + hila::prettyprint(sum), abs(sum), 1e-4);

        // do a combined reduction too
        sum = 0;
        rv = 0;
        onsites (ALL) {
            rv[X.x()] += f[X];
            rv[0] += 1;
            rv[1] += -0.01;

            sum += f[X];
        }

        sum /= lattice.volume();
        Complex<double> sum2 = 0;
        rv[0] -= lattice.volume();
        rv[1] += 0.01 * lattice.volume();

        for (int i = 0; i < lattice.size(e_x); i++) {
            sum2 += expi(2 * M_PI * i / lattice.size(e_x)) -
                    rv[i] / (lattice.volume() / lattice.size(e_x));
        }
        report_pass("Combined reductions, sum " + hila::prettyprint(sum) + ", sum2 " +
                        hila::prettyprint(sum2),
                    abs(sum) + abs(sum2), 1e-4);


        f[ALL] = 1;
        f[ODD] = -1;

        sum = 0;
        onsites (EVEN)
            sum += f[X];
        report_pass("Even site reduction", abs(sum - lattice.volume() / 2), 1e-2);

        sum = 0;
        onsites (ODD)
            sum += f[X];
        report_pass("Odd site reduction", abs(sum + lattice.volume() / 2), 1e-2);
    }

    {
        // reductionvector with long
        Field<int64_t> lf;
        lf[ALL] = X.x();

        ReductionVector<int64_t> rv(lattice.size(e_x));

        onsites (ALL) {
            rv[X.x()] += lf[X];
        }

        long s = 0;
        for (int x = 0; x < rv.size(); x++) {
            s += abs(rv[x] - x * (lattice.volume() / lattice.size(e_x)));
        }

        report_pass("ReductionVector<int64_t>, sum " + hila::prettyprint(s), s, 1e-15);

        int64_t psum = 0;
        onsites (ALL) {
            if (X.parity() == EVEN)
                psum += 1;
        }
        // If volume is odd there is 1 extra even site
        psum -= lattice.volume() / 2 + lattice.volume() % 2;

        report_pass("Reduction with parity", psum, 1e-15);
    }

    {
        // reductionvector with long
        Field<int64_t> lf;
        lf[ALL] = X.x();

        ReductionVector<int64_t> rv(lattice.size(e_x));

        onsites (ALL) {
            rv[X.x()] += lf[X];
        }

        long s = 0;
        for (int x = 0; x < rv.size(); x++) {
            s += abs(rv[x] - x * (lattice.volume() / lattice.size(e_x)));
        }

        report_pass("ReductionVector<int64_t>, sum " + hila::prettyprint(s), s, 1e-15);
    }

#if NDIM > 2
    {

        constexpr int N = 5;

        Field<SU<N, double>> mf;
        mf = 1;
        // onsites(ALL) mf[X] = (mf[X] + mf[X])*0.5;

        ReductionVector<SU<N, double>> rmf(lattice.size(e_x));

        onsites (ALL) {
            rmf[X.x()] += mf[X];
        }


        double diff = 0;
        for (int i = 0; i < lattice.size(e_x); i++) {
            diff += (rmf[i] - lattice.size(e_y) * lattice.size(e_z)).squarenorm();
        }


        report_pass("SU(" + std::to_string(N) + ") ReductionVector", diff, 1e-8);
    }
#endif
}


/**
 * @brief Test site access
 * @details Test simple site access by writing and reading to random field site index c.
 *
 */
void test_site_access() {

    Field<Complex<double>> f = 0;
    Complex<double> sum;

    CoordinateVector c;
    for (int i = 0; i < 3; i++) {
        foralldir (d)
            c[d] = hila::random() * lattice.size(d);
        hila::broadcast(c); // need to broadcast!

        f[c] = 4;   // write to location c
        sum = f[c]; // read back - does mpi comm
        report_pass("Setting and reading a value at " + hila::prettyprint(c.transpose()),
                    abs(sum - 4), 1e-10);
    }
}


/**
 * @brief Test min and max operations on fields
 * @details Test min and max operations on fields with different parities.
 *
 */
void test_minmax() {

    CoordinateVector c, loc;
    Field<double> n;

    for (Parity par : {ALL, EVEN, ODD}) {

        n[par] = hila::random();
        do {
            foralldir (d)
                c[d] = hila::random() * lattice.size(d);
        } while (!(par == ALL || c.parity() == par));
        hila::broadcast(c);
        n[c] = 2;

        if (par != ALL)
            n[opp_parity(par)] = 3;

        auto v = n.max(par, loc);
        report_pass("Maxloc parity " + hila::prettyprint(par) + " is " +
                        hila::prettyprint(loc.transpose()),
                    (c - loc).norm(), 1e-8);
        report_pass("Max parity " + hila::prettyprint(par) + " value " + hila::prettyprint(v),
                    v - 2, 1e-9);

        do {
            foralldir (d)
                c[d] = hila::random() * lattice.size(d);
        } while (!(par == ALL || c.parity() == par));
        hila::broadcast(c);
        n[c] = -1;

        if (par != ALL)
            n[opp_parity(par)] = -3;

        v = n.min(par, loc);
        report_pass("Minloc parity " + hila::prettyprint(par) + " is " +
                        hila::prettyprint(loc.transpose()),
                    (c - loc).norm(), 1e-8);
        report_pass("Min parity " + hila::prettyprint(par) + " value " + hila::prettyprint(v),
                    v + 1, 1e-9);
    }
}

/**
 * @brief Test Gaussian random field generation
 * @details Generate a Gaussian random field and check its average and width^2.
 *
 */
void test_random() {

    constexpr int n_loops = 100;

    Field<double> f;

    double fsum = 0, fsqr = 0;
    for (int i = 0; i < n_loops; i++) {
        f.gaussian_random();
        double s = 0, s2 = 0;
        onsites (ALL) {
            s += f[X];
            s2 += sqr(f[X]);
        }
        fsum += s / lattice.volume();
        fsqr += s2 / lattice.volume();
    }

    fsum /= n_loops;
    fsqr /= n_loops;

    report_pass("Gaussian random average (6 sigma limit) " + hila::prettyprint(fsum), abs(fsum),
                6 / sqrt(((double)n_loops) * lattice.volume()));

    report_pass("Gaussian random width^2 " + hila::prettyprint(fsqr), fsqr - 1,
                6 / sqrt(((double)n_loops) * lattice.volume()));
}


/**
 * @brief Test setting elements and selecting them
 * @details Set elements in a field at random coordinates and then select them using SiteSelect
 * and SiteValueSelect.
 */
void test_set_elements_and_select() {

    Field<Complex<double>> f = 0;
    CoordinateVector c;

    std::vector<CoordinateVector> cvec;
    std::vector<Complex<double>> vals;


    if_rank0 () {
        int k = 1;
        for (int i = 0; i <= 50; i++) {
            foralldir (d)
                c[d] = hila::random() * lattice.size(d);
            bool found = false;

            // don't overwrite previous loc
            for (auto &r : cvec)
                if (r == c)
                    found = true;
            if (!found) {
                cvec.push_back(c);
                vals.push_back(Complex<double>(k, k));
                k++;
            }
        }
    }

    hila::broadcast(vals);
    hila::broadcast(cvec);

    f.set_elements(vals, cvec);

    vals = f.get_elements(cvec, true);

    Complex<double> sum = 0;
    for (int i = 0; i < vals.size(); i++) {
        sum += vals[i] - Complex<double>(i + 1, i + 1);
    }

    report_pass("Field set_elements and get_elements with " + hila::prettyprint(vals.size()) +
                    " coordinates",
                abs(sum), 1e-8);

    // use the same vector for siteselect

    SiteSelect s;
    SiteValueSelect<Complex<double>> sv;

    onsites (ALL) {
        if (squarenorm(f[X]) >= 1) {
            s.select(X);
            sv.select(X, f[X]);
        }
    }

    report_pass("SiteSelect size " + hila::prettyprint(s.size()), s.size() - cvec.size(), 1e-3);
    report_pass("SiteValueSelect size " + hila::prettyprint(sv.size()), sv.size() - cvec.size(),
                1e-3);

    int ok = 1;

    for (auto &c : cvec) {
        bool found = false;
        for (int i = 0; !found && i < s.size(); i++) {
            if (s.coordinates(i) == c)
                found = true;
        }
        if (!found)
            ok = 0;
    }

    report_pass("SiteSelect content", 1 - ok, 1e-6);

    ok = 1;
    for (int k = 0; k < cvec.size(); k++) {

        bool found = false;
        for (int i = 0; !found && i < sv.size(); i++) {
            if (sv.coordinates(i) == cvec[k] && sv.value(i) == vals[k])
                found = true;
        }
        if (!found)
            ok = 0;
    }

    report_pass("SiteValueSelect content", 1 - ok, 1e-6);


    onsites (ALL) {}
}

/**
 * @brief Shifts
 */

void test_shift() {

#if NDIM > 1

    Field<CoordinateVector> cf;

    onsites (ALL) {
        cf[X] = X.coordinates();
    }

    int64_t sumvx = 0, sumvxy = 0;
    onsites (ALL) {
        sumvx += (cf[X + 3 * e_x] - cf[X]).mod(lattice.size()).norm_L1() - 3;
        sumvxy += (cf[X + e_x] - cf[X - 2 * e_y]).mod(lattice.size()).norm_L1() - 3;
    }

    report_pass("Shift fields", abs(sumvx) + abs(sumvxy), 1e-6);

#endif
}


/**
 * @brief Test subvolume operations
 *
 */
void test_subvolumes() {

    bool ok = true;
    for (int i = 0; i < 20; i++) {
        CoordinateVector c;
        foralldir (d)
            c[d] = hila::random() * lattice.size(d);
        auto si = SiteIndex(c);
        if (si.coordinates() != c)
            ok = false;
    }


    report_pass("SiteIndex", ok == false, 1e-2);


    Field<SiteIndex> f;
    f[ALL] = SiteIndex(X.coordinates());

    CoordinateVector c;
    c.fill(-1);
    size_t vol = lattice.volume();
    foralldir (d)
        if (d < NDIM - 1) {
            c[d] = hila::random() * lattice.size(d);
            hila::broadcast(c);

            auto slice = f.get_slice(c);

            vol /= lattice.size(d);
            if_rank0 () {

                report_pass(hila::prettyprint(NDIM - (int)d - 1) + "-dimensional slice size " +
                                hila::prettyprint(vol),
                            slice.size() - vol, 1e-3);


                bool pass = true;

                CoordinateVector mycoord = 0;
                foralldir (d2) {
                    if (c[d2] >= 0)
                        mycoord[d2] = c[d2];
                }
                foralldir (d2) {
                    if (c[d2] < 0) {
                        mycoord[d2] = -1;
                        break;
                    }
                }

                for (auto s : slice) {

                    // get the coordinate which should be here
                    bool add = true;
                    foralldir (d2) {
                        if (add && c[d2] < 0) {
                            mycoord[d2]++;
                            if (mycoord[d2] < lattice.size(d2)) {
                                add = false;
                            } else {
                                mycoord[d2] = 0;
                            }
                        }
                    }

                    if (pass && mycoord != s.coordinates()) {
                        hila::out0 << "Slice coord error, should be "
                                   << hila::prettyprint(mycoord.transpose()) << " is "
                                   << hila::prettyprint(s.coordinates().transpose()) << " slice is "
                                   << hila::prettyprint(c.transpose()) << '\n';
                        pass = false;
                        break;
                    }
                }

                report_pass("slice content", pass == false, 1e-2);
            }
        }
}

/**
 * @brief Test FFT
 * @details Test forward and inverse FFT on fields.
 *
 */


void test_fft() {


    {
        Field<Complex<double>> f, p, p2;

        // Start with unit field
        f = 1;

        // After one FFT the field is 0 except at coord 0
        p2 = 0;

        p2[CoordinateVector(0)] = lattice.volume();

        FFT_field(f, p);

        double eps = squarenorm_relative(p, p2);

        report_pass("FFT Complex<double> constant field", eps, 1e-13 * sqrt(lattice.volume()));

        //-----------------------------------------------------------------
        // After two applications the field should be back to a constant * volume

        FFT_field(p, f, fft_direction::back);

        double sum = 0;
        double tnorm = 0;
        onsites (ALL) {
            sum += (f[X] - lattice.volume()).squarenorm();
            tnorm += f[X].squarenorm();
        }

        eps = fabs(sum / tnorm);

        report_pass("FFT Complex<double> inverse transform", eps, 1e-10);
    }

    {
        Field<Complex<float>> f, p, p2;

        // Start with unit field
        f = 1;

        // After one FFT the field is 0 except at coord 0
        p2 = 0;

        p2[CoordinateVector(0)] = lattice.volume();

        FFT_field(f, p);

        double eps = squarenorm_relative(p, p2);

        // float rounding is larger for lattice sizes with large prime factors (e.g. 122 = 2*61)
        report_pass("FFT Complex<float> constant field", eps, 5e-6 * sqrt(lattice.volume()));

        //-----------------------------------------------------------------
        // After two applications the field should be back to a constant * volume

        FFT_field(p, f, fft_direction::back);

        double sum = 0;
        double tnorm = 0;
        onsites (ALL) {
            sum += (f[X] - lattice.volume()).squarenorm();
            tnorm += f[X].squarenorm();
        }

        eps = fabs(sum / tnorm);

        report_pass("FFT Complex<float> inverse transform", eps, 1e-10);
    }


    //-----------------------------------------------------------------

    {

        Field<Complex<double>> f, p, p2;
        double eps;

        for (int iter = 0; iter < 5; iter++) {
            Vector<NDIM, double> kv;
            CoordinateVector kx;
            foralldir (d) {
                kx[d] = hila::broadcast(hila::random()) * lattice.size(d);
            }

            kv = kx.convert_to_k();


            onsites (ALL) {
                double d = kv.dot(X.coordinates());
                f[X] = expi(d);
            }

            FFT_field(f, p);

            p2 = 0;
            p2[kx] = lattice.volume();

            eps = squarenorm_relative(p, p2);

            report_pass("FFT of wave vector " + hila::prettyprint(kx.transpose()), eps,
                        1e-13 * sqrt(lattice.volume()));
        }

        //-----------------------------------------------------------------

        {
            Field<double> r;
            onsites (ALL)
                r[X] = hila::gaussrand();

            f = r.FFT_real_to_complex();
            p = f.FFT(fft_direction::back) / lattice.volume();
            eps = squarenorm_relative(r, p);

            report_pass("FFT real to complex", eps, 1e-13 * sqrt(lattice.volume()));

            bool is_odd = false;
            foralldir (d)
                is_odd |= (lattice.size(d) % 2 > 0);
            if (!is_odd) {
                auto r2 = f.FFT_complex_to_real(fft_direction::back) / lattice.volume();
                eps = squarenorm_relative(r, r2);

                report_pass("FFT complex to real", eps, 1e-13 * sqrt(lattice.volume()));
            } else {
                hila::out0 << " ...  Skipping FFT complex to real because lattice size is odd\n";
            }
        }

        //-----------------------------------------------------------------
        // Check fft norm


        onsites (ALL) {
            p[X] = hila::random() * exp(-X.coordinates().convert_to_k().squarenorm());
        }
        double np = p.squarenorm();

        bool is_odd = false;
        foralldir (d)
            is_odd |= (lattice.size(d) % 2 > 0);
        if (!is_odd) {
            f = p.FFT(fft_direction::back) / sqrt(lattice.volume());

            double nf = f.squarenorm();
            report_pass("Norm of field = " + hila::prettyprint(nf) +
                            " and FFT = " + hila::prettyprint(np),
                        (nf - np) / nf, 1e-10);
        } else {
            hila::out0 << " ...  Skipping FFT norm check because lattice size is odd\n";
        }


        hila::k_binning b;
        b.k_max(M_PI * sqrt(1.0 * NDIM));

        auto bf = b.bin_k_field(p.conj() * p);

        double s = 0;
        for (auto b : bf) {
            s += abs(b);
        }
        hila::broadcast(s);

        report_pass("Norm of binned FFT = " + hila::prettyprint(s), (s - np) / np, 1e-10);
    }
}

/**
 * @brief Test real Fourier transform FFT_real
 * @details Compare with the complex FFT: element values to one direction, binned power
 * and cross spectra, Parseval, inverse transform, float and multi-component fields.
 */
void test_fft_real() {

    // relative difference of two binned spectra
    auto rel_diff = [](const std::vector<double> &a, const std::vector<double> &b) {
        double d = 0, n = 0;
        for (int i = 0; i < a.size(); i++) {
            d += fabs(a[i] - b[i]);
            n += fabs(b[i]);
        }
        return hila::broadcast(d / n);
    };

    Field<double> r, g;
    onsites (ALL) {
        r[X] = hila::gaussrand();
        g[X] = 0.5 * r[X] + hila::gaussrand();
    }
    double rnorm = r.squarenorm();

    // transform to one direction at a time, compare elementwise with the complex transform:
    // h(k) = sqrt(2) Re F(k) for 0 < k < L/2,  h(L-k) = sqrt(2) Im F(k) = -sqrt(2) Im F(L-k)
    foralldir (d) {
        CoordinateVector dirs;
        dirs.fill(false);
        dirs[d] = true;

        Field<Complex<double>> F;
        onsites (ALL)
            F[X] = Complex<double>(r[X], 0);
        F = F.FFT(dirs);

        Field<double> h = r.FFT_real(dirs);
        int L = lattice.size(d);
        double eps = 0;
        onsites (ALL) {
            int k = X.coordinate(d);
            double expect;
            if (k == 0 || 2 * k == L)
                expect = F[X].real();
            else if (2 * k < L)
                expect = M_SQRT2 * F[X].real();
            else
                expect = -M_SQRT2 * F[X].imag();
            eps += squarenorm(h[X] - expect);
        }
        report_pass("FFT_real to direction " + hila::prettyprint(d) + " vs. complex FFT",
                    eps / rnorm, 1e-26 * L);
    }

    Field<double> hr = r.FFT_real();
    Field<double> hg = g.FFT_real();

    {
        // all directions at once vs. one direction at a time
        Field<double> hs = r;
        foralldir (d) {
            CoordinateVector dirs;
            dirs.fill(false);
            dirs[d] = true;
            hs = hs.FFT_real(dirs);
        }
        report_pass("FFT_real all directions vs. one at a time",
                    squarenorm_relative(hr, hs) / hr.squarenorm(), 1e-26);
    }

    // Parseval
    report_pass("FFT_real Parseval", fabs(hr.squarenorm() / lattice.volume() - rnorm) / rnorm,
                1e-12);

    // inverse
    Field<double> r2 = hr.FFT_real(fft_direction::back) / lattice.volume();
    report_pass("FFT_real inverse transform", squarenorm_relative(r, r2) / rnorm, 1e-26);

    // binned power and cross spectra vs complex FFT
    {
        auto Fr = r.FFT_real_to_complex();
        auto Fg = g.FFT_real_to_complex();

        hila::k_binning b;
        b.k_max(M_PI * sqrt(1.0 * NDIM));

        report_pass("FFT_real binned power spectrum",
                    rel_diff(b.bin_k_field_squarenorm(hr), b.bin_k_field_squarenorm(Fr)), 1e-12);

        Field<double> prod_real = hr * hg;
        Field<double> prod_cmplx;
        onsites (ALL)
            prod_cmplx[X] = (Fr[X].conj() * Fg[X]).real();

        report_pass("FFT_real binned cross spectrum",
                    rel_diff(b.bin_k_field(prod_real), b.bin_k_field(prod_cmplx)), 1e-12);

        // spectraldensity() uses FFT_real for real fields
        Field<Complex<double>> rc;
        onsites (ALL)
            rc[X] = Complex<double>(r[X], 0);
        report_pass("spectraldensity of real field",
                    rel_diff(b.spectraldensity(r), b.spectraldensity(rc)), 1e-12);
    }

    // multi-component field: each component transformed separately
    {
        Field<Vector<2, double>> v;
        onsites (ALL) {
            v[X].e(0) = r[X];
            v[X].e(1) = g[X];
        }
        auto hv = v.FFT_real();
        double eps = 0;
        onsites (ALL) {
            eps += squarenorm(hv[X].e(0) - hr[X]) + squarenorm(hv[X].e(1) - hg[X]);
        }
        report_pass("FFT_real Vector<2,double>", eps / hr.squarenorm(), 1e-26);
    }

    // float
    {
        Field<float> rf;
        onsites (ALL)
            rf[X] = r[X];
        Field<float> hf = rf.FFT_real();
        double eps = 0;
        onsites (ALL)
            eps += squarenorm(hf[X] - hr[X]);
        report_pass("FFT_real float vs double", eps / hr.squarenorm(), 1e-11);

        Field<float> rf2 = hf.FFT_real(fft_direction::back) / lattice.volume();
        report_pass("FFT_real float inverse", squarenorm_relative(rf, rf2) / rnorm, 1e-11);
    }

}

/**
 * @brief Test spectral density
 * @details Test spectral density extraction from field.
 *
 */
void test_spectraldensity() {


    // test spectral density for single waves

    Field<Complex<double>> f, p;

    hila::k_binning b;
    b.k_max(M_PI * sqrt(1.0 * NDIM));
    b.bins(b.bins() - 1); // reduce the number of default bins by 1

    // test std binning first

    for (int iter = 0; iter < 3; iter++) {
        Vector<NDIM, double> kv;
        CoordinateVector kx;
        foralldir (d) {
            kx[d] = hila::broadcast(hila::random()) * lattice.size(d);
        }
        kv = kx.convert_to_k();
        auto absk = kv.norm();

        // test first std. binning (normally )
        f = 0;
        f[kx] = 1;

        auto fb = b.bin_k_field(f);

        double sum = 0;
        for (int i = 0; i < b.bins(); i++) {
            if (b.bin_min(i) <= absk && b.bin_max(i) > absk) {
                // should be one hit here
                sum += squarenorm(fb[i] - 1);
            } else {
                sum += squarenorm(fb[i]);
            }
        }

        report_pass("Binning test at vector " + hila::prettyprint(kx.transpose()), sum, 1e-10);

        // then test the spectral density extraction
        // set single wave
        onsites (ALL) {
            double d = kv.dot(X.coordinates());
            f[X] = expi(d);
        }

        auto fsd = b.spectraldensity(f);

        sum = 0;
        for (int i = 0; i < b.bins(); i++) {
            if (b.bin_min(i) <= absk && b.bin_max(i) > absk) {
                // should be one hit here
                sum += fabs(fsd[i] / pow(lattice.volume(), 2) - 1);
            } else {
                sum += fabs(fsd[i]);
            }
        }

        report_pass("Spectral density test with above vector ", sum, 1e-10);
    }
}

/**
 * @brief Test cross spectrum
 * @details Compare k_binning::cross_spectrum with Re(F^* G) binned from complex FFTs,
 * for real, complex and Vector fields
 */
void test_cross_spectrum() {

    auto rel_diff = [](const std::vector<double> &a, const std::vector<double> &b) {
        double d = 0, n = 0;
        for (int i = 0; i < a.size(); i++) {
            d += fabs(a[i] - b[i]);
            n += fabs(b[i]);
        }
        return hila::broadcast(d / n);
    };

    hila::k_binning b;
    b.k_max(M_PI * sqrt(1.0 * NDIM));

    // real fields, correlated
    Field<double> r, s;
    onsites (ALL) {
        r[X] = hila::gaussrand();
        s[X] = 0.5 * r[X] + hila::gaussrand();
    }

    // reference: Re(F^* G) from complex FFTs
    auto Fr = r.FFT_real_to_complex();
    auto Fs = s.FFT_real_to_complex();
    Field<double> ref;
    onsites (ALL)
        ref[X] = (Fr[X].conj() * Fs[X]).real();

    report_pass("cross_spectrum of real fields",
                rel_diff(b.cross_spectrum(r, s), b.bin_k_field(ref)), 1e-12);

    report_pass("cross_spectrum(f,f) == spectraldensity(f)",
                rel_diff(b.cross_spectrum(r, r), b.spectraldensity(r)), 1e-12);

    // complex fields
    Field<Complex<double>> cr, cs;
    onsites (ALL) {
        cr[X] = Complex<double>(hila::gaussrand(), hila::gaussrand());
        cs[X] = 0.5 * cr[X] + Complex<double>(hila::gaussrand(), hila::gaussrand());
    }
    auto Fcr = cr.FFT();
    auto Fcs = cs.FFT();
    onsites (ALL)
        ref[X] = (Fcr[X].conj() * Fcs[X]).real();

    report_pass("cross_spectrum of complex fields",
                rel_diff(b.cross_spectrum(cr, cs), b.bin_k_field(ref)), 1e-12);

    // Vector field: sum over components
    Field<Vector<2, double>> v, w;
    onsites (ALL) {
        v[X].e(0) = r[X];
        v[X].e(1) = s[X];
        w[X].e(0) = s[X];
        w[X].e(1) = 2 * r[X];
    }
    auto c_vw = b.cross_spectrum(v, w);
    auto c_rs = b.cross_spectrum(r, s);
    auto c_sr = b.cross_spectrum(s, r);
    std::vector<double> c_sum(c_rs.size());
    for (int i = 0; i < c_sum.size(); i++)
        c_sum[i] = c_rs[i] + 2 * c_sr[i];

    report_pass("cross_spectrum of Vector<2,double> fields", rel_diff(c_vw, c_sum), 1e-12);
}

//--------------------------------------------------------------------------------

void test_field_slices() {

    Field<SiteIndex> s;
    onsites (ALL) {
        s[X] = SiteIndex(X.coordinates());
    }
}


/**
 * @brief Test matrix operations
 * @details Test matrix multiplication and addition, as well as operations with imaginary operator.
 *
 */
void test_matrix_operations() {

    Field<Matrix<3, 2, Complex<double>>> mf;

    onsites (ALL)
        mf[X].fill(1 + I);

    Matrix<3, 3, Complex<double>> cm;
    cm.asArray() = 4;
    double sum = 0;
    onsites (ALL) {
        sum += (mf[X] * mf[X].dagger() - cm).squarenorm();
    }

    report_pass("matrix multiply and addition", sum, 1e-8);

    auto dm = cm * I - 2 * I;
    dm.asArray() *= I;
    dm = ((dm - 2).asArray() + 4).asMatrix();
    report_pass("Array and imaginary unit operations", dm.squarenorm(), 1e-8);
}


/**
 * @brief Test matrix decomposition
 * @details Test matrix decompositions such as eigen decomposition and singular value decomposition
 * (SVD).
 */
void test_matrix_algebra() {

    using myMatrix = SquareMatrix<4, Complex<double>>;

    Field<myMatrix> M;
    Field<double> delta;

    M.gaussian_random(2.0);

    // eigenvalue test - show that  M = U D U^*, where D is diagonal eigenvalue matrix and U
    // matrix of eigenvectors

    onsites (ALL) {
        auto H = M[X] * M[X].dagger(); // make hermitean

        auto r = H.eigen_hermitean();
        delta[X] = (H - r.eigenvectors * r.eigenvalues * r.eigenvectors.dagger()).norm();
    }

    auto max_delta = delta.max();

    report_pass("Eigenvalue analysis with " + hila::prettyprint(myMatrix::rows()) + "x" +
                    hila::prettyprint(myMatrix::columns()) + " Hermitean matrix",
                max_delta, 1e-10);

    // Singular value test - non-pivoted

    onsites (ALL) {
        auto r = M[X].svd();
        delta[X] = (M[X] - r.U * r.singularvalues * r.V.dagger()).norm();
    }

    max_delta = delta.max();

    report_pass("SVD with " + hila::prettyprint(myMatrix::rows()) + "x" +
                    hila::prettyprint(myMatrix::columns()) + " Complex matrix",
                max_delta, 1e-10);

    // pivoted singular values

    M.gaussian_random();

    onsites (ALL) {
        auto r = M[X].svd_pivot(hila::sort::ascending);
        delta[X] = (M[X] - r.U * r.singularvalues * r.V.dagger()).norm();
    }

    max_delta = delta.max();

    report_pass("Fully pivoted SVD with " + hila::prettyprint(myMatrix::rows()) + "x" +
                    hila::prettyprint(myMatrix::columns()) + " Complex matrix",
                max_delta, 1e-10);


    M += 2; // ensure that M is invertible
    // one
    auto one = myMatrix(1);
    onsites (ALL) {
        auto inv = M[X].LU_solve(one);
        delta[X] = (M[X] * inv - one).norm();
    }
    max_delta = delta.max();

    report_pass("LU Inversion of " + hila::prettyprint(myMatrix::rows()) + "x" +
                    hila::prettyprint(myMatrix::columns()) + " Complex matrix",
                max_delta, 1e-10);
}

/**
 * @brief Test matrix/vector elementwise operations
 *
 */

void test_element_operations() {

    Field<Matrix<3, 2, Complex<double>>> mf;

    mf = 0;

    double sum = 0, sum2 = 0, sum3 = 0;
    onsites (ALL) {
        mf[X] = hila::elem::exp(mf[X]);
        sum += (hila::elem::sub(mf[X], 1)).squarenorm();
        sum2 += hila::elem::log(mf[X]).squarenorm();
        sum3 += (hila::elem::mul(hila::elem::add(mf[X], -2 * mf[X]), -1) - mf[X]).squarenorm();
    }

    report_pass("hila::elem::exp, log, sub, add, mul", sum + sum2 + sum3, 1e-8);
}


/**
 * @brief Test extended type
 * @details Test extended type for sums that exibit loss in accuracy with double.
 *
 */
void test_extended() {

    // ExtendedPrecision type test with T = double
    Field<double> g;

    ExtendedPrecision e = 2.4;
    e = 5 * e - e - e - 3 * e;
    report_pass("ExtendedPrecision basic arithmetics: " + hila::prettyprint(e), fabs(e.to_double()),
                1e-20);

    e = 1;
    e += 1e-25;
    bool ok = (e > 1) && (e == e);

    report_pass("ExtendedPrecision comparison ops", ok ? 0 : 1, 1e-20);

    double power = 150.0 / lattice.volume();

    onsites (ALL) {
        auto i = SiteIndex(X.coordinates()).value;
        g[X] = exp(-(i * power));
    }

    ExtendedPrecision ev = 0;
    double s = 0;
    float f = 0;

    for (int i = 0; i < 100; i++) {
        onsites (ALL) {
            ev += g[X];
            s += g[X];
            f += g[X];
        }
    }

    double result = ev.to_double();

    for (int i = 0; i < 100; i++) {
        onsites (ALL) {
            ev += -g[X];
            s += -g[X];
            f += -g[X];
        }
    }

    // hila::out0 << ev.value << ' ' << ev.compensation << ' ' << ev.compensation2 << '\n';

    std::stringstream ssres;
    ssres << "ExtendedPrecision reduction: residual " << ev.to_double() / result << " (double "
          << s / result << ", float " << f / result << ")";

    report_pass(ssres.str(), ev.to_double() / result, 1e-25);
}


//--------------------------------------------------------------------------------

void test_clusters() {

    foralldir (d) {
        if (lattice.size(d) == 1) {
            hila::out0 << " ...  Skipping cluster test because lattice extent is 1\n";
            return;
        }
    }

    {
#if NDIM > 2 && !defined(VECTORIZED)

        Field<int> m;

        m = hila::Clusters::background;

        onsites (ALL) {
            auto c = X.coordinates();
            if (c[e_x] == 0 && c[e_y] == 0)
                m[X] = 1;
            else if (c[e_x] == 2 && c[e_z] == 2)
                m[X] = 2;
            else if (c[e_x] == 4 && c[e_y] < 3 && c[e_z] > 1 && c[e_z] <= 4)
                m[X] = 3;
        }

        hila::Clusters cl(m);

        report_pass("Cluster test: number of clusters ", cl.number() - 3, 1e-10);
        if (cl.number() == 3) {
            double sumsize =
                fabs(cl[0].size() - (lattice.volume() / (lattice.size(e_x) * lattice.size(e_y)))) +
                fabs(cl[1].size() - (lattice.volume() / (lattice.size(e_x) * lattice.size(e_z)))) +
                fabs(cl[2].size() - 1 * 3 * 3);
            report_pass("Cluster test: cluster sizes ", sumsize, 1e-10);

            double types = abs(cl[0].type() - 1) + abs(cl[1].type() - 2) + abs(cl[2].type() - 3);
            report_pass("Cluster test: cluster types ", types, 1e-10);

#if NDIM == 3
            double area = abs(cl[0].area() - 4 * lattice.size(e_z)) +
                          abs(cl[1].area() - 4 * lattice.size(e_y)) +
                          abs(cl[2].area() - 2 * (1 * 3 + 1 * 3 + 3 * 3));

            report_pass("Cluster test: cluster area ", area, 1e-10);

#endif
        }

        // test coordinates of cluster
        m = hila::Clusters::background;
        CoordinateVector cv = 0;
        for (int x = 0; x < 4; x++) {
            cv[e_x] = x;
            m[cv] = 1;
        }
        cl.find(m);

        int ok = (cl.number() == 1);
        if (ok) {
            auto c = cl[0].sites();
            ok = (c.size() == 4);

            // the coordinates are in unspecified order
            for (int i = 0; ok && i < 4; i++) {
                auto cvc = c[i].coordinates();
                bool found = false;
                for (int x = 0; !found && x < 4; x++) {
                    cv[e_x] = x;
                    found = (cv == cvc);
                }
                ok = found;
            }
        }
        report_pass("Cluster test: cluster sites ", ok - 1, 1e-10);
    }
    {
#if NDIM == 3
        VectorField<uint8_t> link = 0;
        onsites (ALL) {
            auto c = X.coordinates();
            if ((c[e_x] == 0 || c[e_x] == lattice.size(e_x) / 2) &&
                (c[e_y] == 0 || c[e_y] == lattice.size(e_y) / 2) &&
                c[e_z] < lattice.size(e_z) / 2) {
                link[e_z][X] = 1;
            }

            if ((c[e_x] == 0 || c[e_x] == lattice.size(e_x) / 2) &&
                c[e_y] < lattice.size(e_y) / 2 &&
                (c[e_z] == 0 || c[e_z] == lattice.size(e_z) / 2)) {
                link[e_y][X] = 1;
            }
        }

        hila::Clusters cl(link);

        report_pass("Link cluster test: number of clusters ", cl.number() - 2, 1e-10);

        if (cl.number() == 2) {
            auto s = 2 * (lattice.size(e_y) / 2 + lattice.size(e_z) / 2);
            s = abs(cl[0].size() - s) + abs(cl[1].size() - s);
            report_pass("Link cluster sizes", s, 1e-10);
        }

#endif

#endif
    }
}
//--------------------------------------------------------------------------------

/**
 * @brief Test k_binning on an anisotropically blocked lattice
 * @details Field cos(2 pi n x_d / L_d) g(x_perp), blocked by 2 to direction d, has the same
 * Fourier modes as the unblocked one.  Binned spectra must agree if k is scaled to base
 * lattice units and the spectrum by (blocking volume)^2.
 */
void test_spectraldensity_blocking() {

    Direction db = NDIRECTIONS;
    foralldir (d) {
        if (db == NDIRECTIONS && lattice.size(d) % 2 == 0 && lattice.size(d) >= 8)
            db = d;
    }
    if (db == NDIRECTIONS || !lattice.is_base())
        return;

    CoordinateVector blocking;
    blocking.fill(1);
    blocking[db] = 2;

    hila::print_dashed_line();
    hila::out0 << "Testing spectraldensity on lattice blocked by " << blocking << '\n';

    const int n = 3;
    const int L = lattice.size(db);
    Field<double> f, fb;
    onsites (ALL) {
        double g = 1;
        foralldir (d) {
            if (d != db) {
                int x = X.coordinate(d);
                g *= 1.0 + sin(0.7 * x * x + 1.3 * (d + 1) * x);
            }
        }
        f[X] = cos(2 * M_PI * n * X.coordinate(db) / L) * g;
    }

    hila::k_binning b;
    auto sd = b.spectraldensity(f);

    lattice.block(blocking);
    fb.block_from(f);
    auto sdb = b.spectraldensity(fb);
    auto csb = b.cross_spectrum(fb, fb);
    lattice.unblock();

    double diff = 0, cdiff = 0, norm = 0;
    for (int i = 0; i < sd.size(); i++) {
        diff += fabs(sd[i] - sdb[i]);
        cdiff += fabs(sd[i] - csb[i]);
        norm += fabs(sd[i]);
    }
    report_pass("Spectral density, anisotropic blocking", hila::broadcast(diff / norm), 1e-10);
    report_pass("Cross spectrum, anisotropic blocking", hila::broadcast(cdiff / norm), 1e-10);
}

//--------------------------------------------------------------------------------

void test_blocking() {

    CoordinateVector blocking;
    blocking.fill(2);
    if (lattice.can_block(blocking)) {

        Field<CoordinateVector> cvf, cvfb;
        cvf[ALL] = X.coordinates();

        GaugeField<SU<2, float>> gf;
        gf = -1;

        hila::print_dashed_line();
        lattice.block(blocking);
        hila::out0 << "Testing blocked lattice of size " << lattice.size() << '\n';

        cvfb.block_from(cvf);
        double sum = 0;
        onsites (ALL) {
            sum += (cvfb[X] - 2 * X.coordinates()).squarenorm();
            cvfb[X] *= -1;
        }

        report_pass("Field blocking test", sum, 1e-5);

        gf.block_gauge_to_current_lattice();

        Complex<float> one(1, 0);
        sum = 0;
        foralldir (d) {
            onsites (ALL) {
                sum += (gf[d][X] - one).squarenorm();
            }
        }

        report_pass("GaugeField blocking test", sum, 1e-5);

        test_site_access();
        test_set_elements_and_select();
        test_fft();
        test_fft_real();

        lattice.unblock();
        cvfb.unblock_to(cvf);

        hila::print_dashed_line();
        hila::out0 << "Return lattice to size " << lattice.size() << '\n';


        sum = 0;
        CoordinateVector divisor;
        divisor.fill(2);
        onsites (ALL) {
            if (X.coordinates().is_divisible(divisor)) {
                sum += (X.coordinates() + cvf[X]).squarenorm();
            }
        }

        report_pass("Field unblocking test", sum, 1e-5);

        test_site_access();
        test_set_elements_and_select();
    }
}


//////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {

    hila::initialize(argc, argv);

    hila::input par("parameters");

    CoordinateVector lsize = par.get("lattice size"); // reads NDIM numbers
    long seed = par.get("random seed");

    par.close();

    // setting up the lattice is convenient to do after reading
    // the parameters
    lattice.setup(lsize);

    // We need random number here
    hila::seed_random(seed);


    ///////////////////////////////////////////////////////////////
    // start tests

    test_reductions();
    test_globals();
    test_namespace_field_type();
    test_functions();
    test_site_access();
    test_minmax();
    test_random();
    test_set_elements_and_select();
    test_shift();
    test_subvolumes();
    test_matrix_operations();
    test_element_operations();
    test_fft();
    test_fft_real();
    test_spectraldensity();
    test_cross_spectrum();
    test_matrix_algebra();
    test_extended();
    test_clusters();
    test_blocking();
    test_spectraldensity_blocking();

    hila::finishrun();
}
