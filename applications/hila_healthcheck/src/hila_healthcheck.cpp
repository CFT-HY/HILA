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


    if (hila::myrank() == 0) {
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
            if (hila::myrank() == 0) {

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

        report_pass("FFT Complex<float> constant field", eps, 1e-6 * sqrt(lattice.volume()));

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
        sum += (hila::elem::sub(mf[X],1)).squarenorm();
        sum2 += hila::elem::log(mf[X]).squarenorm();
        sum3 += (hila::elem::mul(hila::elem::add(mf[X],-2*mf[X]),-1) - mf[X]).squarenorm();
    }

    report_pass("hila::elem::exp, log, sub, add, mul", sum + sum2 + sum3, 1e-8);


}


/**
 * @brief Test extended type
 * @details Test extended type for sums that exibit loss in accuracy with double.
 *
 */
void test_extended() {

#if NDIM > 1
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

    for (double mag = 1e8; mag <= 1e+32; mag *= 1e8) {

        size_t nsqr = 0, nmag = 0, n1 = 0;

        onsites (ALL) {
            if (X.x() % 2 == 0) {
                if (X.y() % 2 == 0) {
                    g[X] = sqr(mag);
                    nsqr += 1;
                } else {
                    g[X] = mag;
                    nmag += 1;
                }
            } else {
                g[X] = 1;
                n1 += 1;
            }
        }
        ExtendedPrecision ev = 0;
        double s = 0;
        onsites (ALL) {
            ev += g[X];
            s += g[X];
        }

        double result = n1 + nmag * mag + nsqr * sqr(mag);

        // ev -= lattice.volume() / 2;
        // ev -= mag * (lattice.volume() / 4);
        // ev -= sqr(mag) * (lattice.volume() / 4);

        // above multiplication loses precision!  Subtracting
        // here step-by-step

        ExtendedPrecision res = 0, r;
        r = sqr(mag);
        r *= nsqr;
        res = r;

        r = mag;
        r *= nmag;
        res += r;

        r = 1;
        r *= n1;
        res += r;

        ev -= res;
        s -= res.to_double();


        // for (int i = 0; i < n1; i++) {
        //     ev -= 1;
        //     s -= 1;
        // }


        // for (int i = 0; i < lattice.volume() / 4; i++) {
        //     ev -= mag;
        //     ev -= sqr(mag);
        //     s -= mag;
        //     s -= sqr(mag);
        // }

        // s -= sqr(mag) * lattice.volume() / 4;
        // s -= mag * lattice.volume() / 4;
        // s -= lattice.volume() / 2;

        // hila::out0 << "RES " << ev.value << " + " << ev.compensation << " + " <<
        // ev.compensation2
        //           << '\n';


        std::stringstream ssres;
        ssres << "Extended reduction residual w. delta " << mag << " : " << ev.to_double() / result
              << " (double " << s / result << ")";

        report_pass(ssres.str(), ev.to_double() / result, 1e-20);
    }
#endif
}


//--------------------------------------------------------------------------------

void test_clusters() {

#if NDIM > 2

    Field<int> m;

    m = hila::clusters::background;

    onsites (ALL) {
        auto c = X.coordinates();
        if (c[e_x] == 0 && c[e_y] == 0)
            m[X] = 1;
        else if (c[e_x] == 2 && c[e_z] == 2)
            m[X] = 2;
        else if (c[e_x] == 4 && c[e_y] < 3 && c[e_z] > 1 && c[e_z] <= 4)
            m[X] = 3;
    }

    hila::clusters cl(m);

    report_pass("Cluster test: number of clusters ", cl.number() - 3, 1e-10);
    if (cl.number() == 3) {
        double sumsize =
            fabs(cl.size(0) - (lattice.volume() / (lattice.size(e_x) * lattice.size(e_y)))) +
            fabs(cl.size(1) - (lattice.volume() / (lattice.size(e_x) * lattice.size(e_z)))) +
            fabs(cl.size(2) - 1 * 3 * 3);
        report_pass("Cluster test: cluster sizes ", sumsize, 1e-10);

        double types = abs(cl.type(0) - 1) + abs(cl.type(1) - 2) + abs(cl.type(2) - 3);
        report_pass("Cluster test: cluster types ", types, 1e-10);

#if NDIM == 3
        double area = abs(cl.area(0) - 4 * lattice.size(e_z)) +
                      abs(cl.area(1) - 4 * lattice.size(e_y)) +
                      abs(cl.area(2) - 2 * (1 * 3 + 1 * 3 + 3 * 3));

        report_pass("Cluster test: cluster area ", area, 1e-10);

#endif
    }
#endif
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
    test_functions();
    test_site_access();
    test_minmax();
    test_random();
    test_set_elements_and_select();
    test_subvolumes();
    test_matrix_operations();
    test_element_operations();
    test_fft();
    test_spectraldensity();
    test_matrix_algebra();
    test_extended();
    test_clusters();
    test_blocking();

    hila::finishrun();
}
