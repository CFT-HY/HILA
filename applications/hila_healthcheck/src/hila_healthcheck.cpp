/**
 * @file hila_healthcheck.cpp
 * @author Kari Rummukainen
 * @brief HILA healthcheck application.
 * @details Performs a comprehensive health check on HILA libraries.
 */
#include "hila.h"

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

/////////////////////////////////////////////////////////////////////////////////////

void test_functions() {

    Field<double> df = 0;
    report_pass("Field functions: real exp", exp(df).sum() - lattice.volume(), 1e-8);

    Field<Complex<double>> cf = 0;
    report_pass("Field functions: complex exp", abs(exp(cf).sum() - lattice.volume()), 1e-8);

    df[ALL] = sin(X.x() * 2 * M_PI / lattice.size(e_x));
    report_pass("Field functions: sin", df.sum(), 1e-8);
}


/////////////////////////////////////////////////////////////////////////////////////

void check_reductions() {


    // test reductions

    {
        Field<Complex<double>> f;
        f[ALL] = expi(2 * M_PI * X.coordinate(e_x) / lattice.size(e_x));

        Complex<double> sum = 0;
        onsites(ALL) sum += f[X];

        sum /= lattice.volume();
        report_pass("Complex reduction value " + hila::prettyprint(sum), abs(sum), 1e-4);

        ReductionVector<Complex<double>> rv(lattice.size(e_x));

        onsites(ALL) {
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
        onsites(ALL) {
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
    }

    {
        // reductionvector with long
        Field<long> lf;
        lf[ALL] = X.x();

        ReductionVector<long> rv(lattice.size(e_x));

        onsites(ALL) {
            rv[X.x()] += (lf[X] == X.x());
        }

        long s = 0;
        for (int x = 0; x < rv.size(); x++) {
            s += abs(rv[x] - (lattice.volume() / lattice.size(e_x)));
        }

        report_pass("ReductionVector<long>, sum " + hila::prettyprint(s), s, 1e-15);
    }

    {
        Field<SU<3, double>> mf;
        mf = 1;
        // onsites(ALL) mf[X] = (mf[X] + mf[X])*0.5;

        ReductionVector<SU<3, double>> rmf(lattice.size(e_x));

        onsites(ALL) {
            //    rmf[X.x()] += (mf[X] + mf[X])*0.5;
            rmf[X.x()] += mf[X];
        }

        double diff = 0;
        for (int i = 0; i < lattice.size(e_x); i++) {
            diff += (rmf[i] - lattice.size(e_y) * lattice.size(e_z)).squarenorm();
        }

        report_pass("SU(3) ReductionVector", diff, 1e-8);
    }
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// write something to a location

void test_site_access() {

    Field<Complex<double>> f = 0;
    Complex<double> sum;

    CoordinateVector c;
    for (int i = 0; i < 3; i++) {
        foralldir(d) c[d] = hila::random() * lattice.size(d);
        hila::broadcast(c); // need to broadcast!

        f[c] = 4;   // write to location c
        sum = f[c]; // read back - does mpi comm
        report_pass("Setting and reading a value at " + hila::prettyprint(c.transpose()),
                    abs(sum - 4), 1e-10);
    }
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// test maxloc and minloc operation

void test_minmax() {

    CoordinateVector c, loc;
    Field<double> n;

    for (Parity par : {ALL, EVEN, ODD}) {

        n[par] = hila::random();
        do {
            foralldir(d) c[d] = hila::random() * lattice.size(d);
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
            foralldir(d) c[d] = hila::random() * lattice.size(d);
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

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// test random number properties
// rough test, testing spectrum of gaussians

void test_random() {

    constexpr int n_loops = 100;

    Field<double> f;

    double fsum = 0, fsqr = 0;
    for (int i = 0; i < n_loops; i++) {
        f.gaussian_random();
        double s = 0, s2 = 0;
        onsites(ALL) {
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


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// test access to a list of sites

void test_set_elements_and_select() {

    Field<Complex<double>> f = 0;
    CoordinateVector c;

    std::vector<CoordinateVector> cvec;
    std::vector<Complex<double>> vals;


    if (hila::myrank() == 0) {
        int k = 1;
        for (int i = 0; i <= 50; i++) {
            foralldir(d) c[d] = hila::random() * lattice.size(d);
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

    onsites(ALL) {
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


    onsites(ALL) {}
}

/////////////////////////////////////////////////////////////////////////////////////

void test_subvolumes() {

    bool ok = true;
    for (int i = 0; i < 20; i++) {
        CoordinateVector c;
        foralldir(d) c[d] = hila::random() * lattice.size(d);
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
    foralldir(d) if (d < NDIM - 1) {
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
            foralldir(d2) {
                if (c[d2] >= 0)
                    mycoord[d2] = c[d2];
            }
            foralldir(d2) {
                if (c[d2] < 0) {
                    mycoord[d2] = -1;
                    break;
                }
            }

            for (auto s : slice) {

                // get the coordinate which should be here
                bool add = true;
                foralldir(d2) {
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

/////////////////////////////////////////////////////////////////////////////////////


void fft_test() {


    using T = Complex<double>;

    Field<T> f, p, p2;

    // Start with unit field
    f = 1;

    // After one FFT the field is 0 except at coord 0
    p2 = 0;

    p2[{0, 0, 0}] = lattice.volume();

    FFT_field(f, p);

    double eps = squarenorm_relative(p, p2);

    report_pass("FFT constant field", eps, 1e-13 * sqrt(lattice.volume()));

    //-----------------------------------------------------------------
    // After two applications the field should be back to a constant * volume

    FFT_field(p, f, fft_direction::back);

    double sum = 0;
    double tnorm = 0;
    onsites(ALL) {
        sum += (f[X] - lattice.volume()).squarenorm();
        tnorm += f[X].squarenorm();
    }

    eps = fabs(sum / tnorm);

    report_pass("FFT inverse transform", eps, 1e-10);

    //-----------------------------------------------------------------
    for (int iter = 0; iter < 5; iter++) {
        Vector<NDIM, double> kv;
        CoordinateVector kx;
        foralldir(d) {
            kx[d] = hila::broadcast(hila::random()) * lattice.size(d);
        }

        kv = kx.convert_to_k();


        onsites(ALL) {
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
        onsites(ALL) r[X] = hila::gaussrand();

        f = r.FFT_real_to_complex();
        p = f.FFT(fft_direction::back) / lattice.volume();
        eps = squarenorm_relative(r, p);

        report_pass("FFT real to complex", eps, 1e-13 * sqrt(lattice.volume()));

        auto r2 = f.FFT_complex_to_real(fft_direction::back) / lattice.volume();
        eps = squarenorm_relative(r, r2);

        report_pass("FFT complex to real", eps, 1e-13 * sqrt(lattice.volume()));
    }

    //-----------------------------------------------------------------
    // Check fft norm


    onsites(ALL) {
        p[X] = hila::random() * exp(-X.coordinates().convert_to_k().squarenorm());
    }
    f = p.FFT(fft_direction::back) / sqrt(lattice.volume());

    double nf = f.squarenorm();
    double np = p.squarenorm();
    report_pass("Norm of field = " + hila::prettyprint(nf) + " and FFT = " + hila::prettyprint(np),
                (nf - np) / nf, 1e-10);

    hila::k_binning b;
    b.k_max(M_PI * sqrt(3.0));

    auto bf = b.bin_k_field(p.conj() * p);

    double s = 0;
    for (auto b : bf) {
        s += abs(b);
    }
    hila::broadcast(s);

    report_pass("Norm of binned FFT = " + hila::prettyprint(s), (s - np) / np, 1e-10);
}

//---------------------------------------------------------------------------

void spectraldensity_test() {


    // test spectral density for single waves

    Field<Complex<double>> f, p;

    hila::k_binning b;
    b.k_max(M_PI * sqrt(3.0));

    // test std binning first

    for (int iter = 0; iter < 3; iter++) {
        Vector<NDIM, double> kv;
        CoordinateVector kx;
        foralldir(d) {
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
        onsites(ALL) {
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
    onsites(ALL) {
        s[X] = SiteIndex(X.coordinates());
    }
}

//--------------------------------------------------------------------------------

void test_matrix_operations() {

    Field<Matrix<3, 2, Complex<double>>> mf;

    onsites(ALL) mf[X].fill(1 + I);

    Matrix<3, 3, Complex<double>> cm;
    cm.asArray() = 4;
    double sum = 0;
    onsites(ALL) {
        sum += (mf[X] * mf[X].dagger() - cm).squarenorm();
    }

    report_pass("matrix multiply and addition", sum, 1e-8);

    auto dm = cm * I - 2 * I;
    dm.asArray() *= I;
    dm = ((dm - 2).asArray() + 4).asMatrix();
    report_pass("Array and imaginary unit operations", dm.squarenorm(), 1e-8);
}


//--------------------------------------------------------------------------------

void test_matrix_algebra() {

    using myMatrix = SquareMatrix<4, Complex<double>>;

    Field<myMatrix> M;
    Field<double> delta;

    M.gaussian_random(2.0);

    // eigenvalue test - show that  M = U D U^*, where D is diagonal eigenvalue matrix and U
    // matrix of eigenvectors

    onsites(ALL) {
        auto H = M[X] * M[X].dagger(); // make hermitean

        auto r = H.eigen_hermitean();
        delta[X] = (H - r.eigenvectors * r.eigenvalues * r.eigenvectors.dagger()).norm();
    }

    auto max_delta = delta.max();

    report_pass("Eigenvalue analysis with " + hila::prettyprint(myMatrix::rows()) + "x" +
                    hila::prettyprint(myMatrix::columns()) + " Hermitean matrix",
                max_delta, 1e-10);

    // Singular value test - non-pivoted

    onsites(ALL) {
        auto r = M[X].svd();
        delta[X] = (M[X] - r.U * r.singularvalues * r.V.dagger()).norm();
    }

    max_delta = delta.max();

    report_pass("SVD with " + hila::prettyprint(myMatrix::rows()) + "x" +
                    hila::prettyprint(myMatrix::columns()) + " Complex matrix",
                max_delta, 1e-10);

    // pivoted singular values

    M.gaussian_random();

    onsites(ALL) {
        auto r = M[X].svd_pivot(hila::sort::ascending);
        delta[X] = (M[X] - r.U * r.singularvalues * r.V.dagger()).norm();
    }

    max_delta = delta.max();

    report_pass("Fully pivoted SVD with " + hila::prettyprint(myMatrix::rows()) + "x" +
                    hila::prettyprint(myMatrix::columns()) + " Complex matrix",
                max_delta, 1e-10);
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

    check_reductions();

    test_functions();

    test_site_access();

    test_minmax();

    test_random();

    test_set_elements_and_select();

    test_subvolumes();

    test_matrix_operations();

    fft_test();

    spectraldensity_test();

    test_matrix_algebra();

    hila::finishrun();
}
