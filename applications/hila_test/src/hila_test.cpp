
#include "hila.h"


// report the result of the test -- TODO: nicer formatting?
bool report_pass(std::string message, double eps, double limit) {
    if (eps < limit) {
        hila::out0 << "--- " << message << " passed" << std::endl;
        return true;
    } else {
        hila::out0 << "*** " << message << " FAILED: eps " << eps << " limit " << limit
                   << std::endl;
        return false;
    }
}


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

    report_pass("FFT constant field", eps, 1e-10);

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

        kv = convert_to_k(kx);

        onsites(ALL) {
            double d = kv.dot(X.coordinates());
            f[X] = expi(d);
        }

        FFT_field(f, p);

        p2 = 0;
        p2[kx] = lattice.volume();

        eps = squarenorm_relative(p, p2);

        report_pass("FFT of wave vector" + hila::prettyprint(kx), eps, 1e-10);
    }

    Field<double> r;
    onsites(ALL) r[X] = hila::gaussrand();

    f = r.FFT_real_to_complex();
    p = f.FFT(fft_direction::back) / lattice.volume();
    eps = squarenorm_relative(r, p);

    report_pass("FFT real to complex", eps, 1e-10);

    auto r2 = f.FFT_complex_to_real(fft_direction::back) / lattice.volume();
    eps = squarenorm_relative(r, r2);

    report_pass("FFT complex to real", eps, 1e-10);
}

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
        kv = convert_to_k(kx);
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

        report_pass("Binning test at vector " + hila::prettyprint(kx), sum, 1e-10);

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
    // FFT tests

    fft_test();

    spectraldensity_test();

    hila::finishrun();
}
