/**
 * Timing of typical onsites() kernels, for comparing LAYOUT_ON_THE_FLY against the
 * neighbour/coordinate arrays.  Fields are not changed between iterations of a kernel,
 * so after the first call no communication is done - this measures the loop kernels.
 */

#include "hila.h"

using cf = Complex<float>;
using cd = Complex<double>;
using su3f = Matrix<3, 3, Complex<float>>;

void copy_cf(Field<cf> &a, const Field<cf> &b) {
    onsites (ALL)
        a[X] = b[X];
}

void laplace_cf(Field<cf> &a, const Field<cf> &b, Parity par) {
    onsites (par) {
        cf s = -2.0f * NDIM * b[X];
        foralldir (d)
            s += b[X + d] + b[X - d];
        a[X] = s;
    }
}

void laplace_cd(Field<cd> &a, const Field<cd> &b) {
    onsites (ALL) {
        cd s = -2.0 * NDIM * b[X];
        foralldir (d)
            s += b[X + d] + b[X - d];
        a[X] = s;
    }
}

void nn_sum_su3(Field<su3f> &a, const Field<su3f> &b) {
    onsites (ALL) {
        su3f s = 0;
        foralldir (d)
            s += b[X + d] + b[X - d];
        a[X] = s;
    }
}

void coords_cf(Field<cf> &a) {
    onsites (ALL)
        a[X] = cf(X.coordinate(e_x) + X.coordinate(e_y), X.coordinate(e_z));
}

void parity_cf(Field<cf> &a, const Field<cf> &b) {
    onsites (ALL) {
        if (X.parity() == EVEN)
            a[X] = b[X];
        else
            a[X] = -b[X];
    }
}

template <typename F>
void time_kernel(const char *name, int iters, double bytes_per_site, F &&kernel) {
    kernel(); // warm-up, also does the (one-time) halo exchange
    hila::synchronize();
    double t0 = hila::gettime();
    for (int i = 0; i < iters; i++)
        kernel();
    hila::synchronize();
    double t = (hila::gettime() - t0) / iters;
    double gbs = bytes_per_site * lattice.volume() / t / 1e9 / hila::number_of_nodes();
    hila::out0 << std::left << std::setw(22) << name << std::right << std::fixed
               << std::setprecision(4) << std::setw(10) << t * 1e3 << " ms  " << std::setw(8)
               << std::setprecision(1) << gbs << " GB/s/rank (min. traffic)\n";
}

int main(int argc, char **argv) {

    hila::initialize(argc, argv);

    hila::input par("parameters");
    CoordinateVector lsize = par.get("lattice size");
    int iters = par.get("iterations");
    par.close();

    lattice.setup(lsize);
    hila::seed_random(1);

    Field<cf> a, b;
    Field<cd> ad, bd;
    Field<su3f> ma, mb;

    b.gaussian_random();
    bd.gaussian_random();
    mb.gaussian_random();
    a = 0;
    ad = 0;
    ma = 0;

#ifdef LAYOUT_ON_THE_FLY
    hila::out0 << "\n=== LAYOUT_ON_THE_FLY ===\n";
#else
    hila::out0 << "\n=== neighbour/coordinate arrays (default) ===\n";
#endif
    hila::out0 << "Lattice " << lsize << ", " << hila::number_of_nodes() << " ranks, " << iters
               << " iterations\n";

    time_kernel("copy Complex<float>", iters, 2 * sizeof(cf), [&]() { copy_cf(a, b); });
    time_kernel("laplace Complex<float>", iters, 2 * sizeof(cf),
                [&]() { laplace_cf(a, b, ALL); });
    time_kernel("laplace C<float> EVEN", iters, sizeof(cf), [&]() { laplace_cf(a, b, EVEN); });
    time_kernel("laplace Complex<dbl>", iters, 2 * sizeof(cd), [&]() { laplace_cd(ad, bd); });
    time_kernel("nn-sum SU3<float>", iters, 2 * sizeof(su3f), [&]() { nn_sum_su3(ma, mb); });
    time_kernel("coordinates", iters, sizeof(cf), [&]() { coords_cf(a); });
    time_kernel("parity", iters, 2 * sizeof(cf), [&]() { parity_cf(a, b); });

    // checksum so that the compiler cannot remove anything
    hila::out0 << "checksum " << a.sum() << ' ' << ad.sum() << ' ' << ma.sum().trace() << '\n';

    hila::finishrun();
}
