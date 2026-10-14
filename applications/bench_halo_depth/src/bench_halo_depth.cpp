/**
 * Timing of f[X + k*dir] access.  With HALO_DEPTH >= k hilapp reads the deep halo directly
 * (gather_depth + neighbour_offset), with HALO_DEPTH=1 it falls back to shift-copy.
 * Compile both and compare.
 *
 * "static": source field is not changed between iterations, so a deep halo gather is
 *           done only once (shift-copy still copies every time).
 * "changed": source field is marked changed before each call, so communication is
 *           included each iteration.
 * nn-kernels are timed for reference (should not depend on HALO_DEPTH).
 */

#include "hila.h"

using cf = Complex<float>;
using su3f = Matrix<3, 3, Complex<float>>;

#if NDIM != 3
#error "bench_halo_depth: NDIM=3 only"
#endif

template <typename T>
void nn_sum(Field<T> &a, const Field<T> &b) {
    onsites (ALL) {
        T s = b[X + e_x] + b[X - e_x];
        s += b[X + e_y] + b[X - e_y];
        s += b[X + e_z] + b[X - e_z];
        a[X] = s;
    }
}

// 2-step "laplacian", all axis directions
template <typename T>
void nnn_sum(Field<T> &a, const Field<T> &b) {
    onsites (ALL) {
        T s = b[X + 2 * e_x] + b[X - 2 * e_x];
        s += b[X + 2 * e_y] + b[X - 2 * e_y];
        s += b[X + 2 * e_z] + b[X - 2 * e_z];
        a[X] = s;
    }
}

// typical improved derivative: 1- and 2-step terms in one direction
template <typename T>
void deriv4(Field<T> &a, const Field<T> &b, Direction d) {
    onsites (ALL) {
        a[X] = (b[X + d] - b[X - d]) * (2.0f / 3) - (b[X + 2 * d] - b[X - 2 * d]) * (1.0f / 12);
    }
}

// 3-step along one direction
template <typename T>
void step3(Field<T> &a, const Field<T> &b) {
    onsites (ALL) {
        a[X] = b[X + 3 * e_x] + b[X - 3 * e_x];
    }
}

template <typename F>
void time_kernel(const std::string &name, int iters, F &&kernel) {
    kernel(); // warm-up
    hila::synchronize();
    double t0 = hila::gettime();
    for (int i = 0; i < iters; i++)
        kernel();
    hila::synchronize();
    double t = (hila::gettime() - t0) / iters;
    hila::out0 << std::left << std::setw(30) << name << std::right << std::fixed
               << std::setprecision(4) << std::setw(10) << t * 1e3 << " ms\n";
}

template <typename T>
void bench(const std::string &tname, int iters) {
    Field<T> a, b;
    b.gaussian_random();
    a = 0;

    // static source
    time_kernel(tname + " nn-sum", iters, [&]() { nn_sum(a, b); });
    time_kernel(tname + " 2-step sum", iters, [&]() { nnn_sum(a, b); });
    time_kernel(tname + " deriv4 e_y", iters, [&]() { deriv4(a, b, e_y); });
    if (HALO_DEPTH == 1 || HALO_DEPTH >= 3)
        time_kernel(tname + " 3-step e_x", iters, [&]() { step3(a, b); });

    // source changed every iteration: communication included
    time_kernel(tname + " changed: mark only", iters, [&]() { b.mark_changed(ALL); });
    time_kernel(tname + " changed: nn-sum", iters, [&]() {
        b.mark_changed(ALL);
        nn_sum(a, b);
    });
    time_kernel(tname + " changed: 2-step sum", iters, [&]() {
        b.mark_changed(ALL);
        nnn_sum(a, b);
    });
    time_kernel(tname + " changed: deriv4 e_y", iters, [&]() {
        b.mark_changed(ALL);
        deriv4(a, b, e_y);
    });

    hila::out0 << tname << " checksum " << a.sum() << '\n';
}

int main(int argc, char **argv) {

    hila::initialize(argc, argv);

    hila::input par("parameters");
    CoordinateVector lsize = par.get("lattice size");
    int iters = par.get("iterations");
    par.close();

    lattice.setup(lsize);
    hila::seed_random(1);

    hila::out0 << "\n=== HALO_DEPTH " << HALO_DEPTH
               << (HALO_DEPTH > 1 ? " (deep halo access)" : " (shift-copy)") << " ===\n";
    hila::out0 << "Lattice " << lsize << ", " << hila::number_of_nodes() << " ranks, " << iters
               << " iterations\n";

    bench<cf>("C<float>", iters);
    bench<su3f>("SU3<float>", iters);

    hila::finishrun();
}
