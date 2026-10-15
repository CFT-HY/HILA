/**
 * Prototype of row-structured site loops in the lexicographic layout (EVEN_SITES_FIRST=0),
 * CPU only.
 *
 * Sites are traversed row by row (x fastest).  Within a row every neighbour index is
 * affine in x: for directions other than x it is nb_base[d] + x, where nb_base[d] is the
 * neighbour of the first site of the row (on-node, wrapped on node or in the halo, all
 * have unit stride in x).  For x the neighbours are i+1 and i-1, except at the row ends,
 * which use the general site_layout neighbour().  Thus the inner loop has no branches and
 * only unit-stride accesses, and the compiler can vectorize it.
 *
 * Compares hand-written row loops against the hilapp-generated onsites() loops and checks
 * that they give identical results.
 */

#include "hila.h"

#ifdef EVEN_SITES_FIRST
#error "bench_rows: compile with EVEN_SITES_FIRST=0"
#endif
#if defined(CUDA) || defined(HIP)
#error "bench_rows: CPU only"
#endif

using cf = Complex<float>;
using cd = Complex<double>;
using su3f = Matrix<3, 3, Complex<float>>;

/////////////////////////////////////////////////////////////////////////////////
/// Row loop driver: body(i, nb) for each site i, with nb[d] the neighbour index to
/// direction d.  The interior of each row calls body with nb[d] = nb_base[d] + x.
/////////////////////////////////////////////////////////////////////////////////

template <typename Body>
inline void row_loop(Body &&body) {
    const auto &lay = lattice->mynode.layout;
    const int nx = lay.size[0];
    const unsigned rows = lay.volume / nx;

    for (unsigned r = 0; r < rows; r++) {
        const unsigned i0 = r * nx;
        const hila::layout_site s0 = lay.site(i0);

        // neighbours of the 1st site of the row
        unsigned nbb[NDIRS];
        for (int d = 0; d < NDIRS; d++)
            nbb[d] = lay.neighbour(s0, d);

        // row ends: general neighbours
        auto general = [&](int x) {
            unsigned nb[NDIRS];
            hila::layout_site s = lay.site(i0 + x);
            for (int d = 0; d < NDIRS; d++)
                nb[d] = lay.neighbour(s, d);
            body(i0 + x, nb);
        };

        general(0);
        // interior of the row: unit stride everywhere
        for (int x = 1; x < nx - 1; x++) {
            unsigned nb[NDIRS];
            for (int d = 0; d < NDIRS; d++)
                nb[d] = nbb[d] + x;
            nb[e_x] = i0 + x + 1;
            nb[-e_x] = i0 + x - 1;
            body(i0 + x, nb);
        }
        if (nx > 1)
            general(nx - 1);
    }
}

/// gather all directions (halo valid for the raw loops below)
template <typename T>
void gather_all(const Field<T> &f) {
    for (int d = 0; d < NDIRS; d++)
        f.start_gather((Direction)d, ALL);
    for (int d = 0; d < NDIRS; d++)
        f.wait_gather((Direction)d, ALL);
}

/////////////////////////////////////////////////////////////////////////////////
// hilapp-generated versions

template <typename T>
void laplace_onsites(Field<T> &a, const Field<T> &b) {
    onsites (ALL) {
        T s = -2.0f * NDIM * b[X];
        foralldir (d)
            s += b[X + d] + b[X - d];
        a[X] = s;
    }
}

template <typename T>
void nnsum_onsites(Field<T> &a, const Field<T> &b) {
    onsites (ALL) {
        T s = 0;
        foralldir (d)
            s += b[X + d] + b[X - d];
        a[X] = s;
    }
}

/////////////////////////////////////////////////////////////////////////////////
// row loop versions on raw storage

template <typename T>
void laplace_rows(Field<T> &a, const Field<T> &b) {
    gather_all(b);
    T *RESTRICT A = a.field_buffer();
    const T *RESTRICT B = b.field_buffer();
    row_loop([&](unsigned i, const unsigned *nb) {
        T s = -2.0f * NDIM * B[i];
        for (int d = 0; d < NDIM; d++)
            s += B[nb[d]] + B[nb[NDIRS - 1 - d]];
        A[i] = s;
    });
    a.mark_changed(ALL);
}

template <typename T>
void nnsum_rows(Field<T> &a, const Field<T> &b) {
    gather_all(b);
    T *RESTRICT A = a.field_buffer();
    const T *RESTRICT B = b.field_buffer();
    row_loop([&](unsigned i, const unsigned *nb) {
        T s = B[nb[0]] + B[nb[NDIRS - 1]];
        for (int d = 1; d < NDIM; d++)
            s += B[nb[d]] + B[nb[NDIRS - 1 - d]];
        A[i] = s;
    });
    a.mark_changed(ALL);
}

/////////////////////////////////////////////////////////////////////////////////

template <typename F>
double time_kernel(int iters, F &&kernel) {
    kernel(); // warm-up, halo exchange
    hila::synchronize();
    double t0 = hila::gettime();
    for (int i = 0; i < iters; i++)
        kernel();
    hila::synchronize();
    return (hila::gettime() - t0) / iters * 1e3;
}

template <typename T>
void compare(const std::string &name, int iters, void (*k_on)(Field<T> &, const Field<T> &),
             void (*k_rows)(Field<T> &, const Field<T> &)) {
    Field<T> a1, a2, b;
    b.gaussian_random();
    a1 = 0;
    a2 = 0;
    double t_on = time_kernel(iters, [&]() { k_on(a1, b); });
    double t_rows = time_kernel(iters, [&]() { k_rows(a2, b); });
    double diff = (a1 - a2).squarenorm();
    hila::out0 << std::left << std::setw(24) << name << std::right << std::fixed
               << std::setprecision(4) << std::setw(10) << t_on << " ms " << std::setw(10)
               << t_rows << " ms  speedup " << std::setprecision(2) << std::setw(5)
               << t_on / t_rows << "  diff " << std::scientific << std::setprecision(1) << diff
               << '\n';
}

template <typename T>
void copy_kernel(Field<T> &a, const Field<T> &b) {
    onsites (ALL)
        a[X] = b[X];
}

int main(int argc, char **argv) {

    hila::initialize(argc, argv);

    hila::input par("parameters");
    CoordinateVector lsize = par.get("lattice size");
    int iters = par.get("iterations");
    par.close();

    lattice.setup(lsize);
    hila::seed_random(1);

    hila::out0 << "\n=== row loops vs onsites(), lexicographic layout, lattice " << lsize << ", "
               << hila::number_of_nodes() << " ranks ===\n";
    hila::out0 << std::left << std::setw(24) << "kernel" << std::right << std::setw(13)
               << "onsites" << std::setw(13) << "rows" << '\n';

    {
        Field<cf> a, b;
        b.gaussian_random();
        double t = time_kernel(iters, [&]() { copy_kernel(a, b); });
        hila::out0 << std::left << std::setw(24) << "copy C<float> (ref)" << std::right
                   << std::fixed << std::setprecision(4) << std::setw(10) << t << " ms\n";
    }

    compare<cf>("laplace C<float>", iters, laplace_onsites<cf>, laplace_rows<cf>);
    compare<cd>("laplace C<double>", iters, laplace_onsites<cd>, laplace_rows<cd>);
    compare<su3f>("nn-sum SU3<float>", iters, nnsum_onsites<su3f>, nnsum_rows<su3f>);

    hila::finishrun();
}
