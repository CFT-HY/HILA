// Test of field functions: functions which access fields at X and at the
// neighbours X + dir, called from site loops.  The results are compared with
// the same computation written directly in the loops.

#include "hila.h"

// plain function, constant and argument directions
double gradient(const Field<double> &f, const X_index_type &X, Direction d) {
    return 0.5 * (f[X + d] - f[X - d]);
}

// template, directions from a loop inside the function, read at X
template <typename T>
T laplacian(const Field<T> &f, const X_index_type &X) {
    T r = -2 * NDIM * f[X];
    foralldir(d) r += f[X + d] + f[X - d];
    return r;
}

// two fields and X methods
double weighted(const Field<double> &a, const Field<Complex<double>> &b, const X_index_type &X) {
    return a[X + e_x] * b[X - e_y].re + X.coordinate(e_z) + (X.parity() == EVEN ? 1 : 0);
}

// helper functions, a constructor and a global variable.  On GPUs these must
// become __device__
hila::global<double> scale;

struct Pair {
    double a, b;
    Pair(double x, double y) : a(x), b(y) {}
};

double square_helper(double x) {
    return x * x;
}

double helper2(double x) {
    return square_helper(x) + 1;
}

double with_helpers(const Field<double> &f, const X_index_type &X, Direction d) {
    Pair p(f[X + d], f[X]);
    return helper2(p.a) * scale() + p.b;
}

// declared before the loop, defined after it
double forward_diff(const Field<double> &f, const X_index_type &X, Direction d);

// template with a direction parameter: field values passed as parameters
template <typename T>
T diff(const Field<T> &f, const X_index_type &X, Direction d) {
    return f[X + d] - f[X];
}

// a field read only at X, and a scalar parameter
double add_value(const Field<double> &f, const X_index_type &X, double v) {
    return f[X] + v;
}

// offsets: axis offset X + 2*d (deep halo if HALO_DEPTH >= 2), diagonal with two
// parameters, and a constant diagonal
double second_diff(const Field<double> &f, const X_index_type &X, Direction d) {
    return f[X + 2 * d] - 2 * f[X] + f[X - 2 * d];
}

double diagonal(const Field<double> &f, const X_index_type &X, Direction d1, Direction d2) {
    return f[X + d1 + d2] - f[X + d1 - d2] + f[X + e_x - e_z];
}

// unused Field parameter and a default direction
double default_dir(const Field<double> &unused, const Field<double> &f, const X_index_type &X,
                   Direction d = e_y) {
    return f[X - d];
}

// sum |a - b| / sum |b|, b should be non-zero
double rel_diff(const Field<double> &a, const Field<double> &b) {
    double d = 0, n = 0;
    onsites(ALL) {
        d += abs(a[X] - b[X]);
        n += abs(b[X]);
    }
    if (n == 0)
        return -1; // flag trivial test
    return d / n;
}

int main(int argc, char **argv) {
    hila::initialize(argc, argv);
    lattice.setup({8, 8, 8});
    hila::seed_random(1);

    Field<double> f, g, ref;
    Field<Complex<double>> c, lc, lcref;
    f.gaussian_random();
    c.gaussian_random();

    // 1: constant direction
    onsites(ALL) g[X] = gradient(f, X, e_y);
    onsites(ALL) ref[X] = 0.5 * (f[X + e_y] - f[X - e_y]);
    hila::out0 << "gradient e_y:          " << rel_diff(g, ref) << '\n';

    // 2: loop local direction in the loop
    onsites(ALL) {
        double s = 0;
        foralldir(d) s += gradient(f, X, d) * (d + 1);
        g[X] = s;
    }
    onsites(ALL) {
        double s = 0;
        foralldir(d) s += 0.5 * (f[X + d] - f[X - d]) * (d + 1);
        ref[X] = s;
    }
    hila::out0 << "gradient loop dir:     " << rel_diff(g, ref) << '\n';

    // 3: template, real and complex
    onsites(ALL) g[X] = laplacian(f, X);
    onsites(ALL) {
        double r = -2 * NDIM * f[X];
        foralldir(d) r += f[X + d] + f[X - d];
        ref[X] = r;
    }
    hila::out0 << "laplacian double:      " << rel_diff(g, ref) << '\n';

    onsites(ALL) lc[X] = laplacian(c, X);
    onsites(ALL) {
        Complex<double> r = -2 * NDIM * c[X];
        foralldir(d) r += c[X + d] + c[X - d];
        lcref[X] = r;
    }
    onsites(ALL) {
        g[X] = lc[X].re + lc[X].im;
        ref[X] = lcref[X].re + lcref[X].im;
    }
    hila::out0 << "laplacian complex:     " << rel_diff(g, ref) << '\n';

    // 4: EVEN parity, two fields, X methods
    g = 0;
    ref = 0;
    onsites(EVEN) g[X] = weighted(f, c, X);
    onsites(EVEN) ref[X] = f[X + e_x] * c[X - e_y].re + X.coordinate(e_z) + 1;
    hila::out0 << "weighted EVEN:         " << rel_diff(g, ref) << '\n';

    // 5: function defined after the loop
    onsites(ALL) g[X] = forward_diff(f, X, e_z);
    onsites(ALL) ref[X] = f[X + e_z] - f[X];
    hila::out0 << "forward_diff:          " << rel_diff(g, ref) << '\n';

    // 6: helper functions, constructor, global
    scale = 1.5;
    onsites(ALL) g[X] = with_helpers(f, X, -e_x);
    onsites(ALL) ref[X] = (f[X - e_x] * f[X - e_x] + 1) * 1.5 + f[X];
    hila::out0 << "with_helpers:          " << rel_diff(g, ref) << '\n';

    // 7: same neighbours in the loop and in the function, constant and variable direction
    Direction dd = e_z;
    onsites(ALL) g[X] = gradient(f, X, e_y) + f[X + e_y] + gradient(f, X, dd) * f[X - dd];
    onsites(ALL) ref[X] = 0.5 * (f[X + e_y] - f[X - e_y]) + f[X + e_y] +
                          0.5 * (f[X + dd] - f[X - dd]) * f[X - dd];
    hila::out0 << "shared neighbours:     " << rel_diff(g, ref) << '\n';

    // 8: template, real and complex
    onsites(ALL) g[X] = diff(f, X, -e_x);
    onsites(ALL) ref[X] = f[X - e_x] - f[X];
    hila::out0 << "diff template double:  " << rel_diff(g, ref) << '\n';

    onsites(ALL) lc[X] = diff(c, X, dd);
    onsites(ALL) lcref[X] = c[X + dd] - c[X];
    onsites(ALL) {
        g[X] = lc[X].re + lc[X].im;
        ref[X] = lcref[X].re + lcref[X].im;
    }
    hila::out0 << "diff template complex: " << rel_diff(g, ref) << '\n';

    // 9: field written before the call in the same loop is passed at X
    Field<double> h, href;
    onsites(ALL) {
        h[X] = 2 * f[X];
        g[X] = add_value(h, X, 1.0);
    }
    onsites(ALL) ref[X] = 2 * f[X] + 1;
    hila::out0 << "written field at X:    " << rel_diff(g, ref) << '\n';

    // 10: field function call in an argument of a field function
    onsites(ALL) g[X] = add_value(f, X, gradient(f, X, e_x));
    onsites(ALL) ref[X] = f[X] + 0.5 * (f[X + e_x] - f[X - e_x]);
    hila::out0 << "nested call argument:  " << rel_diff(g, ref) << '\n';

    // 11: unused Field parameter, default direction
    onsites(ALL) g[X] = default_dir(c.real(), f, X) + default_dir(f, f, X, e_x);
    onsites(ALL) ref[X] = f[X - e_y] + f[X - e_x];
    hila::out0 << "default direction:     " << rel_diff(g, ref) << '\n';

    // 12: ODD loop writes h and reads its EVEN neighbours through a function
    h = f;
    href = f;
    onsites(ODD) h[X] = gradient(h, X, e_x);
    onsites(ODD) href[X] = 0.5 * (href[X + e_x] - href[X - e_x]);
    hila::out0 << "ODD write, EVEN read:  " << rel_diff(h, href) << '\n';

    // 13: offsets, also used directly in the loop
    onsites(ALL) g[X] = second_diff(f, X, dd) + f[X + 2 * dd];
    onsites(ALL) ref[X] = f[X + 2 * dd] - 2 * f[X] + f[X - 2 * dd] + f[X + 2 * dd];
    hila::out0 << "offset X + 2*d:        " << rel_diff(g, ref) << '\n';

    onsites(ALL) g[X] = diagonal(f, X, e_y, dd) + diagonal(f, X, -e_x, e_y);
    onsites(ALL) ref[X] = (f[X + e_y + dd] - f[X + e_y - dd] + f[X + e_x - e_z]) +
                          (f[X - e_x + e_y] - f[X - e_x - e_y] + f[X + e_x - e_z]);
    hila::out0 << "offset diagonal:       " << rel_diff(g, ref) << '\n';

    hila::finishrun();
}

double forward_diff(const Field<double> &f, const X_index_type &X, Direction d) {
    return f[X + d] - f[X];
}
