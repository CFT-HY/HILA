
#include "plumbing/defs.h"
#include "datatypes/cmplx.h"

#include "plumbing/field.h"

// extern Field<int> glob;

Complex<double> d(Complex<double> x) { return x; }
Complex<double> e(Complex<double> x) { return d(x); }
// #pragma hila ast dump
Complex<double> f(const Complex<double> &x) { return e(x); }

template <typename T> T xyz(out_only T &v) {
    v = 1.3;
    return sin(v);
}

using ft = Complex<double>;

template <typename T> class v2 {
  public:
    using base_type = hila::number_type<T>;
    Complex<T> a[2];

    void setter() out_only { a[0] = 1; }
};

// template <typename g>
// Field<g> f2(const Field<g> &f, int i);

double dv(double *d) { return *d + 1; }

template <typename g> void pf(Field<g> &s) { s[ALL] = 1; }

int main() {

    Field<Complex<double>> a, b, c;
    int i;
    Field<double> t(1.0), s;
    Field<Complex<float>> kissa;

    auto y = X;

    // pf(t);

    Field<v2<double>> A;

    CoordinateVector v = e_x - 2 * e_y;

    Parity p = ODD;

    // a[ALL] = b[X+2*e_x+e_y];

    a = b.shift(v);

    Direction d = e_x, d2 = e_y;

    A[ALL] = {Complex(1, 0), Complex(0, 0)};

    double dvar;
    onsites(p) {

        kissa[X] = 2;
        A[X].setter();
        auto tv = t[X];
        s[X] = xyz(tv);
        // t[X] = tmp;
    }

    return 0;
}
