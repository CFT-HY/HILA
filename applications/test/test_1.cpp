
#include "plumbing/field.h"

void transformer_control(const char *);

template <typename T> void sub(Field<T> &a, const Field<T> &b, parity p) { a[p] = b[X]; }

// #include "sub.h"

#ifdef kissa
Miau
#endif

    int
    main() {
    int i;
    Field<double> lf;
    Field<double> dd;
    double dp[10], somevar, *b, c;

    // cout << "Starting..\n";

    // a[EVEN];

    parity p = ALL;
    lf[ALL] = dd[X] + dp[1];

    lf = dd + dp[1];

    double t;
    onsites(p) {
        for (int k = 0; k < 10; k++) {
            lf[X] = dd[X + Direction::xup] + *b + c;
        }
        // transformer_control("dump-ast");
        t += dd[X];
        dd[X] *= lf[X];
    }

    somevar = t + 1.0;
    t = 10;

    sub(lf, dd, EVEN);

    std::cout << "After 1st\n";

    return 0;
}
