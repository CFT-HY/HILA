// See test_gpu_globals.cpp: 4 doubles = 32 bytes of globals in this file

#include "hila.h"

hila::global<double> b0, b1, b2, b3;

double file_b_value() {
    b0 = 1 * 100;
    b1 = 2 * 100;
    b2 = 3 * 100;
    b3 = 4 * 100;

    Field<double> f;
    onsites (ALL)
        f[X] = b0() + b1() + b2() + b3();

    return f.sum() / lattice.volume();
}
