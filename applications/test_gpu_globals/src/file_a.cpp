// See test_gpu_globals.cpp: 4 doubles = 32 bytes of globals in this file

#include "hila.h"

hila::global<double> a0, a1, a2, a3;

double file_a_value() {
    a0 = 1 * 10;
    a1 = 2 * 10;
    a2 = 3 * 10;
    a3 = 4 * 10;

    Field<double> f;
    onsites (ALL)
        f[X] = a0() + a1() + a2() + a3();

    return f.sum() / lattice.volume();
}
