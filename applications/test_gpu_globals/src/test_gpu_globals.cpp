/**
 * Test of hila::global variables on GPUs.
 *
 * On GPUs hila::global variables have a device copy, which hilapp places in __device__ memory
 * (default) or in __constant__ memory (GPU_GLOBALS_IN_CONSTANT_MEMORY defined, make CONSTANT=1).
 * The CUDA device linker (seen with CUDA 12.0 and 13.3) does not always lay out __constant__
 * variables of several files correctly: if the constant sections of two object files have equal
 * size, the linked section can be too small, and the program fails at lattice setup with
 *   CUDA error in command: cudaMemcpyToSymbol(_dev_coordinates, ...)
 *   CUDA error string: invalid argument
 *
 * The globals here are chosen to provoke that:
 *   this file:  5 doubles = 40 bytes, same as the hila system constants in hila_gpu.cpp
 *   file_a.cpp: 4 doubles = 32 bytes
 *   file_b.cpp: 4 doubles = 32 bytes, same as file_a.cpp
 *
 * Expected: default build passes.  With CONSTANT=1 the result tells whether the linker of
 * this CUDA version has the problem.  The binary can also be checked without running it:
 *   python3 check_gpu_sections.py build/test_gpu_globals
 */

#include "hila.h"

hila::global<double> m0, m1, m2, m3, m4;

double file_a_value();
double file_b_value();

int main(int argc, char **argv) {

    hila::initialize(argc, argv);

#if defined(CUDA) || defined(HIP)
#ifdef GPU_GLOBALS_IN_CONSTANT_MEMORY
    hila::out0 << "User hila::global variables in __constant__ memory\n";
#else
    hila::out0 << "User hila::global variables in __device__ memory\n";
#endif
#endif

    // fails here if the constant memory layout is broken
    lattice.setup({16, 16, 16});
    hila::out0 << "Lattice setup OK\n";

    m0 = 1;
    m1 = 2;
    m2 = 3;
    m3 = 4;
    m4 = 5;

    Field<double> f;
    onsites (ALL)
        f[X] = m0() + m1() + m2() + m3() + m4();

    double v = f.sum() / lattice.volume();
    double va = file_a_value();
    double vb = file_b_value();

    bool ok = (v == 15) && (va == 100) && (vb == 1000);

    hila::out0 << "Values: main " << v << " (expect 15), file_a " << va << " (expect 100), file_b "
               << vb << " (expect 1000)\n";
    hila::out0 << "TEST hila::global on GPU: " << (ok ? "PASSED" : "FAILED") << '\n';

    hila::finishrun();
}
