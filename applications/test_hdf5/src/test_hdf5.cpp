/**
 * @file test_hdf5.cpp
 * @brief Test of HDF5 output (hila::hdf5_file).
 *
 * Writes test.h5 with datasets of several types (scalar, complex, SU(3), single precision,
 * integer, CoordinateVector, rank 0 arrays), groups and attributes, together with raw binary
 * reference files (*.bin) of the same fields.  Validate with h5py:
 *
 *     mpirun -np 4 ./build/test_hdf5
 *     python3 check_hdf5.py
 *
 * The lattice size 16 x 8 x 12 is hard-coded in check_hdf5.py.
 */

#include "hila.h"

int main(int argc, char **argv) {
    hila::initialize(argc, argv);
    lattice.setup({16, 8, 12});
    hila::seed_random(1);

    Field<double> s;
    Field<Complex<double>> c;
    Field<SU<3, double>> u;
    Field<int> k;
    Field<CoordinateVector> cv;
    s.gaussian_random();
    c.gaussian_random();
    u.gaussian_random();
    onsites (ALL) {
        k[X] = X.coordinate(e_x) + 100 * X.coordinate(e_y) + 10000 * X.coordinate(e_z);
        cv[X] = X.coordinates();
    }

    // reference raw binary files
    s.write("s.bin");
    c.write("c.bin");
    u.write("u.bin");

    // data held on rank 0
    std::vector<double> v;
    std::vector<Complex<float>> vc;
    if_rank0 () {
        for (int i = 0; i < 5 * 3; i++)
            v.push_back(i * 0.5);
        for (int i = 0; i < 4; i++)
            vc.push_back(Complex<float>(i, -i));
    }

    {
        hila::hdf5_file h5("test.h5");
        h5.attribute("beta", 6.0)
            .attribute("ntraj", 100)
            .attribute("lattice size", lattice.size())
            .attribute("name", "hila test")
            .attribute("vec", std::vector<double>{1.5, 2.5})
            .attribute("cplx", Complex<double>(1, 2))
            .attribute("i64", (int64_t)1 << 40);
        h5.write("s", s).attribute("what", "gaussian");
        h5.write("c", c);
        h5.write("gauge/u", u);
        h5.write("gauge/u_float", u, hila::io::to_float);
        h5.write("ints/k", k);
        h5.write("ints/cv", cv);
        h5.write("arrays/v", v, {5, 3});
        h5.write("arrays/vc", vc, {4});
        h5.group("empty/sub").attribute("note", "empty group");
        h5.group("gauge").attribute("group attr", 3.25f);
        h5.write("arrays/" + std::string(300, 'L'), v, {15}, hila::io::to_float);
        // closed at end of scope
    }

    hila::finishrun();
}
