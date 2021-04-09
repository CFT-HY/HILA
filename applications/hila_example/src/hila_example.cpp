#include "hila.h"
#include "fft.h"

// define an alias for a 3x3 complex matrix
using Mtype = SquareMatrix<3,Complex<double>>;

int main(int argc, char **argv) {

    // hila::initialize should be called as early as possible
    hila::initialize(argc, argv);


    // hila provides an input class hila::input, which is 
    // a convenient way to read in parameters from input files.
    // parameters are presented as key - value pairs, as an example
    //  " lattice size  64, 64, 64 "
    // is read below.
    //
    // Values are broadcast to all MPI nodes.
    //
    // .get() -method can read many different input types, 
    // see file "input.h" for documentation

    hila::input par("parameters");

    CoordinateVector lsize;
    lsize              = par.get("lattice size");    // reads NDIM numbers
    int loops          = par.get("smear loops");
    double smear_coeff = par.get("smear coefficient");
    int taylor_order   = par.get("expansion order");

    long seed          = par.get("random seed");

    par.close();    // file is closed also when par goes out of scope

    // setting up the lattice is convenient to do after reading
    // the parameter
    lattice->setup(lsize);

    // We need random number here
    hila::seed_random(seed);



    // 2 matrix fields
    Field<Mtype> f, g;

    // set g to gaussian rnd matrix
    onsites(ALL) g[X].gaussian();

    output0 << "Smearing a Gaussian random "
            << Mtype::rows() << "x" << Mtype::columns() 
            << " complex matrix field " << loops << " times\n";


    double c1 = 1 - 6*smear_coeff;

    // Use a timer to time periodic events
    // Automatically reported at the end
    // Good idea to define these static
    static hila::timer smear_timer("Smear");

    for (int l=0; l<loops; l++) {
        smear_timer.start();
        onsites(ALL) {
            f[X] = g[X+e_x] + g[X-e_x] +
                   g[X+e_y] + g[X-e_y] +
                   g[X+e_z] + g[X-e_z];
        }
        g[ALL] = c1*g[X] + smear_coeff*f[X];
        smear_timer.stop();
    }

    output0 << "field g at (0,0,0) after smearing:\n";
    output0 << g[{0,0,0}] << '\n';


    output0 << "Calculating exp(g) using Taylor expansin to order "
            << taylor_order << '\n';

    // another way to time, using gettime
    double t = hila::gettime();

    onsites(ALL) {
        Mtype product;
        product = 1;
        f[X] = 1;
        int64_t fac = 1;
        for (int i=1; i<=taylor_order; i++) {
            product *= g[X];
            fac *= i;
            f[X] += product/fac;
        }
    }

    output0 << "Taylor expansion, time " << hila::gettime() - t << " seconds\n";

    output0 << "Result at (0,0,0):\n";
    output0 << f[{0,0,0}] << '\n';


    t = hila::gettime();

    FFT_field(f, g);

    output0 << "FFT time " << hila::gettime() - t << '\n';


    hila::finishrun();
}
