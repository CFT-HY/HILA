#include "hila.h"

#include <random>

static_assert(NDIM == 3, "NDIM must be 3 here");


using MyType = Complex<double>;

int main(int argc, char *argv[]) {

    // initialize system
    hila::initialize(argc, argv);

    // set up 32^3 lattice
    lattice->setup({256, 256, 256});

    // Random numbers are used here - use time to seed
    hila::seed_random(0);

    // lattice field
    Field<double> g = 2.0;

    g[{1, 1, 1}] = 0.9;
    g[{4, 0, 0}] = 1.0;
    g[{1, 1, 31}] = 2.1;
    g[{24,0, 2}] = 2.3;

    double val1, val2, val3, val4, val5, val6;
    CoordinateVector loc1, loc2, loc3, loc4, loc5, loc6;
    for (auto i = 0; i < 1000; i++)
    {
        val1 = g.min(ODD, loc1);
        val2 = g.min(EVEN,loc2);
        val3 = g.max(ODD, loc3);
        val4 = g.max(EVEN,loc4);
        val5 = g.min(ALL, loc5);
        val6 = g.max(ALL, loc6);
    }
    
    // val = g.max(loc);
    output0 << "Min value of ODD sites " << val1 << " at location: " << loc1 << '\n';
    output0 << "Min value of EVEN sites  " << val2 << " at location: " << loc2 << '\n';
    output0 << "Max value of ODD sites " << val3 << " at location: " << loc3 << '\n';
    output0 << "Max value of EVEN sites  " << val4 << " at location: " << loc4 << '\n';
    output0 << "Min value of ALL sites  " << val5 << " at location: " << loc5 << '\n';
    output0 << "Max value of ALL sites  " << val6 << " at location: " << loc6 << '\n';


    //output0 << "Reduction test " << reduced << "\n";
    hila::finishrun();
    return 0;
}
