#include "hila.h"
//#include <roctx.h>
//#include <roctracer_ext.h>
int main(int argc, char *argv[]) {


    hila::initialize(argc, argv);
    hila::input par("parameters-matrix-test");
    CoordinateVector lsize;

    lsize = par.get("volume"); // reads NDIM numbers

    int N = par.get("iter");

    lattice.setup(lsize);
    // Random numbers are used here - use time to seed
    hila::seed_random(0);

    Field<SU<2, double>> f_2;
    // Field<SU<2, double>> d_2;
    // Field<SU<8, float>> f_8;
    // Field<SU<8, double>> d_8;
    hila::timer f_2_timer("flaot SU(8)");
    // hila::timer d_2_timer("double SU(2)");
    // hila::timer f_8_timer("float SU(8)");
    // hila::timer d_8_timer("double SU(8)");

    f_2 = 1;
    // d_2 = 1;
    // f_8 = 1;
    // d_8 = 1;

    double reduce;
    
    f_2_timer.start();
    // for (int i=0; i<N; i++) {
    //     foralldir(d) {
    //         onsites(EVEN) {
    //             f_2[X] = f_2[X + d]*f_2[X - d];
    //             reduce += trace(f_2[X]).real();
    //         }
    //         onsites(ODD) {
    //             f_2[X] = f_2[X + d]*f_2[X - d];
    //             reduce += trace(f_2[X]).real();
    //         }
    //     }
    // }

    for (int i=0; i<N; i++) {
        onsites(EVEN) {
            f_2[X] = f_2[X + e_z]*f_2[X - e_z];
        }
        onsites(ODD) {
            f_2[X] = f_2[X + e_z]*f_2[X - e_z];
        }
    }

    hila::synchronize();
    f_2_timer.stop();

    Reduction<Complex<double>> result;
    result.allreduce(false);
    onsites(ALL) {
        result += trace(f_2[X]);
    }    
    hila::out0 << f_2[{0,0,0}] << std::endl;
    hila::out0 << "Result is " << result.value()/lattice.volume() << "\n";
    hila::finishrun();
    return (0);
}