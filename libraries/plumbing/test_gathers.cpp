//////////////////////////////////////////////////////////////////////////////
/// Test the standard gather here
//////////////////////////////////////////////////////////////////////////////

#include "hila.h"

void gather_test() {

    int64_t s = 0;
    onsites(ALL) {
        s += 1;
    }

    if (s != lattice.volume()) {
        hila::out0 << " Reduction test error!  Sum " << s << " should be "
                << lattice.volume() << '\n';
        hila::terminate(1);
    }


    foralldir (d) {

        CoordinateVector dif1 = 0, dif2 = 0;
        Field<CoordinateVector> f1, f2;
#pragma hila novector
        onsites(ALL) {
            f1[X] = X.coordinates();
            f2[X] = (X.coordinates() + d).mod(lattice.size());
        }

        onsites(ALL) {
            dif1 += abs(f1[X + d] - f2[X]);
            dif2 += abs(f1[X] - f2[X - d]);
        }

        if (dif1.squarenorm() != 0) {
            hila::out0 << " Std up-gather test error! Node " << hila::myrank()
                    << " direction " << (unsigned)d << " dif1 " << dif1 << '\n';
            hila::terminate(1);
        }

        if (dif2.squarenorm() != 0) {
            hila::out0 << " Std down-gather test error! Node " << hila::myrank()
                    << " direction " << (unsigned)d << " dif2 " << dif2 << '\n';
            hila::terminate(1);
        }

#if 0 && defined(SPECIAL_BOUNDARY_CONDITIONS)
        // test antiperiodic b.c. to one direction
        if (next_direction(d) == NDIM) {
            f2.set_boundary_condition(d, hila::bc::ANTIPERIODIC);

            onsites(ALL) {
                if (X.coordinate(d) == lattice.size(d) - 1)
                    f2[X] = -f1[X];
                else
                    f2[X] = f1[X];
            }

            dif1 = 0;
            onsites(ALL) { dif1 += f1[X] - f2[X - d]; }

            if (dif1 != 0) {
                hila::out0 << " Antiperiodic up-gather test error! Node " << hila::myrank()
                        << " direction " << (unsigned)d << '\n';
                hila::terminate(1);
            }
        }
#endif
    }
}

void test_std_gathers() {
    // gather_test<int>();
    gather_test();

#if defined(CUDA) || defined(HIP)
    gpuMemPoolPurge();
#endif

    hila::timestamp("Communication tests done");
    print_dashed_line();

    if (hila::myrank() == 0)
        hila::out.flush();
}
