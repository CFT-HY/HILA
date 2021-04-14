//////////////////////////////////////////////////////////////////////////////
/// Test the standard gather here
//////////////////////////////////////////////////////////////////////////////

#include "defs.h"
#include "coordinates.h"
#include "lattice.h"

#include "field.h"
#include "datatypes/matrix.h"

template <typename T> void gather_test() {

    foralldir(d) {

        T dif1 = 0, dif2 = 0;
        Field<T> f1, f2;
        onsites(ALL) {
            f1[X] = X.coordinate(d);
            f2[X] = mod(X.coordinate(d) + 1, lattice->size(d));
        }

        onsites(ALL) {
            dif1 += f1[X + d] - f2[X];
            dif2 += f1[X] - f2[X - d];
        }

        if (dif1 != 0) {
            output0 << " Std up-gather test error! Node " << hila::myrank()
                    << " direction " << (unsigned)d << '\n';
            exit(1);
        }
        if (dif2 != 0) {
            output0 << " Std down-gather test error! Node " << hila::myrank()
                    << " direction " << (unsigned)d << '\n';
            exit(1);
        }

#if 0 && defined(SPECIAL_BOUNDARY_CONDITIONS)
        // test antiperiodic b.c. to one direction
        if (next_direction(d) == NDIM) {
            f2.set_boundary_condition(d, BoundaryCondition::ANTIPERIODIC);

            onsites(ALL) {
                if (X.coordinate(d) == lattice->size(d) - 1)
                    f2[X] = -f1[X];
                else
                    f2[X] = f1[X];
            }

            dif1 = 0;
            onsites(ALL) { dif1 += f1[X] - f2[X - d]; }

            if (dif1 != 0) {
                output0 << " Antiperiodic up-gather test error! Node " << hila::myrank()
                        << " direction " << (unsigned)d << '\n';
                exit(1);
            }
        }
#endif
    }
}

void test_std_gathers() {
    // gather_test<int>();
    gather_test<double>();

    hila::timestamp("Communication tests done");
    print_dashed_line();

    if (hila::myrank() == 0)
        hila::output.flush();
}
