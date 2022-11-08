#ifndef POLYAKOV_H_
#define POLYAKOV_H_

#include "hila.h"

/////////////////////////////////////////////////////////////////////////////
/// Measure Polyakov lines to direction dir
///
/// Naive implementation, includes extra communications

template <typename T>
Complex<double> measure_polyakov(const GaugeField<T> &U, Direction dir) {

    Field<T> polyakov = U[dir];

    // mult links so that polyakov[X.dir == 0] contains the polyakov loop
    for (int plane = lattice.size(dir) - 2; plane >= 0; plane--) {

        // safe_access(polyakov) pragma allows the expression below, otherwise
        // hilapp would reject it because X and X+dir can refer to the same
        // site on different "iterations" of the loop.  However, here this
        // is restricted on single dir-plane so it works but we must tell it to hilapp.

#pragma hila safe_access(polyakov)
        onsites(ALL) {
            if (X.coordinate(dir) == plane) {
                polyakov[X] = U[dir][X] * polyakov[X + dir];
            }
        }
    }

    Complex<double> ploop = 0;

    onsites(ALL) if (X.coordinate(dir) == 0) {
        ploop += trace(polyakov[X]);
    }

    // return average polyakov
    return ploop / (lattice.volume() / lattice.size(dir));
}

#endif