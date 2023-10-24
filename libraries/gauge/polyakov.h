/** @file polyakov.h */

#ifndef POLYAKOV_H_
#define POLYAKOV_H_

#include "hila.h"

/**
 * @brief Measure Polyakov lines to direction dir
 * @details Naive implementation, includes extra communication
 * @tparam T GaugeField Group
 * @param U GaugeField to measure
 * @param dir Direction
 * @return Complex<double>
 */
template <typename T>
Complex<double> measure_polyakov(const GaugeField<T> &U, Direction dir = Direction(NDIM - 1)) {

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

// template <typename T>
// Complex<double> measure_polyakov_twist(const GaugeField<T> &U, Direction dir = Direction(NDIM - 1)) {

//     Field<T> polyakov = U[dir];

//     // mult links so that polyakov[X.dir == 0] contains the polyakov loop
//     for (int plane = lattice.size(dir) - 2; plane >= 0; plane--) {

//         // safe_access(polyakov) pragma allows the expression below, otherwise
//         // hilapp would reject it because X and X+dir can refer to the same
//         // site on different "iterations" of the loop.  However, here this
//         // is restricted on single dir-plane so it works but we must tell it to hilapp.

// #pragma hila safe_access(polyakov)
//         onsites(ALL) {
//             if (X.coordinate(dir) == plane) {
//                 polyakov[X] = U[dir][X] * polyakov[X + dir];
//             }
//         }
//     }

//     Complex<double> ploop = 0;
//     ReductionVector<Complex<double>> ploop_z(lattice.size(e_z) + 1);
//     ploop_z = 0;
//     onsites(ALL) if (X.coordinate(dir) == 0) {
//         Complex<double> p = trace(polyakov[X]);
//         ploop += p;
//         //ploop_z[X.z()] += p;
//     }
//     ploop_z[lattice.size(e_z)] = ploop / (lattice.volume() / lattice.size(dir));
//     for (int i = 0; i < ploop_z.size() - 1; i++) {
//         ploop_z[i] /= (lattice.volume() / lattice.size(dir));
//     }
//     return ploop_z.vector();
// //    return ploop;
// }

#endif