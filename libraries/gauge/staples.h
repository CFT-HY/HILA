#ifndef STAPLESUM_H_
#define STAPLESUM_H_

#include "hila.h"

/////////////////////////////////////////////////////////////////////////////
/// Sum the staples of link matrices to direction dir
/// we could do
/// foralldir(d2) if (d2 != d1)
///    stapes[par] += U[d2][X]*U[d1][X+d2]*U[d2][X+d1].dagger()  +
///                   U[d2][X-d2].dagger()*U[d1][X-d2]*U[d2][X-d2+d1]
/// but let us do it a bit more optimized way

template <typename T>
void staplesum(const GaugeField<T> &U, Field<T> &staples, Direction d1, Parity par = ALL) {

    Field<T> lower;

    bool first = true;
    foralldir(d2) if (d2 != d1) {

        // anticipate that these are needed
        // not really necessary, but may be faster
        U[d2].start_gather(d1, ALL);
        U[d1].start_gather(d2, par);

        // calculate first lower 'U' of the staple sum
        // do it on opp parity
        onsites(opp_parity(par)) {
            lower[X] = U[d2][X].dagger() * U[d1][X] * U[d2][X + d1];
        }

        // calculate then the upper 'n', and add the lower
        // lower could also be added on a separate loop
        if (first) {
            onsites(par) {
                staples[X] = U[d2][X] * U[d1][X + d2] * U[d2][X + d1].dagger() + lower[X - d2];
            }
            first = false;
        } else {
            onsites(par) {
                staples[X] += U[d2][X] * U[d1][X + d2] * U[d2][X + d1].dagger() + lower[X - d2];
            }
        }
    }
}

template <typename T>
void staplesum_db(const GaugeField<T> &U, Field<T> &staples, Direction d1, Parity par,
                  double deltab) {

    Field<T> lower;

    // zero staples just in case
    staples[par] = 0;

    foralldir(d2) if (d2 != d1) {

        // anticipate that these are needed
        // not really necessary, but may be faster
        U[d2].start_gather(d1, ALL);
        U[d1].start_gather(d2, par);

        // calculate first lower 'U' of the staple sum
        // do it on opp parity
        onsites(opp_parity(par)) {
            double m;
            if (2 * X.z() < lattice.size(e_z))
                m = 1.0 - deltab;
            else
                m = 1.0 + deltab;

            lower[X] = m * U[d2][X].dagger() * U[d1][X] * U[d2][X + d1];
        }

        // calculate then the upper 'n', and add the lower
        // lower could also be added on a separate loop
        onsites(par) {
            double m;
            if (2 * X.z() < lattice.size(e_z))
                m = 1.0 - deltab;
            else
                m = 1.0 + deltab;

            staples[X] += m * U[d2][X] * U[d1][X + d2] * U[d2][X + d1].dagger() + lower[X - d2];
        }
    }
}

#endif