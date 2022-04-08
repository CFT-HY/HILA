#ifndef STAPLESUM_H_
#define STAPLESUM_H_

#include "hila.h"
#include "vector_field.h"

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

    // zero staples just inc case
    staples = 0;
    foralldir(d2) if (d2 != d1) {

        // anticipate that these are needed
        // not really necessary, but may be faster
        U[d2].start_gather(d1, ALL);
        U[d1].start_gather(d2, par);

        // calculate first lower 'U' of the staple sum
        // do it on opp parity
        onsites(~par) {
            lower[X] = U[d2][X].dagger() * U[d1][X] * U[d2][X + d1];
        }

        // calculate then the upper 'n', and add the lower
        // lower could also be added on a separate loop 
        onsites(par) {
            staples[X] +=
                U[d2][X] * U[d1][X + d2] * U[d2][X + d1].dagger() + lower[X - d2];
        }
    }
}

#endif