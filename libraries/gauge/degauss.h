#ifndef DEGAUSS_H_
#define DEGAUSS_H_

#include "hila.h"

/////////////////////////////////////////////////////////////////////////////
/// Remove Gauss' law violation from a gauge field + E-field - assuming no charged
/// fields
/// Call:
///  degauss( const GaugeField<grp> &U, VectorField<Algebra<grp>> &E, double tolerance);
/// tolerance is the max violation / site
// Gauss' law violation is
//
//  G^a(x) = D_i E_i^a
//
// Restoration is done by evolving the fields 'downhill' towards
// G == 0


template <typename group>
void get_gauss_violation(const GaugeField<group> &U,
                         const VectorField<Algebra<group>> &E,
                         Field<Algebra<group>> &gauss) {


    // Measure Gauss
    gauss = 0;
    foralldir(d) {
        gauss[ALL] +=
            E[d][X] - (U[d][X - d].dagger() * E[d][X - d].expand() * U[d][X - d])
                          .project_to_algebra();
    }
}


template <typename group>
void gauss_fix_step(const GaugeField<group> &U, VectorField<Algebra<group>> &E,
                    const Field<Algebra<group>> &violation,
                    const hila::scalar_type<group> mag) {

    foralldir(dir) {
        onsites(ALL) {
            // get violation from "up"
            auto ta = (U[dir][X] * violation[X + dir].expand() * U[dir][X].dagger())
                          .project_to_algebra();
            E[dir][X] -= mag * (violation[X] - ta);
        }
    }
}


template <typename group>
auto degauss_step(const GaugeField<group> &U, VectorField<Algebra<group>> &E) {

    Field<Algebra<group>> violation;

    //  Use here the trick of 2 step sizes, the smaller one guarantees
    //  that the UV is stable and the larger one makes gradient flow
    //  faster.

    constexpr double aleph = 1.25;
    constexpr hila::scalar_type<group> onealeph = aleph / 12.0;
    constexpr hila::scalar_type<group> twoaleph = 2.0 * onealeph;

    //  Now what's the initial violation
    get_gauss_violation(U, E, violation);
    // Apply as a correction
    gauss_fix_step(U, E, violation, onealeph);

    // and repeat with twoaleph
    get_gauss_violation(U, E, violation);
    gauss_fix_step(U, E, violation, twoaleph);

    // we return not the last violation, thus the real violation
    // should be smaller than this
    // Return violation per site
    return violation.squarenorm() / lattice.volume();
}


/////////////////////////////////////////////////////////////////////
/// De-Gauss the gauge + E(momenta) -system.  quality is the violation/site
///
template <typename group>
int degauss(const GaugeField<group> &U, VectorField<Algebra<group>> &E,
            double quality) {

    constexpr int loop_max = 10000;
    double viol;
    int loop;

    static hila::timer degauss_timer("Degauss");

    /* loop over Degauss_step until violation small enough */

    degauss_timer.start();

    loop = 0;
    do {
        loop++;
        viol = degauss_step(U, E);
    } while (viol > quality && loop <= loop_max);

    if (loop > loop_max && hila::myrank() == 0) {
        hila::out << " ********** LOOP_MAX " << loop_max
                     << " reached in degauss, violation/site " << viol << '\n';
    }

    degauss_timer.stop();

    return loop;
}


#endif