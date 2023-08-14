#ifndef TWIST_SPECIFIC_METHODS_H_
#define TWIST_SPECIFIC_METHODS_H_

#include "hila.h"

/**
 * @brief Sum the staples of link matrices to direction dir by inducing twist at z=0 t=0 on x-y
 * plane
 *
 * @tparam T
 * @param U GaugeField to compute staples for
 * @param staples Filed to compute staplesum into at each lattice point
 * @param d1 Direction to compute staplesum for
 * @param par Parity to compute staplesum for
 * @param twist_coeff Integer to rotate phase with
 */
template <typename T>
void staplesum_twist(const GaugeField<T> &U, Field<T> &staples, Direction d1, int twist_coeff,
                     Parity par = ALL) {

    Field<T> lower;
    // GaugeField<double> twist = 0;
    // onsites(ALL) {
    //     if (X.z() == 0 && X.t() == 0) {
    //         twist[e_z][X] = twist_coeff;
    //         twist[e_t][X] = -twist_coeff;
    //     }
    // }

    bool first = true;
    staples = 0;

    foralldir(d2) if (d2 != d1) {

        // anticipate that these are needed
        // not really necessary, but may be faster
        U[d2].start_gather(d1, ALL);
        U[d1].start_gather(d2, par);

        double twist;
        if (d1 == e_z && d2 == e_t)
            twist = twist_coeff;
        else if (d1 == e_t && d2 == e_z)
            twist = -twist_coeff;
        else
            twist = 0;

        twist /= NCOLOR;

        // calculate first lower 'U' of the staple sum
        // do it on opp parity
        onsites(opp_parity(par)) {
            lower[X] = U[d2][X].dagger() * U[d1][X] * U[d2][X + d1];
            if (X.z() == 0 && X.t() == 0)
                lower[X] *= expi(-2 * M_PI * twist);
        }

        // calculate then the upper 'n', and add the lower
        // lower could also be added on a separate loop
        onsites(par) {
            auto upper = U[d2][X] * U[d1][X + d2] * U[d2][X + d1].dagger();
            if (X.z() == 0 && X.t() == 0)
                upper *= expi(2 * M_PI * twist);

            staples[X] += upper + lower[X - d2];
        }
    }
}

/**
 * @brief Computes Wilson action
 * @details \f{align}{ S &=  \beta\sum_{\textbf{dir}_1 < \textbf{dir}_2}\sum_{X} \frac{1}{N}
 * \Re\mathrm{Tr}\left[ 1- U_{\textbf{dir}_1 \textbf{dir}_2}(X) \right] \f} Where \f$\beta =
 * 2N/g^2\f$
 *
 * @return double
 */
template <typename T>
std::vector<double> measure_plaq_with_z(GaugeField<T> U, int twist_coeff) {
    Reduction<double> plaq;
    ReductionVector<double> plaq_vec(lattice.size(e_z) + 1);
    plaq_vec = 0;
    plaq.allreduce(false);
    plaq_vec.allreduce(false);
    GaugeField<double> twist = 0;
    onsites(ALL) {
        if (X.z() == 0 && X.t() == 0) {
            twist[e_z][X] = -twist_coeff;
        }
    }
    foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {

        onsites(ALL) {
            double p;
            p = 1.0 -
                real(expi(2 * M_PI * (twist[dir1][X] / NCOLOR))*trace(U[dir1][X] * U[dir2][X + dir1] * U[dir1][X + dir2].dagger() *
                           U[dir2][X].dagger() )) /
                    T::size();
            plaq += p;
            plaq_vec[X.z()] += p;
        }
    }
    plaq_vec[lattice.size(e_z)] = plaq.value() / (lattice.volume() * NDIM * (NDIM - 1) / 2);
    for (int i = 0; i < plaq_vec.size() - 1; i++) {
        plaq_vec[i] /= (lattice.volume() * NDIM * (NDIM - 1)) / 2;
    }

    return plaq_vec.vector();
}

/**
 * @brief Measure Polyakov lines to direction dir bining based on z index
 * @details Naive implementation, includes extra communication
 * @tparam T GaugeField Group
 * @param U GaugeField to measure
 * @param dir Direction
 * @return Complex<double>
 */
template <typename T>
std::vector<Complex<double>> measure_polyakov_twist(const GaugeField<T> &U, Direction dir = Direction(NDIM - 1)) {

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
    ReductionVector<Complex<double>> ploop_z(lattice.size(e_z) + 1);
    ploop_z = 0;
    onsites(ALL) if (X.coordinate(dir) == 0) {
        Complex<double> p = trace(polyakov[X]);
        ploop += p;
        ploop_z[X.z()] += p;
    }
    ploop_z[lattice.size(e_z)] = ploop / (lattice.volume() / lattice.size(dir));
    for (int i = 0; i < ploop_z.size() - 1; i++) {
        ploop_z[i] /= (lattice.volume() / lattice.size(dir));
    }
    //return average polyakov
    return ploop_z.vector();
    //return ploop;
}

#endif
