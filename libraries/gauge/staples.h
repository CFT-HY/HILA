/** @file staples.h */

#ifndef STAPLESUM_H_
#define STAPLESUM_H_

#include "hila.h"

/**
 * @brief Sum the staples of link matrices to direction dir
 *
 * Naive method is to compute:
 *
 * \code {.cpp}
 * foralldir(d2) if (d2 != d1)
 *     stapes[par] += U[d2][X]*U[d1][X+d2]*U[d2][X+d1].dagger()  +
 *                    U[d2][X-d2].dagger()*U[d1][X-d2]*U[d2][X-d2+d1]
 * \endcode
 *
 * But the method is computed in a slightly more optimized way
 *
 * @tparam T
 * @param U GaugeField to compute staples for
 * @param staples Filed to compute staplesum into at each lattice point
 * @param d1 Direction to compute staplesum for
 * @param par Parity to compute staplesum for
 */
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
        hila::out0 << lower.get_value_at(0) << '\n';
        hila::out0 << "Direction d2: ";
        hila::out0 << d2 << std::endl;
        {
            auto v = lower.get_value_at(0);
            std::ostringstream _oss;
            _oss << v;
            std::string _s = _oss.str();
            if (_s.find("nan") != std::string::npos || _s.find("NaN") != std::string::npos) {
                hila::out0 << "FATAL: NaN detected in staples[0]; aborting\n";
                hila::finishrun();
                std::exit(EXIT_FAILURE);
            }
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
        hila::out0 << staples.get_value_at(0) << '\n';
        hila::out0 << "Direction d2: ";
        hila::out0 << d2 << std::endl;
        {
            auto v = staples.get_value_at(0);
            std::ostringstream _oss;
            _oss << v;
            std::string _s = _oss.str();
            if (_s.find("nan") != std::string::npos || _s.find("NaN") != std::string::npos) {
                hila::out0 << "FATAL: NaN detected in staples[0]; aborting\n";
                hila::finishrun();
                std::exit(EXIT_FAILURE);
            }
        }
    }
    hila::out0 << "Finished staplesum for direction d1: ";  
    hila::out0 << d1 << std::endl;
}

template <typename T>
using plaqw_t = std::array<std::array<Field<T>, NDIM>, NDIM>;
/**
 * @brief Sum the staples of link matrices to direction dir taking into account plaquette weights
 *
 * Naive method is to compute:
 *
 * \code {.cpp}
 * foralldir(d2) if (d2 != d1)
 *     stapes[par] += U[d2][X]*U[d1][X+d2]*U[d2][X+d1].dagger()  +
 *                    U[d2][X-d2].dagger()*U[d1][X-d2]*U[d2][X-d2+d1]
 * \endcode
 *
 * But the method is computed in a slightly more optimized way
 *
 * @tparam T
 * @param U GaugeField to compute staples for
 * @param staples Filed to compute staplesum into at each lattice point
 * @param d1 Direction to compute staplesum for
 * @param plaqw plaquette weights 
 * @param par Parity to compute staplesum for
 */
template <typename T>
void staplesum(const GaugeField<T> &U, Field<T> &staples, Direction d1,
               const plaqw_t<hila::arithmetic_type<T>> &plaqw, Parity par = ALL) {

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
            lower[X] = plaqw[d1][d2][X] * U[d2][X].dagger() * U[d1][X] * U[d2][X + d1];
        }

        // calculate then the upper 'n', and add the lower
        // lower could also be added on a separate loop
        if (first) {
            onsites(par) {
                staples[X] = plaqw[d1][d2][X] * U[d2][X] * U[d1][X + d2] * U[d2][X + d1].dagger() + lower[X - d2];
            }
            first = false;
        } else {
            onsites(par) {
                staples[X] += plaqw[d1][d2][X] * U[d2][X] * U[d1][X + d2] * U[d2][X + d1].dagger() + lower[X - d2];
            }
        }

    }

}

#endif