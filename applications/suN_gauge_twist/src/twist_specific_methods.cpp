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
void staplesum_twist(const GaugeField<T> &U, Field<T> &staples, Direction d1, Parity par = ALL,
                     int twist_coeff = 1) {

    Field<T> lower;
    GaugeField<double> twist = 0;
    onsites(ALL) {
        if (X.z() == 0 && X.t() == 0) {
            twist[e_z][X] = twist_coeff;
            twist[e_t][X] = -twist_coeff;
        }
    }

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
                staples[X] = U[d2][X] * U[d1][X + d2] * U[d2][X + d1].dagger() *
                                 expi(2 * M_PI * (twist[d1][X] / NCOLOR)) +
                             lower[X - d2] * expi(2 * M_PI * (twist[d1][X] / NCOLOR)).conj();
            }
            first = false;
        } else {
            onsites(par) {
                staples[X] += U[d2][X] * U[d1][X + d2] * U[d2][X + d1].dagger() *
                                  expi(2 * M_PI * (twist[d1][X] / NCOLOR)) +
                              lower[X - d2] * expi(2 * M_PI * (twist[d1][X] / NCOLOR)).conj();
            }
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
double measure_plaq_twist(GaugeField<T> U, int twist_coeff = 1) {
    Reduction<double> plaq;
    plaq.allreduce(false);

    foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {

        onsites(ALL) {
            if (dir1 == e_z && dir2 == e_t && X.z() == 0 && X.t() == 0) {
                plaq += 1.0 -
                        real(trace(U[dir1][X] * U[dir2][X + dir1] * U[dir1][X + dir2].dagger() *
                                   U[dir2][X].dagger() * expi(2 * M_PI * (twist_coeff / NCOLOR)))) /
                            T::size();
            } else {
                plaq += 1.0 - real(trace(U[dir1][X] * U[dir2][X + dir1] *
                                         U[dir1][X + dir2].dagger() * U[dir2][X].dagger())) /
                                  T::size();
            }
        }
    }

    return plaq.value();
}
