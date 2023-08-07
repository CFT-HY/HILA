#include <tuple>
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
                     int twist_coeff = 2) {

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
                             lower[X - d2] * expi(-2 * M_PI * (twist[d1][X] / NCOLOR));
            }
            first = false;
        } else {
            onsites(par) {
                staples[X] += U[d2][X] * U[d1][X + d2] * U[d2][X + d1].dagger() *
                                  expi(2 * M_PI * (twist[d1][X] / NCOLOR)) +
                              lower[X - d2] * expi(-2 * M_PI * (twist[d1][X] / NCOLOR));
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
// template <typename T>
// double measure_plaq_twist(GaugeField<T> U, int twist_coeff = 1) {
//     Reduction<double> plaq;
//     plaq.allreduce(false);

//     foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {

//         onsites(ALL) {
//             if (dir1 == e_z && X.z() == 0 && X.t() == 0) {
//                 plaq += 1.0 -
//                         real(trace(U[dir1][X] * U[dir2][X + dir1] * U[dir1][X + dir2].dagger() *
//                                    U[dir2][X].dagger() * expi(2 * M_PI * (twist_coeff /
//                                    NCOLOR)))) /
//                             T::size();
//             } else {
//                 plaq += 1.0 - real(trace(U[dir1][X] * U[dir2][X + dir1] *
//                                          U[dir1][X + dir2].dagger() * U[dir2][X].dagger())) /
//                                   T::size();
//             }
//         }
//     }

//     return plaq.value();
// }

// template <typename T>
// double measure_plaq_twist_better(GaugeField<T> U, int twist_coeff = 1) {
//     Reduction<double> plaq;
//     plaq.allreduce(false);
//     GaugeField<double> twist = 0;
//     onsites(ALL) {
//         if (X.z() == 0 && X.t() == 0) {
//             twist[e_z][X] = twist_coeff;
//             twist[e_t][X] = -twist_coeff;
//         }
//     }

//     foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {
//         onsites(ALL) {
//             plaq += 1.0 -
//                     real(trace(U[dir1][X] * U[dir2][X + dir1] * U[dir1][X + dir2].dagger() *
//                                 U[dir2][X].dagger() * expi(2 * M_PI * (twist[dir1][X] /
//                                 NCOLOR)))) /
//                         T::size();
//         }
//     }

//     return plaq.value();
// }

template <typename T>
std::vector<double> measure_plaq_with_z(GaugeField<T> U, int twist_coeff = 1) {
    Reduction<double> plaq;
    ReductionVector<double> plaq_vec(lattice.size(e_z) + 1);
    plaq.allreduce(false);
    plaq_vec.allreduce(false);

    foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {

        onsites(ALL) {
            double p;
            if (dir1 == e_z && X.z() == 0 && X.t() == 0) {
                p = 1.0 -
                    real(trace(U[dir1][X] * U[dir2][X + dir1] * U[dir1][X + dir2].dagger() *
                               U[dir2][X].dagger() * expi(2 * M_PI * (twist_coeff / NCOLOR)))) /
                        T::size();
                plaq += p;
                plaq_vec[X.z()] += p;
            } else {
                p = 1.0 - real(trace(U[dir1][X] * U[dir2][X + dir1] * U[dir1][X + dir2].dagger() *
                                     U[dir2][X].dagger())) /
                              T::size();
                plaq += p;
                plaq_vec[X.z()] += p;
            }
        }
    }
    plaq_vec[lattice.size(e_z)] = plaq.value() / (lattice.volume() * NDIM * (NDIM - 1) / 2);
    for (int i = 0; i < plaq_vec.size() - 1; i++) {
        plaq_vec[i] /= (lattice.volume() * NDIM * (NDIM - 1)) / 2;
    }

    return plaq_vec.vector();
}
//     Reduction<double> plaq;
//     //ReductionVector<double> p_z(lattice.size(e_z));
//     // p_z = 0;
//     plaq.allreduce(false);
//     // p_z.allreduce(false);
//     foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {

//         onsites(ALL) {
//             //double p;
//             if (dir1 == e_z && X.z() == 0 && X.t() == 0) {
//                 // p = 1.0 -
//                 //         real(trace(U[dir1][X] * U[dir2][X + dir1] * U[dir1][X + dir2].dagger()
//                 *
//                 //                    U[dir2][X].dagger() * expi(2 * M_PI * (twist_coeff /
//                 NCOLOR)))) /
//                 //             T::size();
//                 plaq += 1.0 -
//                         real(trace(U[dir1][X] * U[dir2][X + dir1] * U[dir1][X + dir2].dagger() *
//                                    U[dir2][X].dagger() * expi(2 * M_PI * (twist_coeff /
//                                    NCOLOR)))) /
//                             T::size();
//                 //p_z[X.z()] += p;

//             } else {
//                 // p = 1.0 - real(trace(U[dir1][X] * U[dir2][X + dir1] *
//                 //                          U[dir1][X + dir2].dagger() * U[dir2][X].dagger())) /
//                 //                   T::size();
//                 plaq += 1.0 -
//                         real(trace(U[dir1][X] * U[dir2][X + dir1] * U[dir1][X + dir2].dagger() *
//                                    U[dir2][X].dagger() * expi(2 * M_PI * (twist_coeff /
//                                    NCOLOR)))) /
//                             T::size();
//                 //p_z[X.z()] += p;
//             }
//         }
//     }

//     return plaq.value();
// }