#ifndef STOUT_SMEAR_NCHM_H_
#define STOUT_SMEAR_NCHM_H_

#include "hila.h"

#include "gauge/wilson_line_and_force.h"

/////////////////////////////////////////////////////////////////////////////
/// Do stout smearing for the gauge fields
///  U' = exp( -coeff * Proj_to_Alg(U * Staples) ) * U

/**
 * @brief Compute sum of staples around all links
 * @details Compute staple sums around all positively oriented (i.e. corresponding
 * plaquette sums are obtained by multiplying the staple sums by the corresponding
 * positively oriented link variables)

 * @tparam T Matrix element type
 * @param U input gauge field of type T
 * @param staples output vector/gauge field of type T
 * @return void
 */
template <typename T>
void nchm_staplesums(const GaugeField<T> &U, out_only GaugeField<T> &staples) {
    Field<T> tF1, lw1;
    foralldir(d1) {
        foralldir(d2) if (d2 != d1) {
            U[d2].start_gather(d1, ALL);
        }
        onsites(ALL) staples[d1][X] = 0;
    }

    foralldir(d1) foralldir(d2) if (d2 != d1) {
        onsites(ALL) {
            // tF1[X] = (U[d1][X] * U[d2][X + d1]).dagger();
            mult_aa(U[d1][X], U[d2][X + d1], tF1[X]);
            // lw1[X] = tF1[X] * U[d2][X];
            mult(tF1[X], U[d2][X], lw1[X]);
        }

        lw1.start_gather(-d2, ALL);

        onsites(ALL) {
            // staple[d2][X] += U[d1][X + d2] * tF1[X];
            mult_add(U[d1][X + d2], tF1[X], staples[d2][X]);
        }
        onsites(ALL) {
            staples[d1][X] += lw1[X - d2];
        }
    }
}

template <typename T, typename atype = hila::arithmetic_type<T>>
void nchm_stout_smear1(const GaugeField<T> &U, out_only GaugeField<T> &stout,
                       out_only std::array<Field<int>, NDIM> &niter, atype coeff) {
    nchm_staplesums(U, stout);
    foralldir(d) {
        onsites(ALL) {
            stout[d][X] = altexp((U[d][X] * stout[d][X]).project_to_algebra_scaled(-coeff).expand(),
                                niter[d][X]) * U[d][X];
        }
    }
}

template <typename T, typename atype = hila::arithmetic_type<T>>
void nchm_stout_smear1(const GaugeField<T> &U, out_only GaugeField<T> &stout,
                       out_only VectorField<T> &stap,
                       out_only std::array<Field<int>, NDIM> &niter, atype coeff) {
    nchm_staplesums(U, stap);
    foralldir(d) {
        onsites(ALL) {
            stout[d][X] = altexp((U[d][X] * stap[d][X]).project_to_algebra_scaled(-coeff).expand(),
                              niter[d][X]) *
                          U[d][X];
        }
    }
}

template <typename T, typename atype = hila::arithmetic_type<T>>
void nchm_stout_smear(const GaugeField<T> &U, out_only std::vector<GaugeField<T>> &stoutlist,
                      out_only std::vector<GaugeField<T>> &staplist,
                      out_only std::vector<std::array<Field<int>, NDIM>> &niterlist, atype coeff) {
    // performs nst stout smearing steps on the gauge field U and stores the gague field obtained
    // after the i-th smearing step in stoutlist[i]
    stoutlist[0] = U; // original field
    for (int i = 1; i < stoutlist.size(); ++i) {
        nchm_stout_smear1(stoutlist[i - 1], stoutlist[i], staplist[i - 1], niterlist[i - 1], coeff);
    }
}

template <typename T, typename atype = hila::arithmetic_type<T>>
void nchm_stout_smear1(const GaugeField<T> &U, out_only GaugeField<T> &stout, atype coeff) {
    nchm_staplesums(U, stout);
    foralldir(d) {
        onsites(ALL) {
            stout[d][X] =
                altexp((U[d][X] * stout[d][X]).project_to_algebra_scaled(-coeff).expand()) *
                U[d][X];
        }
    }
}

template <typename T, typename atype = hila::arithmetic_type<T>>
void nchm_stout_smear1(const GaugeField<T> &U, out_only GaugeField<T> &stout,
                       out_only GaugeField<T> &stap, atype coeff) {
    nchm_staplesums(U, stap);
    foralldir(d) {
        onsites(ALL) {
            stout[d][X] =
                altexp((U[d][X] * stap[d][X]).project_to_algebra_scaled(-coeff).expand()) * U[d][X];
        }
    }
}

template <typename T, typename atype = hila::arithmetic_type<T>>
void nchm_stout_smear(const GaugeField<T> &U, GaugeField<T> &stout, atype coeff, int iter) {
    GaugeField<T> tmp;
    if (iter % 2 == 0) {
        stout = U;
        for (int i = 0; i < iter / 2; ++i) {
            nchm_stout_smear1(stout, tmp, coeff);
            nchm_stout_smear1(tmp, stout, coeff);
        }
    } else {
        tmp = U;
        for (int i = 0; i < iter / 2; ++i) {
            nchm_stout_smear1(tmp, stout, coeff);
            nchm_stout_smear1(stout, tmp, coeff);
        }
        nchm_stout_smear1(tmp, stout, coeff);
    }
}

template <typename T, typename atype = hila::arithmetic_type<T>>
void nchm_stout_smear(const GaugeField<T> &U, GaugeField<T> &stout,
                      std::vector<std::array<Field<int>, NDIM>> &niterlist, atype coeff) {
    GaugeField<T> tmp;
    int iter = niterlist.size();
    if (iter % 2 == 0) {
        stout = U;
        for (int i = 0; i < iter / 2; ++i) {
            nchm_stout_smear1(stout, tmp, niterlist[2 * i], coeff);
            nchm_stout_smear1(tmp, stout, niterlist[2 * i + 1], coeff);
        }
    } else {
        tmp = U;
        for (int i = 0; i < iter / 2; ++i) {
            nchm_stout_smear1(tmp, stout, niterlist[2 * i], coeff);
            nchm_stout_smear1(stout, tmp, niterlist[2 * i + 1], coeff);
        }
        nchm_stout_smear1(tmp, stout, niterlist[iter - 1], coeff);
    }
}


template <typename T, int N = T::rows(), int NSTP = 2 * (NDIM - 1),
          typename atype = hila::arithmetic_type<T>>
void nchm_stout_smear_force(const std::vector<GaugeField<T>> &stoutlist,
                            const std::vector<VectorField<T>> &staplist,
                            const VectorField<Algebra<T>> &K, out_only VectorField<Algebra<T>> &KS,
                            std::vector<std::array<Field<int>, NDIM>> &niterlist, atype coeff) {
    // uses the list stoutlist[] of smeared gauge fields to compute the pullback of the
    // algebra-valued force field K under the smearing and returns the force acting on the unsmeared
    // link variables as algebra-valued field KS
    // Note: our definition of the force field is different from the one used in
    // [arXiv:hep-lat/0311018v1], in order to match the force field representation used by our HMC
    // implementation.
    VectorField<T> K1;
    Field<T> K21, K22, K23, K24;

    foralldir(d1) onsites(ALL) KS[d1][X] = K[d1][X];

    for (int i = stoutlist.size() - 2; i >= 0; --i) {
        const GaugeField<T> &U0 = stoutlist[i + 1];
        const GaugeField<T> &U = stoutlist[i];
        const VectorField<T> &staps = staplist[i];
        const std::array<Field<int>, NDIM> &niter = niterlist[i];

        foralldir(d1) {
            //get_stout_staples(U, d1, stapl, staps);
            onsites(ALL) {
                // compute stout smearing operator and its derivatives:

                // temp. variables:
                T mtexp;
                T mdtexp;
                // turn staple sum into plaquette sum by multiplying with link variable:
                //T tplaqs = U[d1][X] * staps[X];
                T tplaqs = U[d1][X] * staps[d1][X];
                // the following function computes first for X = -coeff * tplaqs the smearing
                // operator Q = exp(X) and its derivatives dQ/dX[][], and uses these to
                // compute the two matrices:
                // mtexp = Q.dagger() * KS[d1][X].expand() * Q
                // and
                // mdtexp[i][j] = trace(Q.dagger * KS[d1][X].expand() * dQ/dX[j][i]) :
                mult_exp(tplaqs.project_to_algebra_scaled(-coeff).expand(),
                         U[d1][X] * U0[d1][X].dagger() * KS[d1][X].expand(), mtexp, mdtexp, niter[d1][X]);

                // set K1[d1][X] to be the equivalent of the \Lambda matrix from eq.(73) in
                // [arXiv:hep-lat/0311018v1]:
                K1[d1][X] = mdtexp.project_to_algebra_scaled(coeff).expand();

                // equivalent of first line and first term on second line of eq.(75) in
                // [arXiv:hep-lat/0311018v1]:
                KS[d1][X] = (mtexp - tplaqs * K1[d1][X]).project_to_algebra();

                // multiply K1[d1] by U[d1]:
                K1[d1][X] *= U[d1][X];
            }
        }

        // equivalent of remaining terms of eq.(75) in [arXiv:hep-lat/0311018v1]:
        foralldir(d1) foralldir(d2) if (d1 != d2) {
            onsites(ALL) {
                T U2, U4, tM1;

                U2 = U[d2][X + d1];
                U4 = U[d2][X].dagger();

                tM1 = U2 * U[d1][X + d2].dagger();
                K21[X] = U4 * K1[d1][X] * tM1;
                tM1 *= U4;
                K22[X] = tM1 * K1[d1][X];
                K24[X] = K1[d1][X] * tM1;

                tM1 = U2 * K1[d1][X + d2].dagger() * U4;
                K23[X] = U[d1][X] * tM1;
                K22[X] += tM1 * U[d1][X];

                K24[X] += K23[X];
            }

            onsites(ALL) {
                KS[d2][X] -= (K22[X - d1] - K24[X]).project_to_algebra();
                KS[d1][X] -= (K23[X] - K21[X - d2]).project_to_algebra();
            }
        }
    }
}

#endif