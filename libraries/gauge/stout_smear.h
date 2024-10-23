#ifndef STOUT_SMEAR_H_
#define STOUT_SMEAR_H_

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
void staplesums(const GaugeField<T> &U, out_only GaugeField<T> &staples) {
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
void stout_smear1(const GaugeField<T> &U, out_only GaugeField<T> &stout, atype coeff) {
    staplesums(U, stout);
    foralldir(d) {
        onsites(ALL) {
            stout[d][X] = chexp((U[d][X] * stout[d][X]).project_to_algebra_scaled(-coeff).expand()) * U[d][X];
        }
    }
}

/**
 * @brief performs iter stout smearing steps with smearing coefficient coeff on the gauge field U
 * and stores the resulting gague field in stout
 * @tparam T Matrix element type
 * @tparam atype arithmetic_type of T
 * @param U input gauge field of type T
 * @param stout Matrix of type T, containing resulting smeared gauge field
 * @param coeff atype number specifying the smearing coefficient
 * @param iter int number specifying the number of stout smearing steps to perform
 * @return void
 */
template <typename T, typename atype = hila::arithmetic_type<T>>
void stout_smear(const GaugeField<T> &U, GaugeField<T> &stout, atype coeff, int iter) {
    GaugeField<T> tmp;
    if(iter % 2 == 0) {
        stout = U;
        for (int i = 0; i < iter / 2; ++i) {
            stout_smear1(stout, tmp, coeff);
            stout_smear1(tmp, stout, coeff);
        }
    } else {
        tmp = U;
        for (int i = 0; i < iter / 2; ++i) {
            stout_smear1(tmp, stout, coeff);
            stout_smear1(stout, tmp, coeff);
        }
        stout_smear1(tmp, stout, coeff);
    }
}

/**
 * @brief Compute staples and staple sums around all links
 * @details Compute staples and their sums around all positively oriented (i.e.
 * corresponding plaquettes are obtained by multiplying the staple matrices by
 * their corresponding positively oriented missing link variables)

 * @tparam T Matrix element type
 * @tparam nstap number of staples around each link
 * @param U input gauge field of type T
 * @param S array of nstap vector/gauge fields of type T
 * @param Ssum vector/gauge field of type T, containing the staple sums
 * @return void
 */
template <typename T, const int nstap = 2 * (NDIM - 1)>
void get_stout_staples(const GaugeField<T> &U, VectorField<T>(out_only &S)[nstap],
                         out_only VectorField<T> &Ssum) {
    Field<T> tF1, lw1;
    foralldir(d1) {
        foralldir(d2) if (d2 != d1) {
            U[d2].start_gather(d1, ALL);
        }
    }
    int l, rl, k, rk;
    l = 0;
    foralldir(d1) {
        k = 0;
        foralldir(d2) if (d2 != d1) {
            onsites(ALL) {
                // tF1[X] = (U[d1][X] * U[d2][X + d1]).dagger();
                mult_aa(U[d1][X], U[d2][X + d1], tF1[X]);
                // lw1[X] = tF1[X] * U[d2][X];
                mult(tF1[X], U[d2][X], lw1[X]);
            }

            lw1.start_gather(-d2, ALL);
            if(d1<d2) {
                rl = l;
            } else {
                rl = l - 1;
            }
            onsites(ALL) {
                // S[rl][d2][X] = U[d1][X + d2] * tF1[X];
                mult(U[d1][X + d2], tF1[X], S[rl][d2][X]);
            }
            rk = nstap - 1 - k;
            onsites(ALL) {
                S[rk][d1][X] = lw1[X - d2];
            }
            ++k;
        }
        ++l;
    }
    foralldir(d1) {
        onsites(ALL) Ssum[d1][X] = S[0][d1][X];
        for (k = 1; k < nstap; ++k) {
            onsites(ALL) Ssum[d1][X] += S[k][d1][X];
        }
    }
}

template <typename T, typename atype = hila::arithmetic_type<T>>
void stout_smear(const GaugeField<T> &U, out_only std::vector<GaugeField<T>> &stoutlist,
                 atype coeff) {
    // performs as many stout smearing steps on the gauge field U as necessary to fill stoutlist
    // storing the gague field obtained after the i-th smearing step in stoutlist[i]
    stoutlist[0] = U; // original field
    for (int i = 1; i < stoutlist.size(); ++i) {
        stout_smear1(stoutlist[i - 1], stoutlist[i], coeff);
    }
}

template <typename T, int N = T::rows(), int NSTP = 2 * (NDIM - 1),
          typename atype = hila::arithmetic_type<T>>
void stout_smear_force(const std::vector<GaugeField<T>> &stoutlist,
                       const VectorField<Algebra<T>> &K, out_only VectorField<Algebra<T>> &KS,
                       atype coeff) {
    // uses the list stoutlist[] of smeared gauge fields to compute the pullback of the
    // algebra-valued force field K under the smearing and returns the force acting on the unsmeared
    // link variables as algebra-valued field KS
    // Note: our definition of the force field is different from the one used in [arXiv:hep-lat/0311018v1],
    // in order to match the force field representation used by our HMC implementation.
    VectorField<T> stapl[NSTP];
    //Field<T> staps;
    VectorField<T> staps;
    VectorField<T> K1;

    foralldir(d1) onsites(ALL) KS[d1][X] = K[d1][X];

    for (int i = stoutlist.size() - 2; i >= 0; --i) {
        const GaugeField<T> &U0 = stoutlist[i + 1];
        const GaugeField<T> &U = stoutlist[i];

        get_stout_staples(U, stapl, staps);
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
                // mdtexp[i][j] = trace(Q.dagger() * KS[d1][X].expand() * dQ/dX[j][i]) :
                mult_chexp(tplaqs.project_to_algebra_scaled(-coeff).expand(), KS[d1][X].expand(),
                           mtexp, mdtexp);
                
                // set K1[d1][X] to be the equivalent of the \Lambda matrix from eq.(73) in
                // [arXiv:hep-lat/0311018v1]:
                K1[d1][X] = mdtexp.project_to_algebra_scaled(coeff).expand();

                // equivalent of first line and first term on second line of eq.(75) in
                // [arXiv:hep-lat/0311018v1]:
                KS[d1][X] = (mtexp - tplaqs * K1[d1][X]).project_to_algebra();

                // multiply K1[d1] by U[d1]:
                K1[d1][X] *= U[d1][X];
                // (the latter is done so that the parallel transport of
                // U[d1][X] * stapl[][d1][X] * K1[d1][X] in negative d1-direction, required
                // in the next section, can be computed simply as stapl[][d1][X-d1] * K1[d1][X-d1]
            }

            // multiply all the stapl[][d1] by K1[d1] and start sending the resulting matrices
            // to the neighboring sites in negative d1-direction: 
            // (this is done in the same order in which they will be used below)
            int istp = 0;
            int istn;
            foralldir(d2) if (d2 != d1) {
                onsites(ALL) stapl[istp][d1][X] *= K1[d1][X];
                stapl[istp][d1].start_gather(-d1, ALL);
                istn = NSTP - 1 - istp;
                onsites(ALL) stapl[istn][d1][X] *= K1[d1][X];
                stapl[istn][d1].start_gather(-d1, ALL);
                ++istp;
            }
        }

        foralldir(d1) {
            // compute the equivalent of the terms in the curly bracket of eq.(75) in
            // [arXiv:hep-lat/0311018v1]:
            int istp = 0;
            int istn;
            foralldir(d2) if (d2 != d1) {
                // term from d1-d2-plaquette for positive d2-direction:
                onsites(ALL) {
                    //staps[X] = stapl[istp][d1][X - d1];
                    staps[d1][X] = stapl[istp][d1][X - d1];
                }
                // define for each modified stout staple matrix the correspondin path of gauge
                // links:
                std::vector<Direction> pathp = {d2, -d1, -d2};
                //get_wloop_force_from_wl_add(U, pathp, staps, 1.0, KS);
                get_wloop_force_from_wl_add(U, pathp, staps[d1], 1.0, KS);

                // term from d1-d2-plaquette for negative d2-direction:
                istn = NSTP - 1 - istp;
                onsites(ALL) {
                    //staps[X] = stapl[istn][d1][X - d1];
                    staps[d1][X] = stapl[istn][d1][X - d1];
                }
                // define for each modified stout staple matrix the correspondin path of gauge
                // links:
                std::vector<Direction> pathn = {-d2, -d1, d2};
                //get_wloop_force_from_wl_add(U, pathn, staps, 1.0, KS);
                get_wloop_force_from_wl_add(U, pathn, staps[d1], 1.0, KS);
                ++istp;
            }
        }
    }
}

template <typename T, typename atype = hila::arithmetic_type<T>>
void stout_smear1k(const GaugeField<T> &U, out_only GaugeField<T> &stout,
                  out_only GaugeField<T> &stoutk, atype coeff) {
    staplesums(U, stout);
    foralldir(d) {
        onsites(ALL) {
            stout[d][X] = chexpk((U[d][X] * stout[d][X]).project_to_algebra_scaled(-coeff).expand(),
                                 stoutk[d][X]) * U[d][X];
        }
    }
}

template <typename T, typename atype = hila::arithmetic_type<T>>
void stout_smeark(const GaugeField<T> &U, out_only std::vector<GaugeField<T>> &stoutlist,
                 out_only std::vector<GaugeField<T>> &stoutklist, atype coeff) {
    // performs nst stout smearing steps on the gauge field U and stores the gague field obtained
    // after the i-th smearing step in stoutlist[i]
    stoutlist[0] = U; // original field
    for (int i = 1; i < stoutlist.size(); ++i) {
        stout_smear1k(stoutlist[i - 1], stoutlist[i], stoutklist[i - 1], coeff);
    }
}

template <typename T, int N = T::rows(), int NSTP = 2 * (NDIM - 1),
          typename atype = hila::arithmetic_type<T>>
void stout_smeark_force(const std::vector<GaugeField<T>> &stoutlist,
                       const std::vector<GaugeField<T>> &stoutklist,
                       const VectorField<Algebra<T>> &K, out_only VectorField<Algebra<T>> &KS,
                       atype coeff) {
    // uses the list stoutlist[] of smeared gauge fields to compute the pullback of the
    // algebra-valued force field K under the smearing and returns the force acting on the unsmeared
    // link variables as algebra-valued field KS
    // Note: our definition of the force field is different from the one used in
    // [arXiv:hep-lat/0311018v1], in order to match the force field representation used by our HMC
    // implementation.
    VectorField<T> stapl[NSTP];
    // Field<T> staps;
    VectorField<T> staps;
    VectorField<T> K1;

    foralldir(d1) onsites(ALL) KS[d1][X] = K[d1][X];

    for (int i = stoutlist.size() - 2; i >= 0; --i) {
        const GaugeField<T> &U0 = stoutlist[i + 1];
        const GaugeField<T> &U = stoutlist[i];
        const GaugeField<T> &dUK = stoutklist[i];

        get_stout_staples(U, stapl, staps);
        foralldir(d1) {
            // get_stout_staples(U, d1, stapl, staps);
            onsites(ALL) {
                // compute stout smearing operator and its derivatives:

                // temp. variables:
                T mtexp;
                T mdtexp;
                // turn staple sum into plaquette sum by multiplying with link variable:
                // T tplaqs = U[d1][X] * staps[X];
                T tplaqs = U[d1][X] * staps[d1][X];
                // the following function computes first for X = -coeff * tplaqs the smearing
                // operator Q = exp(X) and its derivatives dQ/dX[][], and uses these to
                // compute the two matrices:
                // mtexp = Q.dagger() * KS[d1][X].expand() * Q
                // and
                // mdtexp[i][j] = trace(Q.dagger() * KS[d1][X].expand() * dQ/dX[j][i]) :
                mult_chexpk_fast(tplaqs.project_to_algebra_scaled(-coeff).expand(),
                                 U0[d1][X] * U[d1][X].dagger(), dUK[d1][X], KS[d1][X].expand(),
                                 mtexp, mdtexp);

                // set K1[d1][X] to be the equivalent of the \Lambda matrix from eq.(73) in
                // [arXiv:hep-lat/0311018v1]:
                K1[d1][X] = mdtexp.project_to_algebra_scaled(coeff).expand();

                // equivalent of first line and first term on second line of eq.(75) in
                // [arXiv:hep-lat/0311018v1]:
                KS[d1][X] = (mtexp - tplaqs * K1[d1][X]).project_to_algebra();

                // multiply K1[d1] by U[d1]:
                K1[d1][X] *= U[d1][X];
                // (the latter is done so that the parallel transport of
                // U[d1][X] * stapl[][d1][X] * K1[d1][X] in negative d1-direction, required
                // in the next section, can be computed simply as stapl[][d1][X-d1] * K1[d1][X-d1]
            }

            // multiply all the stapl[][d1] by K1[d1] and start sending the resulting matrices
            // to the neighboring sites in negative d1-direction:
            // (this is done in the same order in which they will be used below)
            int istp = 0;
            int istn;
            foralldir(d2) if (d2 != d1) {
                onsites(ALL) stapl[istp][d1][X] *= K1[d1][X];
                stapl[istp][d1].start_gather(-d1, ALL);
                istn = NSTP - 1 - istp;
                onsites(ALL) stapl[istn][d1][X] *= K1[d1][X];
                stapl[istn][d1].start_gather(-d1, ALL);
                ++istp;
            }
        }

        foralldir(d1) {
            // compute the equivalent of the terms in the curly bracket of eq.(75) in
            // [arXiv:hep-lat/0311018v1]:
            int istp = 0;
            int istn;
            foralldir(d2) if (d2 != d1) {
                // term from d1-d2-plaquette for positive d2-direction:
                onsites(ALL) {
                    // staps[X] = stapl[istp][d1][X - d1];
                    staps[d1][X] = stapl[istp][d1][X - d1];
                }
                // define for each modified stout staple matrix the correspondin path of gauge
                // links:
                std::vector<Direction> pathp = {d2, -d1, -d2};
                // get_wloop_force_from_wl_add(U, pathp, staps, 1.0, KS);
                get_wloop_force_from_wl_add(U, pathp, staps[d1], 1.0, KS);

                // term from d1-d2-plaquette for negative d2-direction:
                istn = NSTP - 1 - istp;
                onsites(ALL) {
                    // staps[X] = stapl[istn][d1][X - d1];
                    staps[d1][X] = stapl[istn][d1][X - d1];
                }
                // define for each modified stout staple matrix the correspondin path of gauge
                // links:
                std::vector<Direction> pathn = {-d2, -d1, d2};
                // get_wloop_force_from_wl_add(U, pathn, staps, 1.0, KS);
                get_wloop_force_from_wl_add(U, pathn, staps[d1], 1.0, KS);
                ++istp;
            }
        }
    }
}

#endif