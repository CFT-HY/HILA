#ifndef STOUT_SMEAR_H_
#define STOUT_SMEAR_H_

#include "hila.h"

#include "gauge/wilson_line_and_force.h"

/////////////////////////////////////////////////////////////////////////////
/// Do stout smearing for the gauge fields
///  U' = exp( -coeff * Proj_to_Alg(U * Staples) ) * U

template <typename T>
void sstaplesum(const GaugeField<T> &U, out_only Field<T> &staples, Direction d1) {

    Field<T> lower;
    bool first = true;
    foralldir(d2) if (d2 != d1) {

        // anticipate that these are needed
        // not really necessary, but may be faster
        U[d2].start_gather(d1, ALL);
        U[d1].start_gather(d2, ALL);

        // calculate first lower 'u' of the staple sum
        // do it on opp parity
        onsites(ALL) {
            lower[X] = (U[d1][X] * U[d2][X + d1]).dagger() * U[d2][X];
        }
        lower.start_gather(-d2, ALL);

        // calculate then the upper 'n', and add the lower
        if (first) {
            onsites(ALL) {
                staples[X] = U[d2][X + d1] * (U[d2][X] * U[d1][X + d2]).dagger();
            }
            first = false;
        } else {
            onsites(ALL) {
                staples[X] += U[d2][X + d1] * (U[d2][X] * U[d1][X + d2]).dagger();
            }
        }

        // add lower
        onsites(ALL) {
            staples[X] += lower[X - d2];
        }
    }
}


template <typename T, typename atype = hila::arithmetic_type<T>>
void stout_smear1(const GaugeField<T> &U, Field<T> &s, atype coeff, Direction d) {
    sstaplesum(U, s, d);
    onsites(ALL) {
        s[X] = chexp((U[d][X] * s[X]).project_to_algebra_scaled(-coeff).expand()) * U[d][X];
    }
}

template <typename T, typename atype = hila::arithmetic_type<T>>
void stout_smear1(const GaugeField<T> &U, GaugeField<T> &stout, atype coeff) {
    foralldir(d) stout_smear1(U,stout[d],coeff,d);
}


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

// functions used to do stout smearing and perform the computation of the pullback of the force for
// smeared actions:

template <typename group, const int nstap = 2 * (NDIM - 1)>
void get_stout_staples(const GaugeField<group> &U, Direction d1,
                          VectorField<group>(out_only &S)[nstap], out_only Field<group> &Ssum) {
    // compute the staples spanned by directions (d1,d2), where d2 runs through all positive and
    // negative directions which are orthogonal to d1. Individual staple fields get stored in
    // S[d1][d2'] and their sum over d2' in Ssum.

    foralldir(d2) if (d2 != d1) {
        U[d2].start_gather(d1, ALL);
        U[d1].start_gather(d2, ALL);
    }
    onsites(ALL) {
        Ssum[X] = 0;
    }
    Field<group> tF;

    int k = 0;
    int rk;
    foralldir(d2) if (d2 != d1) {
        // shifted lower staple
        onsites(ALL) {
            tF[X] = (U[d1][X] * U[d2][X + d1]).dagger() * U[d2][X];
        }
        tF.start_gather(-d2, ALL);

        // upper staple
        onsites(ALL) {
            S[k][d1][X] = U[d2][X + d1] * (U[d2][X] * U[d1][X + d2]).dagger();
            Ssum[X] += S[k][d1][X];
        }

        rk = nstap - 1 - k;
        // lower staple
        onsites(ALL) {
            S[rk][d1][X] = tF[X - d2];
            Ssum[X] += S[rk][d1][X];
        }
        ++k;
    }
}

template <typename T, typename atype = hila::arithmetic_type<T>>
void stout_smear(const GaugeField<T> &U, out_only std::vector<GaugeField<T>> &stoutlist, atype coeff) {
    // performs nst stout smearing steps on the gauge field U and stores the gague field obtained
    // after the i-th smearing step in stoutlist[i]
    stoutlist[0] = U; // original field
    for (int i = 1; i < stoutlist.size(); ++i) {
        stout_smear1(stoutlist[i-1], stoutlist[i], coeff);
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
    Field<T> staps;
    VectorField<T> K1;

    foralldir(d1) onsites(ALL) KS[d1][X] = K[d1][X];

    for (int i = stoutlist.size() - 2; i >= 0; --i) {
        const GaugeField<T> &U0 = stoutlist[i + 1];
        const GaugeField<T> &U = stoutlist[i];

        foralldir(d1) {
            get_stout_staples(U, d1, stapl, staps);
            onsites(ALL) {
                // compute stout smearing operator and its derivatives:

               // temp. variables:
                T mtexp;
                T mdtexp;
                // turn staple sum into plaquette sum by multiplying with link variable:
                T tplaqs = U[d1][X] * staps[X];
                // the following function computes first for X = -coeff * tplaqs the smearing
                // operator Q = exp(X) and its derivatives dQ/dX[][], and uses these to
                // compute the two matrices: 
                // mtexp = Q.dagger() * KS[d1][X].expand() * Q 
                // and 
                // mdtexp[i][j] = trace(Q.dagger * KS[d1][X].expand() * dQ/dX[j][i]) :
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
                    staps[X] = stapl[istp][d1][X - d1];
                }
                // define for each modified stout staple matrix the correspondin path of gauge
                // links:
                std::vector<Direction> pathp = {d2, -d1, -d2};
                get_wloop_force_from_wl_add(U, pathp, staps, 1.0, KS);

                // term from d1-d2-plaquette for negative d2-direction:
                istn = NSTP - 1 - istp;
                onsites(ALL) {
                    staps[X] = stapl[istn][d1][X - d1];
                }
                // define for each modified stout staple matrix the correspondin path of gauge
                // links:
                std::vector<Direction> pathn = {-d2, -d1, d2};
                get_wloop_force_from_wl_add(U, pathn, staps, 1.0, KS);
                ++istp;
            }
        }
    }
}

#endif