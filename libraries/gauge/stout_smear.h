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

        // calculate first lower 'U' of the staple sum
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

template <typename group, const int nstap=2*(NDIM-1)>
void get_stout_plaquettes(const GaugeField<group> &U, Direction d1,
                       VectorField<group>(out_only &S)[nstap], out_only Field<group> &Ssum) {
    // compute the plaquettes spanned by directions (d1,d2), where d2 runs through all positive and
    // negative directions which are orthogonal to d1. Individual plauqette fields get stored in
    // S[d1][d2'] and their sum over d2' in Ssum.

    foralldir(d2) if (d2 != d1) {
        U[d2].start_gather(d1, ALL);
        U[d1].start_gather(d2, ALL);
    }
    onsites(ALL){
        Ssum[X] = 0;
    }
    Field<group> tF;

    int k = 0;
    int rk;
    foralldir(d2) if (d2 != d1) {
        onsites(ALL) {
            tF[X] = U[d1][X + d2] * (U[d1][X] * U[d2][X + d1]).dagger();
            // d1-d2-plaquette with positive d2 that starts and ends at X
            S[k][d1][X] = (U[d2][X] * tF[X]).dagger();
            tF[X] = tF[X] * U[d2][X];
        }
        tF.start_gather(-d2, ALL);
        onsites(ALL) {
            Ssum[X] += S[k][d1][X];
        }
        rk = nstap - 1 - k;
        onsites(ALL) {
            // d1-d2-plaquette with negative d2 that starts and ends at X
            S[rk][d1][X] = tF[X - d2];
            Ssum[X] += S[rk][d1][X];
        }
        ++k;
    }
}

template <typename T, int nst, typename atype = hila::arithmetic_type<T>>
void stout_smear(const GaugeField<T> &U, GaugeField<T>(out_only &stoutlist)[nst],
                 atype coeff) {
    // performs nst stout smearing steps on the gauge field U and stores the gague field obtained
    // after the i-th smearing step in stoutlist[i]
    stoutlist[0] = U; // original field
    for (int i = 1; i < nst; ++i) {
        stout_smear1(stoutlist[i-1], stoutlist[i], coeff);
    }
}

template <typename T, int nst, int N = T::rows(), int NSTP = 2 * (NDIM - 1),
          typename atype = hila::arithmetic_type<T>>
void stout_smear_force(const GaugeField<T> (&stoutlist)[nst], const VectorField<Algebra<T>> &K,
                       out_only VectorField<Algebra<T>> &KS, atype coeff) {
    // uses the list stoutlist[] of smeared gauge fields to compute the pullback of the
    // algebra-valued force field K under the smearing and returs the force acting on the unsmeared
    // link variables as algebra-valued field KS:
    VectorField<T> plaql[NSTP];
    Field<T> plaqs;
    VectorField<T> K1;

    foralldir(d1) onsites(ALL) KS[d1][X] = K[d1][X];

    for (int i = nst - 2; i >= 0; --i) {
        const GaugeField<T> &U0 = stoutlist[i + 1];
        const GaugeField<T> &U = stoutlist[i];

        foralldir(d1) {
            get_stout_plaquettes(U,d1,plaql,plaqs);
            onsites(ALL) {
                // compute stout smearing operator and its derivatives:
                /*
                T texp;
                T dtexp[N][N];
                chexp(plaqs[X].project_to_algebra_scaled(-coeff).expand(), texp, dtexp);

                T m0 = texp.dagger() * KS[d1][X].expand();
                KS[d1][X] = (m0 * texp).project_to_algebra();
                T m1; // equivalent of Gamma matrix, eq. (74) in [arXiv:hep-lat/0311018v1]:
                for (int ic1 = 0; ic1 < N; ++ic1) {
                    for (int ic2 = 0; ic2 < N; ++ic2) {
                        m1.e(ic1, ic2) = mul_trace(m0, dtexp[ic2][ic1]);
                    }
                }
                // equivalent of Lambda matrix, eq. (73) in [arXiv:hep-lat/0311018v1], 
                // multiplied by coeff:
                K1[d1][X] = m1.project_to_algebra_scaled(coeff).expand();
                */
                T mtexp;
                T mdtexp;
                mult_chexp(plaqs[X].project_to_algebra_scaled(-coeff).expand(), KS[d1][X].expand(),
                           mtexp, mdtexp);
                KS[d1][X] = mtexp.project_to_algebra();
                K1[d1][X] = mdtexp.project_to_algebra_scaled(coeff).expand();
            }
        }
        foralldir(d1) {
            int istp = 0;
            foralldir(d2) if (d2 != d1) {
                // d1-d2-plaquette for positive d2-direction:
                 onsites(ALL) {
                    plaqs[X] = plaql[istp][d1][X] * K1[d1][X];
                }
                // define for each stout plaquette matrix the correspondin path of gauge links
                std::vector<Direction> pathp = {d1, d2, -d1, -d2};
                get_wloop_force_from_wl_add(U, pathp, plaqs, 1.0, KS);

                // d1-d2-plaquette for negative d2-direction:
                int istn = NSTP - 1 - istp;
                onsites(ALL) {
                    plaqs[X] = plaql[istn][d1][X] * K1[d1][X];
                }
                // define for each stout plaquette matrix the correspondin path of gauge links
                std::vector<Direction> pathn = {d1, -d2, -d1, d2};
                get_wloop_force_from_wl_add(U, pathn, plaqs, 1.0, KS);
                ++istp;
            }
        }

    }
}

#endif