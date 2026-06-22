/** @file gradient_flow.h */

#ifndef GRADIENT_FLOW_H_
#define GRADIENT_FLOW_H_

#include "hila.h"
#include "gauge/energy_and_topo_charge_clover.h"
#include "gauge/energy_and_topo_charge_log.h"
#include "gauge/wilson_plaquette_action.h"
#include "gauge/clover_action.h"
#include "gauge/log_clover_action.h"

#include "gauge/bulk_prevention_action.h"
#include "gauge/improved_action.h"
#include "gauge/log_plaquette_action.h"

#include "tools/string_format.h"

enum GradientFlowAction {
    Wilson = 0,
    Bulk_Prevention = 1,
    Luscher_Weisz = 2,
    Iwasaki = 3,
    DBW2 = 4,
    Log_Plaq = 5,
    Zeuthen = 6,
    N_GradientFlowAction,
};

// Make it possible to select the gradient flow action either at compile time or at run time.
//
// if GFLOWACTION is defined action is selected at compile time.
//
// Otherwise the user code must declare and set the global variable `g_GF_action` somewhere at run
// time to select the used action.
#ifdef GFLOWACTION
constexpr GradientFlowAction g_GF_action = static_cast<GradientFlowAction>(GFLOWACTION);
static_assert((0 < GFLOWACTION) && (GFLOWACTION < N_GradientFlowAction),
              "Unknown GradientFlowAction");
#else
extern GradientFlowAction g_GF_action;
#endif

/// print info about used flow action to stdout from rank 0.
static inline void print0_gf_info() {
    switch (g_GF_action) {
    case GradientFlowAction::Wilson: // = 0
        hila::out0 << "GFINFO using Wilson's plaquette action\n";
        break;
    case GradientFlowAction::Bulk_Prevention: // = 1
        hila::out0 << "GFINFO using bulk-prevention action\n";
        break;
    case GradientFlowAction::Luscher_Weisz: // = 2
        hila::out0 << "GFINFO using Luscher-Weisz action\n";
        break;
    case GradientFlowAction::Iwasaki: // = 3
        hila::out0 << "GFINFO using Iwasaki action\n";
        break;
    case GradientFlowAction::DBW2: // = 4
        hila::out0 << "GFINFO using DBW2 action\n";
        break;
    case GradientFlowAction::Log_Plaq: // = 5
        hila::out0 << "GFINFO using log-plaquette action\n";
        break;
    case GradientFlowAction::Zeuthen: // = 6
        hila::out0 << "GFINFO using Zeuthen flow\n";
        break;
    default:
        assert(false && "Unknown GradientFlowAction");
    }
}

template <typename group, typename atype = hila::arithmetic_type<group>>
void get_zeuthen_flow_force(const GaugeField<group> &U, VectorField<Algebra<group>> &E,
                            atype eps = 1.0) {

    // first get the Luscher-Weisz improved force
    VectorField<Algebra<group>> F_LW;
    atype c12 = -1.0 / 12.0; // rectangle weight
    atype c11 = 5.0 / 3.0;   // plaquette weight
    get_force_impr(U, F_LW, eps * c11, eps * c12);

    // Zeuthen flow force = (1 + 1/12 D_\mu D_\mu^* ) F_LW

    foralldir (d) {
        onsites (ALL) {
            E[d][X] =
                F_LW[d][X] +
                (1.0 / 12.0) * ( // calculate gradient term improvement
                                 // (\mu=d, no sum) D_\mu D_\mu^* F_LW(X)
                                   left_conjugation(U[d][X], F_LW[d][X])            // U F U^\dagger
                                   - 2.0 * F_LW[d][X]                               //
                                   + right_conjugation(U[d][X - d], F_LW[d][X - d]) // U^dagger F U
                               );
        }
    }
}

template <typename group, typename atype = hila::arithmetic_type<group>>
void get_gf_force(const GaugeField<group> &U, VectorField<Algebra<group>> &E) {
    // wrapper for force computation routine to be used for gradient flow

    atype eps = 1.0; // in principle need factor 2.0 here to switch from unoriented to oriented
                     // plaquettes (factor is usually absorbed in \beta, but gradient flow force
                     // is computed from action term with \beta-factor stripped off)
                     // however; it seems that in practice factor 1.0 is used.
                     // note: when switching to factor 2.0, remember to change also the stability
                     // limit in the do_gradient_flow_adapt() below

    switch (g_GF_action) {
    case GradientFlowAction::Wilson: { // = 0
        get_force_wplaq(U, E, eps);
    } break;
    case GradientFlowAction::Bulk_Prevention: { // = 1
        get_force_bp(U, E, eps);
    } break;
    case GradientFlowAction::Luscher_Weisz: { // = 2
        atype c12 = -1.0 / 12.0;              // rectangle weight
        atype c11 = 1.0 - 8.0 * c12;          // plaquette weight
        get_force_impr(U, E, eps * c11, eps * c12);
    } break;
    case GradientFlowAction::Iwasaki: { // = 3
        atype c12 = -0.331;             // rectangle weight
        atype c11 = 1.0 - 8.0 * c12;    // plaquette weight
        get_force_impr(U, E, eps * c11, eps * c12);
    } break;
    case GradientFlowAction::DBW2: { // = 4
        atype c12 = -1.4088;         // rectangle weight
        atype c11 = 1.0 - 8.0 * c12; // plaquette weight
        get_force_impr(U, E, eps * c11, eps * c12);
    } break;
    case GradientFlowAction::Log_Plaq: { // = 5
        get_force_log_plaq(U, E, eps);
    } break;
    case GradientFlowAction::Zeuthen: { // = 6
        get_zeuthen_flow_force(U, E, eps);
    } break;
    default:
        assert(false && "Unknown GradientFlowAction");
    }
    // 2/9*BP+7/9*Wilson
    // get_force_bp(U,E,eps*2.0/9.0);
    // get_force_wplaq_add(U,E,eps*7.0/9.0);
};

template <typename group, typename atype = hila::arithmetic_type<group>>
atype measure_gf_s(const GaugeField<group> &U) {
    // wrapper for gauge action computation routine
    // (for gauge action that is used to evolve the flow)
    atype res;
    switch (g_GF_action) {
    case GradientFlowAction::Wilson: { // = 0
        res = measure_s_wplaq(U);
    } break;
    case GradientFlowAction::Bulk_Prevention: { // = 1
        res = measure_s_bp(U);
    } break;
    case GradientFlowAction::Zeuthen:         // = 6
    case GradientFlowAction::Luscher_Weisz: { // = 2
        atype c12 = -1.0 / 12.0;              // rectangle weight
        atype c11 = 1.0 - 8.0 * c12;          // plaquette weight
        res = measure_s_impr(U, c11, c12);
    } break;
    case GradientFlowAction::Iwasaki: { // = 3
        atype c12 = -0.331;             // rectangle weight
        atype c11 = 1.0 - 8.0 * c12;    // plaquette weight
        res = measure_s_impr(U, c11, c12);
    } break;
    case GradientFlowAction::DBW2: { // = 4
        atype c12 = -1.4088;         // rectangle weight
        atype c11 = 1.0 - 8.0 * c12; // plaquette weight
        res = measure_s_impr(U, c11, c12);
    } break;
    case GradientFlowAction::Log_Plaq: { // = 5
        res = measure_s_log_plaq(U);
    } break;
    default:
        assert(false && "Unknown GradientFlowAction");
    }
    // 2/9*BP+7/9*Wilson
    // atype res=measure_s_bp(U,E,eps*2.0/9.0)+measure_s_wplaq(U,E,eps*7.0/9.0);

    return res;
}

template <typename group, typename atype = hila::arithmetic_type<group>>
atype measure_dE_wplaq_dt(const GaugeField<group> &U) {
    Reduction<double> de = 0;
    de.allreduce(false).delayed(true);
    VectorField<Algebra<group>> K, Kc;
    get_gf_force(U, K);
    get_force_wplaq(U, Kc, -2.0);
    foralldir (d) {
        onsites (ALL) {
            de += Kc[d][X].dot(K[d][X]);
        }
    }
    return (atype)de.value();
}

template <typename group, typename atype = hila::arithmetic_type<group>>
atype measure_dE_clov_dt(const GaugeField<group> &U) {
    Reduction<double> de = 0;
    de.allreduce(false).delayed(true);
    VectorField<Algebra<group>> K, Kc;
    get_gf_force(U, K);
    get_force_clover(U, Kc, -2.0);
    foralldir (d) {
        onsites (ALL) {
            de += Kc[d][X].dot(K[d][X]);
        }
    }
    return (atype)de.value();
}

template <typename group, typename atype = hila::arithmetic_type<group>>
atype measure_dE_log_dt(const GaugeField<group> &U) {
    Reduction<double> de = 0;
    de.allreduce(false).delayed(true);
    VectorField<Algebra<group>> K, Kc;
    get_gf_force(U, K);
    get_force_log(U, Kc, -2.0);
    foralldir (d) {
        onsites (ALL) {
            de += Kc[d][X].dot(K[d][X]);
        }
    }
    return (atype)de.value();
}

template <typename group, typename atype = hila::arithmetic_type<group>>
void measure_gradient_flow_stuff(const GaugeField<group> &V, atype flow_l, atype t_step) {
    // perform measurements on flowed gauge configuration V at flow scale flow_l
    // [t_step is the flow time integration step size used in last gradient flow step]
    // and print results in formatted form to standard output
    static bool first = true;
    if (first) {
        print0_gf_info();
        // print legend for flow measurement output
        hila::out0 << "LGFLMEAS  l(ambda)        S-flow        S-plaq        E_plaq    dE_plaq/dl  "
                      "       E_clv     dE_clv/dl     Qtopo_clv         E_log     dE_log/dl     "
                      "Qtopo_log   [t_step_size]   [max_S-plaq]\n";
        first = false;
    }
    atype slocal = measure_gf_s(V) /
                   (lattice.volume() * NDIM * (NDIM - 1) / 2); // average action per plaquette

    atype max_plaq = 0;
    atype plaq = measure_s_wplaq(V, max_plaq) /
                 (lattice.volume() * NDIM * (NDIM - 1) / 2); // average wilson plaquette action
    atype eplaq = plaq * NDIM * (NDIM - 1) *
                  group::size(); // naive energy density (based on wilson plaquette action)

    // average energy density and toplogical charge from
    // clover definition of field strength tensor :
    atype qtopocl, ecl;
    measure_topo_charge_and_energy_clover(V, qtopocl, ecl);
    ecl /= lattice.volume();

    // average energy density and toplogical charge from
    // symmetric log definition of field strength tensor :
    atype qtopolog, elog;
    measure_topo_charge_and_energy_log(V, qtopolog, elog);
    elog /= lattice.volume();

    // derivative of plaquette energy density w.r.t. to flow time :
    atype deplaqdt = measure_dE_wplaq_dt(V) / lattice.volume();

    // derivative of clover energy density w.r.t. to flow time :
    atype declovdt = measure_dE_clov_dt(V) / lattice.volume();

    // derivative of log energy density w.r.t. to flow time :
    atype delogdt = measure_dE_log_dt(V) / lattice.volume();

    // print formatted results to standard output :
    hila::out0 << string_format("GFLMEAS  % 9.3f % 0.6e % 0.6e % 0.6e % 0.6e % 0.6e % 0.6e % 0.6e "
                                "% 0.6e % 0.6e % 0.6e     [%0.3e]    [%0.3e]",
                                flow_l, slocal, plaq, eplaq, 0.25 * flow_l * deplaqdt, ecl,
                                0.25 * flow_l * declovdt, qtopocl, elog, 0.25 * flow_l * delogdt,
                                qtopolog, t_step, max_plaq)
               << '\n';
}

template <typename group, typename atype = hila::arithmetic_type<group>>
atype do_gradient_flow_adapt(GaugeField<group> &V, atype l_start, atype l_end, atype atol = 1.0e-6,
                             atype rtol = 1.0e-4, atype tstep = 0.0) {
    // wilson flow integration from flow scale l_start to l_end using 3rd order
    // 3-step Runge-Kutta (RK3) from arXiv:1006.4518 (cf. appendix C of
    // arXiv:2101.05320 for derivation of this Runge-Kutta method)
    // and embedded RK2 for adaptive step size

    atype esp = 3.0;        // expected single step error scaling power: err ~ step^(esp)
                            //   - for RK3 with embedded RK2: esp \approx 3.0
    atype iesp = 1.0 / esp; // inverse single step error scaling power

    atype stepmf = 1.0;
    atype maxstepmf = 10.0; // max. growth factor of adaptive step size
    atype minstepmf = 0.1;  // min. growth factor of adaptive step size

    // translate flow scale interval [l_start,l_end] to corresponding
    // flow time interval [t,tmax] :
    atype t = l_start * l_start / 8.0;
    atype tmax = l_end * l_end / 8.0;

    atype ubstep = (tmax - t) / 2.0; // max. allowed time step

    atype tatol = atol * sqrt(2.0);

    // hila::out0<<"t: "<<t<<" , tmax: "<<tmax<<" , step: "<<tstep<<" , minmaxreldiff:
    // "<<minmaxreldiff<<"\n";

    // temporary variables :
    VectorField<Algebra<group>> k1, k2, tk;
    GaugeField<group> V2, V0;
    Field<atype> reldiff;

    // RK3 coefficients from arXiv:1006.4518 :
    // correspond to standard RK3 with Butcher-tableau
    // (cf. arXiv:2101.05320, Appendix C)
    //  0  |   0     0     0
    //  #  |  1/4    0     0
    //  #  | -2/9   8/9    0
    // -------------------------
    //     |  1/4    0    3/4
    //
    atype a11 = 0.25;
    atype a21 = -17.0 / 36.0, a22 = 8.0 / 9.0;
    atype a33 = 0.75;

    // RK2 coefficients :
    // cf. Alg. 6 and Eqn. (13)-(14) in arXiv:2101.05320 to see
    // how these are obtained from standard RK2 with Butcher-tableau
    //  0  |   0     0
    //  #  |  1/4    0
    // -----------------
    //     |  -1     2
    //
    atype b21 = -1.25, b22 = 2.0;

    atype step = min(tstep, ubstep); // initial step size

    if (t == 0 || step == 0) {
        // when using a gauge action for gradient flow that is different from
        // the one used to sample the gauge cofingurations, the initial force
        // can be huge. Therefore, if no t_step is provided as input, the inital
        // value for step is here adjustet so that
        // step * <largest local force> = maxstk
        atype maxstk = 1.0e-1;

        // get max. local gauge force:
        get_gf_force(V, tk);
        atype maxtk = 0.0;
        foralldir (d) {
            onsites (ALL) {
                reldiff[X] = (tk[d][X].squarenorm());
            }
            atype tmaxtk = reldiff.max();
            if (tmaxtk > maxtk) {
                maxtk = tmaxtk;
            }
        }
        maxtk = sqrt(0.5 * maxtk);

        if (step == 0) {
            if (maxtk > maxstk) {
                step = min(maxstk / maxtk,
                           ubstep); // adjust initial step size based on max. force magnitude
                hila::out0 << "GFINFO using max. gauge force (max_X |F(X)|=" << maxtk
                           << ") to set initial flow time step size: " << step << "\n";
            } else {
                step = min((atype)1.0, ubstep);
            }
        } else if (step * maxtk > maxstk) {
            step = min(maxstk / maxtk,
                       ubstep); // adjust initial step size based on max. force magnitude
            hila::out0 << "GFINFO using max. gauge force (max_X |F(X)|=" << maxtk
                       << ") to set initial flow time step size: " << step << "\n";
        }
    }


    V0 = V;
    bool stop = false;
    while (t < tmax && !stop) {
        tstep = step;
        if (t + step >= tmax) {
            step = tmax - t;
            stop = true;
        }

        get_gf_force(V, k1);
        foralldir (d)
            onsites (ALL) {
                // first steps of RK3 and RK2 are the same :
                V[d][X] = chexp(k1[d][X] * (step * a11)) * V[d][X];
            }

        get_gf_force(V, k2);
        foralldir (d)
            onsites (ALL) {
                // second step of RK2 :
                // (tk[d][X] will be used for rel. error computation)
                tk[d][X] = k2[d][X];
                tk[d][X] *= (step * b22);
                tk[d][X] += k1[d][X] * (step * b21);
                V2[d][X] = chexp(tk[d][X]) * V[d][X];

                // second step of RK3 :
                k2[d][X] *= (step * a22);
                k2[d][X] += k1[d][X] * (step * a21);
                V[d][X] = chexp(k2[d][X]) * V[d][X];
            }

        get_gf_force(V, k1);
        foralldir (d)
            onsites (ALL) {
                // third step of RK3 :
                k1[d][X] *= (step * a33);
                k1[d][X] -= k2[d][X];
                V[d][X] = chexp(k1[d][X]) * V[d][X];
            }

        // determine maximum difference between RK3 and RK2,
        // relative to desired accuracy :
        atype relerr = 0.0;
        foralldir (d) {
            onsites (ALL) {
                reldiff[X] = (V2[d][X] * V[d][X].dagger()).project_to_algebra().norm() /
                             (tatol + rtol * tk[d][X].norm() / step);
                // note: we divide tk.norm() by step to have consistent leading stepsize dependency
                // no mather whether relative or absolute error tollerance dominates
            }
            atype trelerr = reldiff.max();
            if (trelerr > relerr) {
                relerr = trelerr;
            }
        }

        if (relerr < 1.0) {
            // proceed to next iteration
            t += step;
            V.reunitarize_gauge();
            V0 = V;
        } else {
            // repeat current iteration if single step error was too large
            V = V0;
            stop = false;
        }

        // determine step size to achieve desired accuracy goal :
        stepmf = pow(relerr, -iesp);
        if (stepmf <= minstepmf) {
            stepmf = minstepmf;
        } else if (stepmf >= maxstepmf) {
            stepmf = maxstepmf;
        }

        // adjust step size :
        step = min((atype)0.9 * stepmf * step, ubstep);
    }

    return tstep;
}

#endif