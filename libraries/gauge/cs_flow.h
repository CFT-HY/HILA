/** @file cs_flow.h */

#ifndef CS_FLOW_H_
#define CS_FLOW_H_

#include "hila.h"
#include "gauge/energy_and_topo_charge_clover.h"
#include "gauge/energy_and_topo_charge_log.h"
#include "gauge/wilson_plaquette_action.h"
#include "gauge/clover_action.h"
#include "gauge/log_clover_action.h"

#ifdef CSFLOWACTION
#define CSFLOWS CSFLOWACTION
#else
#define CSFLOWS 3
#endif

#if CSFLOWS == 1
#include "gauge/bulk_prevention_action.h"
#elif CSFLOWS > 1
#if CSFLOWS < 5
#include "gauge/improved_action.h"
#else
#include "gauge/log_plaquette_action.h"
#endif
#endif

#include "tools/string_format.h"


template <typename group, typename atype = hila::arithmetic_type<group>>
void get_cs_force(const GaugeField<group> &U, VectorField<Algebra<group>> &K) {
    // wrapper for force computation routine to be used for constant action flow

    atype eps = 1.0; // in principle need factor 2.0 here to switch from unoriented to oriented
                     // plaquettes (factor is usually absorbed in \beta, but gradient flow force
                     // is computed from action term with \beta-factor stripped off)
                     // however; it seems that in practice factor 1.0 is used.
                     // note: when switching to factor 2.0, remember to change also the stability
                     // limit in the do_gradient_flow_adapt() below

    VectorField<Algebra<group>> E;

    atype isqrt2 = 1.0/sqrt(2.0);

#if CSFLOWS == 1 // Bulk-Prevention (BP)
    get_force_bp(U, E, eps);
#elif CSFLOWS == 2 // Luscher-Weisz (LW)
        atype c12 = -1.0 / 12.0; // rectangle weight
    atype c11 = 1.0 - 8.0 * c12; // plaquette weight
    get_force_impr(U, E, eps * c11, eps * c12);
#elif CSFLOWS == 3 // IWASAKI
        atype c12 = -0.331;      // rectangle weight
    atype c11 = 1.0 - 8.0 * c12; // plaquette weight
    get_force_impr(U, E, eps * c11, eps * c12);
#elif CSFLOWS == 4 // DBW2
        atype c12 = -1.4088;     // rectangle weight
    atype c11 = 1.0 - 8.0 * c12; // plaquette weight
    get_force_impr(U, E, eps * c11, eps * c12);
#elif CSFLOWS == 5 // LOG-PLAQUETTE
        get_force_log_plaq(U, E, eps);
#else             // WILSON
        get_force_wplaq(U, E, eps);
#endif

    foralldir(d) {
        onsites(ALL) {
            E[d][X] /= isqrt2 * E[d][X].norm();
            K[d][X] -= E[d][X].dot(K[d][X]) * E[d][X];
            K[d][X] /= isqrt2 * K[d][X].norm();
        }
    }

    // 2/9*BP+7/9*Wilson
    // get_force_bp(U,E,eps*2.0/9.0);
    // get_force_wplaq_add(U,E,eps*7.0/9.0);
};

template <typename group, typename atype = hila::arithmetic_type<group>>
atype measure_cs_s(const GaugeField<group> &U) {
    // wrapper for gauge action computation routine
    // (for gauge action that is used to evolve the flow)
#if CSFLOWS == 1 // BP
    atype res = measure_s_bp(U);
#elif CSFLOWS == 2 // LW
    atype c12 = -1.0 / 12.0;     // rectangle weight
    atype c11 = 1.0 - 8.0 * c12; // plaquette weight
    atype res = measure_s_impr(U, c11, c12);
#elif CSFLOWS == 3 // IWASAKI
    atype c12 = -0.331;          // rectangle weight
    atype c11 = 1.0 - 8.0 * c12; // plaquette weight
    atype res = measure_s_impr(U, c11, c12);
#elif CSFLOWS == 4 // DBW2
    atype c12 = -1.4088;         // rectangle weight
    atype c11 = 1.0 - 8.0 * c12; // plaquette weight
    atype res = measure_s_impr(U, c11, c12);
#elif CSFLOWS == 5 // LOG-PLAQUETTE
    atype res = measure_s_log_plaq(U);
#else             // WILSON
    atype res = measure_s_wplaq(U);
#endif

    // 2/9*BP+7/9*Wilson
    // atype res=measure_s_bp(U,E,eps*2.0/9.0)+measure_s_wplaq(U,E,eps*7.0/9.0);

    return res;
}

template <typename group, typename atype = hila::arithmetic_type<group>>
atype measure_dE_wplaq_dt(const GaugeField<group> &U, const VectorField<Algebra<group>>& E) {
    Reduction<double> de = 0;
    de.allreduce(false).delayed(true);
    VectorField<Algebra<group>> K, Kc;
    K = E;
    get_cs_force(U, K);
    get_force_wplaq(U, Kc, -2.0);
    foralldir(d) {
        onsites(ALL) {
            de += Kc[d][X].dot(K[d][X]);
        }
    }
    return (atype)de.value();
}

template <typename group, typename atype = hila::arithmetic_type<group>>
atype measure_dE_clov_dt(const GaugeField<group> &U, const VectorField<Algebra<group>> &E) {
    Reduction<double> de = 0;
    de.allreduce(false).delayed(true);
    VectorField<Algebra<group>> K, Kc;
    K = E;
    get_cs_force(U, K);
    get_force_clover(U, Kc, -2.0);
    foralldir(d) {
        onsites(ALL) {
            de += Kc[d][X].dot(K[d][X]);
        }
    }
    return (atype)de.value();
}

template <typename group, typename atype = hila::arithmetic_type<group>>
atype measure_dE_log_dt(const GaugeField<group> &U, const VectorField<Algebra<group>> &E) {
    Reduction<double> de = 0;
    de.allreduce(false).delayed(true);
    VectorField<Algebra<group>> K, Kc;
    K = E;
    get_cs_force(U, K);
    get_force_log(U, Kc, -2.0);
    foralldir(d) {
        onsites(ALL) {
            de += Kc[d][X].dot(K[d][X]);
        }
    }
    return (atype)de.value();
}

template <typename group, typename atype = hila::arithmetic_type<group>>
void measure_cs_flow_stuff(const GaugeField<group> &V, const VectorField<Algebra<group>> &K,
                           atype flow_t, atype t_step) {
    // perform measurements on flowed gauge configuration V at flow time flow_t
    // [t_step is the flow time integration step size used in last cs flow step]
    // and print results in formatted form to standard output
    static bool first = true;
    if (first) {
        // print info about used flow action
#if CSFLOWS == 1 // BP
        hila::out0 << "CSINFO using bulk-prevention action\n";
#elif CSFLOWS == 2 // LW
        hila::out0 << "CSINFO using Luscher-Weisz action\n";
#elif CSFLOWS == 3 // IWASAKI
        hila::out0 << "CSINFO using Iwasaki action\n";
#elif CSFLOWS == 4 // DBW2
        hila::out0 << "CSINFO using DBW2 action\n";
#elif CSFLOWS == 5 // LOG-PLAQUETTE
        hila::out0 << "cSINFO using log-plaquette action\n";
#else             // WILSON
        hila::out0 << "CSINFO using Wilson's plaquette action\n";
#endif
        // print legend for flow measurement output
        hila::out0 << "LCSFLMEAS         t        S-flow        S-plaq        E_plaq    dE_plaq/dt  "
                      "       E_clv     dE_clv/dt     Qtopo_clv         E_log     dE_log/dt     "
                      "Qtopo_log   [t_step_size]   [max_S-plaq]\n";
        first = false;
    }
    atype slocal = measure_cs_s(V) /
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
    atype deplaqdt = measure_dE_wplaq_dt(V, K) / lattice.volume();

    // derivative of clover energy density w.r.t. to flow time :
    atype declovdt = measure_dE_clov_dt(V, K) / lattice.volume();

    // derivative of log energy density w.r.t. to flow time :
    atype delogdt = measure_dE_log_dt(V, K) / lattice.volume();

    // print formatted results to standard output :
    hila::out0 << string_format("CSFLMEAS  % 9.3f % 0.6e % 0.6e % 0.6e % 0.6e % 0.6e % 0.6e % 0.6e "
                                "% 0.6e % 0.6e % 0.6e     [%0.3e]    [%0.3e]",
                                flow_t, slocal, plaq, eplaq, deplaqdt, ecl, declovdt, qtopocl, elog,
                                delogdt, qtopolog, t_step, max_plaq)
               << '\n';
}

template <typename group, typename atype = hila::arithmetic_type<group>>
atype do_cs_flow_adapt(GaugeField<group> &V, VectorField<Algebra<group>> &E,
                       atype t_len, atype atol = 1.0e-6, atype rtol = 1.0e-6, atype tstep = 0.0) {
    // constant action flow integration from flow scale l_start to l_end using 3rd order
    // 3-step Runge-Kutta (RK3) from arXiv:1006.4518 (cf. appendix C of
    // arXiv:2101.05320 for derivation of this Runge-Kutta method)
    // and embedded RK2 for adaptive step size

    atype esp = 3.0; // expected error scaling power: err ~ step^(esp)
                     //   - for global error: esp\approx 2.0
                     //   - for single step error: esp\approx 3.0
    atype iesp = 1.0 / esp;

    atype maxstepmf = 1.0e2; // max. growth factor of adaptive step size
    atype minmaxreldiff = pow(maxstepmf, -esp);

    // translate flow scale interval [l_start,l_end] to corresponding
    // flow time interval [t,tmax] :
    atype t = 0;
    atype tmax = t_len;

    atype ubstep = (tmax - t) / 3.0; // max. allowed time step

    atype tatol = atol * sqrt(2.0);

    // hila::out0<<"t: "<<t<<" , tmax: "<<tmax<<" , step: "<<tstep<<" , minmaxreldiff:
    // "<<minmaxreldiff<<"\n";

    // temporary variables :
    VectorField<Algebra<group>> k1, k2, tk;
    GaugeField<group> V2, V0;
    Field<atype> reldiff;
    atype maxstep;

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

    if (step == 0) {
        // when using a gauge action for cs flow that is different from
        // the one used to sample the gauge cofingurations, the initial force
        // can be huge. Therefore, if no t_step is provided as input, the inital
        // value for step is here adjustet so that
        // step * <largest local force> = maxstk
        atype maxstk = 1.0e-1;

        // get max. local gauge force:
        
        get_cs_force(V, E);
        tk = E;
        atype maxtk = 0.0;
        foralldir(d) {
            onsites(ALL) {
                reldiff[X] = (tk[d][X].squarenorm());
            }
            atype tmaxtk = reldiff.max();
            if(tmaxtk>maxtk) {
                maxtk = tmaxtk;
            }
        }
        maxtk = sqrt(0.5 * maxtk);

        if (step == 0) {
            if (maxtk > maxstk) {
                step = min(maxstk / maxtk,
                           ubstep); // adjust initial step size based on max. force magnitude
                hila::out0 << "CSFINFO using max. gauge force (max_X |F(X)|=" << maxtk
                           << ") to set initial flow time step size: " << step << "\n";
            } else {
                step = min((atype)1.0, ubstep);
            }
        } else if (step * maxtk > maxstk) {
            step = min(maxstk / maxtk,
                       ubstep); // adjust initial step size based on max. force magnitude
            hila::out0 << "CSFINFO using max. gauge force (max_X |F(X)|=" << maxtk
                       << ") to set initial flow time step size: " << step << "\n";
        }
    }


    V0 = V;
    bool stop = false;
    tstep = step;
    while (t < tmax && !stop) {
        if (t + step >= tmax) {
            step = tmax - t;
            stop = true;
        } else {
            tstep = step;
            if (t + 1.5 * step >= tmax) {
                step = 0.501 * (tmax - t);
            }
        }
        get_cs_force(V, E);
        k1 = E;
        foralldir(d) onsites(ALL) {
            // first steps of RK3 and RK2 are the same :
            V[d][X] = chexp(k1[d][X] * (step * a11)) * V[d][X];
        }
        get_cs_force(V, E);
        k2 = E;
        foralldir(d) onsites(ALL) {
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
        get_cs_force(V, E);
        k1 = E;
        foralldir(d) onsites(ALL) {
            // third step of RK3 :
            k1[d][X] *= (step * a33);
            k1[d][X] -= k2[d][X];
            V[d][X] = chexp(k1[d][X]) * V[d][X];
        }

        // determine maximum difference between RK3 and RK2,
        // relative to desired accuracy :
        atype maxreldiff = 0.0;
        foralldir(d) {
            onsites(ALL) {
                reldiff[X] = (V2[d][X] * V[d][X].dagger()).project_to_algebra().norm() /
                             (tatol + rtol * tk[d][X].norm() / step);
            }
            atype tmaxreldiff = reldiff.max();
            if (tmaxreldiff > maxreldiff) {
                maxreldiff = tmaxreldiff;
            }
        }

        if (maxreldiff < 1.0) {
            // proceed to next iteration
            t += step;
            V.reunitarize_gauge();
            V0 = V;
        } else {
            // repeat current iteration if step size was larger than maxstep
            V = V0;
            stop = false;
        }

        // determine step size to achieve desired accuracy goal :
        if (maxreldiff > minmaxreldiff) {
            maxstep = step * pow(maxreldiff, -iesp);
        } else {
            maxstep = step * maxstepmf;
        }

        // adjust step size :
        step = min((atype)0.9 * maxstep, ubstep);
        //hila::out0<<"t: "<<t<<" , maxreldiff: "<<maxreldiff<<" , maxstep: "
        //            <<maxstep<<" gf , step:"<<step<<"\n";
    }

    return tstep;
}

#endif