#ifndef FERMION_FIELD_H
#define FERMION_FIELD_H

#include "gauge_field.h"
#include "dirac/Hasenbusch.h"
#include "dirac/conjugate_gradient.h"
#include "MRE_guess.h"
#include <cmath>

/// Define the action of a pseudofermion for HMC
///
/// Implements methods for calculating the current action
/// and the force (derivative with respect to the gauge
/// field).
///
/// Includes an implementation of the MRE initial guess,
/// which is calculated in the base of a few previous
/// solutions. Using this requires a higher accuracy,
/// since the initial guess is not time reversible.
///
template <typename gauge_field, typename DIRAC_OP>
class fermion_action : public action_base {
  public:
    using vector_type = typename DIRAC_OP::vector_type;
    using momtype = SquareMatrix<gauge_field::N, Cmplx<typename gauge_field::basetype>>;
    gauge_field &gauge;
    DIRAC_OP &D;
    Field<vector_type> chi;

    /// We save a few previous invertions to build an initial guess.
    /// old_chi contains a list of these
    int MRE_size = 0;
    std::vector<Field<vector_type>> old_chi_inv;

    void setup(int mre_guess_size) {
#if NDIM > 3
        chi.set_boundary_condition(e_t, BoundaryCondition::ANTIPERIODIC);
        chi.set_boundary_condition(-e_t, BoundaryCondition::ANTIPERIODIC);
#endif
        MRE_size = mre_guess_size;
        old_chi_inv.resize(MRE_size);
        for (int i = 0; i < MRE_size; i++) {
            old_chi_inv[i][ALL] = 0;
        }
    }

    fermion_action(DIRAC_OP &d, gauge_field &g) : D(d), gauge(g) {
        chi = 0.0; // Allocates chi and sets it to zero
        setup(0);
    }

    fermion_action(DIRAC_OP &d, gauge_field &g, int mre_guess_size) : D(d), gauge(g) {
        chi = 0.0; // Allocates chi and sets it to zero
        setup(mre_guess_size);
    }

    fermion_action(fermion_action &fa) : gauge(fa.gauge), D(fa.D) {
        chi = fa.chi; // Copies the field
        setup(fa.MRE_size);
    }

    /// Build an initial guess for the fermion matrix inversion
    /// by inverting first in the limited space of a few previous
    /// solutions. These are saved in old_chi.
    void initial_guess(Field<vector_type> &chi, Field<vector_type> &psi) {
        psi[ALL] = 0;
        if (MRE_size > 0) {
            MRE_guess(psi, chi, D, old_chi_inv);
        }
        // If the gauge type is double precision, solve first in single precision
        if constexpr (std::is_same<double, typename gauge_field::basetype>::value) {
            output0 << "Starting with single precision inversion\n";

            auto single_precision = gauge.get_single_precision();
            typename DIRAC_OP::type_flt D_flt(D, single_precision);
            Field<typename DIRAC_OP::type_flt::vector_type> c, p, t1, t2;
            c[ALL] = chi[X];
            p[ALL] = psi[X];
            CG inverse(D_flt);
            inverse.apply(c, p);

            D_flt.apply(p, t1);
            D_flt.dagger(t1, t2);
            psi[ALL] = p[X];
        }
    }

    /// Return the value of the action with the current
    /// field configuration
    double action() {
        Field<vector_type> psi;
        psi.copy_boundary_condition(chi);
        CG<DIRAC_OP> inverse(D);
        double action = 0;

        gauge.refresh();

        psi = 0;
        initial_guess(chi, psi);
        inverse.apply(chi, psi);
        onsites(D.par) { action += chi[X].rdot(psi[X]); }
        return action;
    }

    /// Calculate the action as a field of double precision numbers
    void action(Field<double> &S) {
        Field<vector_type> psi;
        psi.copy_boundary_condition(chi);
        CG<DIRAC_OP> inverse(D);

        gauge.refresh();

        psi = 0;
        initial_guess(chi, psi);
        inverse.apply(chi, psi);
        onsites(D.par) {
            if (disable_avx[X] == 0) {
            };
            S[X] += chi[X].rdot(psi[X]);
        }
    }

    /// Generate a pseudofermion field with a distribution given
    /// by the action chi 1/(D_dagger D) chi
    void draw_gaussian_fields() {
        Field<vector_type> psi;
        psi.copy_boundary_condition(chi);
        gauge.refresh();

        onsites(D.par) {
            if (disable_avx[X] == 0) {
            };
            psi[X].gaussian();
        }
        D.dagger(psi, chi);
    }

    /// Add new solution to the list for MRE
    void save_new_solution(Field<vector_type> &psi) {
        if (MRE_size > 0) {
            for (int i = 1; i < MRE_size; i++) {
                old_chi_inv[i] = old_chi_inv[i - 1];
            }
            old_chi_inv[0] = psi;
        }
    }

    /// Update the momentum with the derivative of the fermion
    /// action
    void force_step(double eps) {
        Field<vector_type> psi, Mpsi;
        psi.copy_boundary_condition(chi);
        Mpsi.copy_boundary_condition(chi);
        Field<momtype> force[NDIM], force2[NDIM];

        CG<DIRAC_OP> inverse(D);
        gauge.refresh();

        output0 << "base force\n";
        initial_guess(chi, psi);
        inverse.apply(chi, psi);
        save_new_solution(psi);

        D.apply(psi, Mpsi);

        D.force(Mpsi, psi, force, 1);
        D.force(psi, Mpsi, force2, -1);

        foralldir(dir) { force[dir][ALL] = -eps * (force[dir][X] + force2[dir][X]); }
        gauge.add_momentum(force);
    }
};

/// The Hasenbusch method for updating fermion fields:
/// Split the Dirac determinant into two parts,
/// D_h1 = D + mh and
/// D_h2 = D * 1 / (D + mh)^dagger
///
///
/// This is the first action term, with D_h1 = D + mh.
/// Since the only real difference here is an addition
/// to the original operator, we can use fermion_action
/// with a different operator.
///
template <typename gauge_field, typename DIRAC_OP>
class Hasenbusch_action_1 : public action_base {
  public:
    Hasenbusch_operator<DIRAC_OP> D_h;
    fermion_action<gauge_field, Hasenbusch_operator<DIRAC_OP>> base_action;
    double _mh;

    Hasenbusch_action_1(DIRAC_OP &d, gauge_field &g, double mh)
        : _mh(mh), D_h(d, mh), base_action(D_h, g) {}
    Hasenbusch_action_1(DIRAC_OP &d, gauge_field &g, double mh, int mre_guess_size)
        : _mh(mh), D_h(d, mh), base_action(D_h, g, mre_guess_size) {}

    Hasenbusch_action_1(Hasenbusch_action_1 &fa)
        : _mh(fa._mh), D_h(fa.D_h), base_action(fa.base_action) {}

    double action() { return (base_action.action()); }
    void action(Field<double> &S) { base_action.action(S); }
    void draw_gaussian_fields() { base_action.draw_gaussian_fields(); }
    void force_step(double eps) { base_action.force_step(eps); }
};

/// The second Hasenbusch action term, D_h2 = D/(D^dagger + mh).
/// The force and action of the second term are significantly
/// different from the standard fermion action and are implemented
/// here.
template <typename gauge_field, typename DIRAC_OP>
class Hasenbusch_action_2 : public action_base {
  public:
    using vector_type = typename DIRAC_OP::vector_type;
    using momtype = SquareMatrix<gauge_field::N, Cmplx<typename gauge_field::basetype>>;
    gauge_field &gauge;
    DIRAC_OP D;
    Hasenbusch_operator<DIRAC_OP> D_h;
    double mh;
    Field<vector_type> chi;

    // We save a few previous invertions to build an initial guess.
    // old_chi contains a list of these
    int MRE_size = 0;
    std::vector<Field<vector_type>> old_chi_inv;

    void setup(int mre_guess_size) {
#if NDIM > 3
        chi.set_boundary_condition(e_t, BoundaryCondition::ANTIPERIODIC);
        chi.set_boundary_condition(-e_t, BoundaryCondition::ANTIPERIODIC);
#endif
        MRE_size = mre_guess_size;
        old_chi_inv.resize(MRE_size);
        for (int i = 0; i < MRE_size; i++) {
            old_chi_inv[i][ALL] = 0;
        }
    }

    Hasenbusch_action_2(DIRAC_OP &d, gauge_field &g, double _mh)
        : mh(_mh), D(d), D_h(d, _mh), gauge(g) {
        chi = 0.0; // Allocates chi and sets it to zero
        setup(0);
    }
    Hasenbusch_action_2(DIRAC_OP &d, gauge_field &g, double _mh, int mre_guess_size)
        : mh(_mh), D(d), D_h(d, _mh), gauge(g) {
        chi = 0.0; // Allocates chi and sets it to zero
        setup(mre_guess_size);
    }

    Hasenbusch_action_2(Hasenbusch_action_2 &fa)
        : mh(fa.mh), D(fa.D), D_h(fa.D_h), gauge(fa.gauge) {
        chi = fa.chi; // Copies the field
        setup(fa.MRE_size);
    }

    /// Return the value of the action with the current
    /// Field configuration
    double action() {
        Field<vector_type> psi;
        Field<vector_type> v;
        psi.copy_boundary_condition(chi);
        v.copy_boundary_condition(chi);
        double action = 0;

        gauge.refresh();
        CG<DIRAC_OP> inverse(D);

        v[ALL] = 0;
        D_h.dagger(chi, psi);
        inverse.apply(psi, v);
        D.apply(v, psi);
        onsites(D.par) {
            if (disable_avx[X] == 0) {
            };
            action += norm_squared(psi[X]);
        }

        return action;
    }

    /// Return the action as a field of double precision numbers
    void action(Field<double> &S) {
        Field<vector_type> psi;
        Field<vector_type> v;
        psi.copy_boundary_condition(chi);
        v.copy_boundary_condition(chi);

        gauge.refresh();
        CG inverse(D);

        v[ALL] = 0;
        D.dagger(chi, psi);
        inverse.apply(psi, v);
        D.apply(v, psi);
        onsites(EVEN) { S[X] += norm_squared(psi[X]); }
    }

    /// Generate a pseudofermion field with a distribution given
    /// by the action chi 1/(D_dagger D) chi
    void draw_gaussian_fields() {
        Field<vector_type> psi;
        Field<vector_type> v;
        psi.copy_boundary_condition(chi);
        v.copy_boundary_condition(chi);
        CG inverse_h(D_h); // Applies 1/(D_h^dagger D_h)
        gauge.refresh();

        psi[ALL] = 0;
        onsites(D.par) {
            if (disable_avx[X] == 0) {
            };
            psi[X].gaussian();
        }
        D.dagger(psi, v);
        psi[ALL] = 0;
        inverse_h.apply(v, psi);
        D_h.apply(psi, chi);
    }

    /// Build an initial guess for the fermion matrix inversion
    /// by inverting first in the limited space of a few previous
    /// solutions. These are saved in old_chi.
    void initial_guess(Field<vector_type> &chi, Field<vector_type> &psi) {
        psi[ALL] = 0;
        if (MRE_size > 0) {
            MRE_guess(psi, chi, D, old_chi_inv);
        }
        // If the gauge type is double precision, solve first in single precision
        if constexpr (std::is_same<double, typename gauge_field::basetype>::value) {
            output0 << "Starting with single precision inversion\n";

            auto single_precision = gauge.get_single_precision();
            typename DIRAC_OP::type_flt D_flt(D, single_precision);
            Field<typename DIRAC_OP::type_flt::vector_type> c, p, t1, t2;
            c[ALL] = chi[X];
            p[ALL] = psi[X];
            CG inverse(D_flt);
            inverse.apply(c, p);

            D_flt.apply(p, t1);
            D_flt.dagger(t1, t2);
            psi[ALL] = p[X];
        }
    }

    /// Add new solution to the list
    void save_new_solution(Field<vector_type> &psi) {
        if (MRE_size > 0) {
            for (int i = 1; i < MRE_size; i++) {
                old_chi_inv[i] = old_chi_inv[i - 1];
            }
            old_chi_inv[0] = psi;
        }
    }

    /// Update the momentum with the derivative of the fermion
    /// action
    void force_step(double eps) {
        Field<vector_type> psi, Mpsi;
        Field<vector_type> Dhchi;
        psi.copy_boundary_condition(chi);
        Mpsi.copy_boundary_condition(chi);
        Dhchi.copy_boundary_condition(chi);
        Field<momtype> force[NDIM], force2[NDIM];

        CG<DIRAC_OP> inverse(D);
        gauge.refresh();

        D_h.dagger(chi, Dhchi);

        initial_guess(Dhchi, psi);
        inverse.apply(Dhchi, psi);
        save_new_solution(psi);

        D.apply(psi, Mpsi);

        Mpsi[D.par] = Mpsi[X] - chi[X];

        D.force(Mpsi, psi, force, 1);
        D.force(psi, Mpsi, force2, -1);

        foralldir(dir) { force[dir][ALL] = -eps * (force[dir][X] + force2[dir][X]); }
        gauge.add_momentum(force);
    }
};

#endif