#ifndef __DIRAC_STAGGERED_H__
#define __DIRAC_STAGGERED_H__

#include "../plumbing/defs.h"
#include "../datatypes/cmplx.h"
#include "../datatypes/matrix.h"
#include "../datatypes/sun.h"
#include "../plumbing/field.h"
#include "../../libraries/hmc/gauge_field.h"

template <typename vector> Field<vector> staggered_dirac_temp[NDIM];

/// Initialize the staggered eta field
inline void init_staggered_eta(Field<double> (&staggered_eta)[NDIM]) {
    // Initialize the staggered eta field
    foralldir(d) {
        onsites(ALL) {
            element<CoordinateVector> l = X.coordinates();
            element<int> sumcoord = 0;
            for (int d2 = e_x; d2 < d; d2++) {
                sumcoord += l[d2];
            }
            // +1 if sumcoord divisible by 2, -1 otherwise
            // If statements not yet implemented for vectors
            staggered_eta[d][X] = (sumcoord % 2) * 2 - 1;
        }
    }
}

/// Apply the mass term v_out = m*v_in
template <typename vtype>
void dirac_staggered_diag(const double mass, const Field<vtype> &v_in,
                          Field<vtype> &v_out, Parity par) {
    v_out[par] = v_out[X] + mass * v_in[X];
}

/// Apply the inverse of the diagonal part,
/// v_out = 1/m * v_in
template <typename vtype>
void dirac_staggered_diag_inverse(const double mass, Field<vtype> &v_out, Parity par) {
    v_out[par] = (1.0 / mass) * v_out[X];
}

/// Apply the derivative part
template <typename mtype, typename vtype>
void dirac_staggered_hop(const Field<mtype> *gauge, const Field<vtype> &v_in,
                         Field<vtype> &v_out, Field<double> (&staggered_eta)[NDIM],
                         Parity par, int sign) {
    Field<vtype>(&vtemp)[NDIM] = staggered_dirac_temp<vtype>;
    foralldir(dir) {
        vtemp[dir].copy_boundary_condition(v_in);
        v_in.start_fetch(dir, par);
    }

    // First multiply the by conjugate before communicating the vector
    foralldir(dir) {
        vtemp[dir][opp_parity(par)] = gauge[dir][X].adjoint() * v_in[X];
        vtemp[dir].start_fetch(-dir, par);
    }

    // Run neighbour fetches and multiplications
    foralldir(dir) {
        v_out[par] = v_out[X] + 0.5 * sign * staggered_eta[dir][X] *
                                    (gauge[dir][X] * v_in[X + dir] - vtemp[dir][X - dir]);
    }
}

/// Calculate derivative  d/dA_x,mu (chi D psi)
/// Necessary for the HMC force calculation.
template <typename gaugetype, typename momtype, typename vtype>
void dirac_staggered_calc_force(const Field<gaugetype> *gauge, const Field<vtype> &chi,
                                const Field<vtype> &psi, Field<momtype> (&out)[NDIM],
                                Field<double> (&staggered_eta)[NDIM], int sign,
                                Parity par) {
    foralldir(dir) {
        out[dir][ALL] = 0;
        out[dir][par] =
            -sign * 0.5 * staggered_eta[dir][X] * chi[X + dir].outer_product(psi[X]);
        out[dir][opp_parity(par)] = out[dir][X] + sign * 0.5 *
                                                      staggered_eta[dir][X + dir] *
                                                      psi[X + dir].outer_product(chi[X]);
    }
}

/// An operator class that applies the staggered Dirac operator
/// D.apply(in, out) aplies the operator
/// D.dagger(int out) aplies the conjugate of the operator
///
/// This is useful for defining inverters as composite
/// operators. For example the conjugate gradient inverter
/// is CG<dirac_staggered>.
template <typename matrix> class dirac_staggered {
  private:
    /// The eta Field in the staggered operator, eta_x,\nu -1^(sum_mu<nu x_\mu)
    Field<double> staggered_eta[NDIM];

  public:
    /// the fermion mass
    double mass;
    /// The SU(N) vector type
    using vector_type = SU_vector<matrix::size, hila::number_type<matrix>>;
    /// The matrix type
    using matrix_type = matrix;
    /// A reference to the gauge links used in the dirac operator
    Field<matrix> (&gauge)[NDIM];

    /// Single precision type in case the base type is double precision.
    /// This is used to precondition the inversion of this operator
    using type_flt = dirac_staggered<typename gauge_field_base<matrix>::gauge_type_flt>;

    /// The parity this operator applies to
    Parity par = ALL;

    // Constructor: initialize mass, gauge and eta
    dirac_staggered(dirac_staggered &d) : gauge(d.gauge), mass(d.mass) {
        // Initialize the eta field (Share this?)
        init_staggered_eta(staggered_eta);
    }
    // Constructor: initialize mass, gauge and eta
    dirac_staggered(double m, Field<matrix> (&g)[NDIM]) : gauge(g), mass(m) {
        init_staggered_eta(staggered_eta);
    }
    // Constructor: initialize mass, gauge and eta
    dirac_staggered(double m, gauge_field_base<matrix> &g) : gauge(g.gauge), mass(m) {
        init_staggered_eta(staggered_eta);
    }

    /// Construct from another Dirac_Wilson operator of a different type.
    template <typename M>
    dirac_staggered(dirac_staggered<M> &d, gauge_field_base<matrix> &g)
        : gauge(g.gauge), mass(d.mass) {
        init_staggered_eta(staggered_eta);
    }

    /// Applies the operator to in
    void apply(const Field<vector_type> &in, Field<vector_type> &out) {
        out[ALL] = 0;
        dirac_staggered_diag(mass, in, out, ALL);
        dirac_staggered_hop(gauge, in, out, staggered_eta, ALL, 1);
    }

    /// Applies the conjugate of the operator
    void dagger(const Field<vector_type> &in, Field<vector_type> &out) {
        out[ALL] = 0;
        dirac_staggered_diag(mass, in, out, ALL);
        dirac_staggered_hop(gauge, in, out, staggered_eta, ALL, -1);
    }

    /// Applies the derivative of the Dirac operator with respect
    /// to the gauge field
    template <typename momtype>
    void force(const Field<vector_type> &chi, const Field<vector_type> &psi,
               Field<momtype> (&force)[NDIM], int sign = 1) {
        dirac_staggered_calc_force(gauge, chi, psi, force, staggered_eta, sign, ALL);
    }
};

/// Multiplying from the left applies the standard Dirac operator
template <typename vector, typename matrix>
Field<vector> operator*(dirac_staggered<matrix> D, const Field<vector> &in) {
    Field<vector> out;
    out.copy_boundary_condition(in);
    D.apply(in, out);
    return out;
}

/// Multiplying from the right applies the conjugate
template <typename vector, typename matrix>
Field<vector> operator*(const Field<vector> &in, dirac_staggered<vector> D) {
    Field<vector> out;
    out.copy_boundary_condition(in);
    D.dagger(in, out);
    return out;
}

/// An even-odd decomposed Wilson Dirac operator. Applies
/// D_{even to odd} D_{diag}^{-1} D_{odd to even} on the even
/// sites of the vector.
///
/// The fermion partition function is
///   det(D) = det(D_eveneodd) + det(D_{diag odd}).
/// Dirac_Wilson_evenodd can be used to replace D_Wilson
/// in the HMC action, as long as the diagonal odd to odd
/// part is accounted for.
///
/// This is useful for defining inverters as composite
/// operators. For example the conjugate gradient inverter
/// is CG<Dirac_Wilson_evenodd>.
///
/// As a side effect, the output Field becomes
/// out = D_{diag}^{-1} D_{odd to even} in
///
template <typename matrix> class dirac_staggered_evenodd {
  private:
    /// The eta Field in the staggered operator, eta_x,\nu -1^(sum_mu<nu x_\mu)
    Field<double> staggered_eta[NDIM];

    /// A reference to the gauge links used in the dirac operator
    Field<matrix> (&gauge)[NDIM];

  public:
    /// the fermion mass
    double mass;
    /// The SU(N) vector type
    using vector_type = SU_vector<matrix::size, typename matrix::base::type>;
    /// The matrix type
    using matrix_type = matrix;

    /// Single precision type in case the base type is double precision.
    /// This is used to precondition the inversion of this operator
    using type_flt =
        dirac_staggered_evenodd<typename gauge_field_base<matrix>::gauge_type_flt>;

    /// The parity this operator applies to
    Parity par = EVEN;

    /// Constructor: initialize mass, gauge and eta
    dirac_staggered_evenodd(dirac_staggered_evenodd &d) : gauge(d.gauge), mass(d.mass) {
        init_staggered_eta(staggered_eta);
    }
    /// Constructor: initialize mass, gauge and eta
    dirac_staggered_evenodd(double m, Field<matrix> (&U)[NDIM]) : gauge(U), mass(m) {
        init_staggered_eta(staggered_eta);
    }
    /// Constructor: initialize mass, gauge and eta
    dirac_staggered_evenodd(double m, gauge_field_base<matrix> &g)
        : gauge(g.gauge), mass(m) {
        init_staggered_eta(staggered_eta);
    }

    /// Construct from another Dirac_Wilson operator of a different type.
    template <typename M>
    dirac_staggered_evenodd(dirac_staggered_evenodd<M> &d, gauge_field_base<matrix> &g)
        : gauge(g.gauge), mass(d.mass) {
        init_staggered_eta(staggered_eta);
    }

    /// Applies the operator to in
    inline void apply(Field<vector_type> &in, Field<vector_type> &out) {
        out[ALL] = 0;
        dirac_staggered_diag(mass, in, out, EVEN);

        dirac_staggered_hop(gauge, in, out, staggered_eta, ODD, 1);
        dirac_staggered_diag_inverse(mass, out, ODD);
        dirac_staggered_hop(gauge, out, out, staggered_eta, EVEN, 1);
    }

    /// Applies the conjugate of the operator
    inline void dagger(Field<vector_type> &in, Field<vector_type> &out) {
        out[ALL] = 0;
        dirac_staggered_diag(mass, in, out, EVEN);

        dirac_staggered_hop(gauge, in, out, staggered_eta, ODD, -1);
        dirac_staggered_diag_inverse(mass, out, ODD);
        dirac_staggered_hop(gauge, out, out, staggered_eta, EVEN, -1);
    }

    /// Applies the derivative of the Dirac operator with respect
    /// to the gauge Field
    template <typename momtype>
    inline void force(const Field<vector_type> &chi, const Field<vector_type> &psi,
                      Field<momtype> (&force)[NDIM], int sign) {
        Field<momtype> force2[NDIM];
        Field<vector_type> tmp;
        tmp.copy_boundary_condition(chi);

        tmp[ALL] = 0;
        dirac_staggered_hop(gauge, chi, tmp, staggered_eta, ODD, -sign);
        dirac_staggered_diag_inverse(mass, tmp, ODD);
        dirac_staggered_calc_force(gauge, tmp, psi, force, staggered_eta, sign, EVEN);

        tmp[ALL] = 0;
        dirac_staggered_hop(gauge, psi, tmp, staggered_eta, ODD, sign);
        dirac_staggered_diag_inverse(mass, tmp, ODD);
        dirac_staggered_calc_force(gauge, chi, tmp, force2, staggered_eta, sign, ODD);

        foralldir(dir) { force[dir][ALL] = force[dir][X] + force2[dir][X]; }
    }
};

#endif