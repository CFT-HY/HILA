#ifndef __DIRAC_WILSON_H__
#define __DIRAC_WILSON_H__

#include "plumbing/defs.h"
#include "datatypes/cmplx.h"
#include "datatypes/matrix.h"
#include "datatypes/wilson_vector.h"
#include "plumbing/field.h"
#include "hmc/gauge_field.h"

template <int N, typename radix>
Field<half_Wilson_vector<N, radix>> wilson_dirac_temp_vector[2 * NDIM];

/// Apply the hopping term to v_out and add to v_in
template <int N, typename radix, typename matrix>
inline void Dirac_Wilson_hop(const Field<matrix> *gauge, const double kappa,
                             const Field<Wilson_vector<N, radix>> &v_in,
                             Field<Wilson_vector<N, radix>> &v_out, parity par,
                             int sign) {
    Field<half_Wilson_vector<N, radix>>(&vtemp)[2 * NDIM] =
        wilson_dirac_temp_vector<N, radix>;
    for (int dir = 0; dir < 2 * NDIM; dir++) {
        vtemp[dir].copy_boundary_condition(v_in);
    }

    // Run neighbour fetches and multiplications
    foralldir(dir) {
        // First multiply the by conjugate before communicating
        onsites(opp_parity(par)) {
            half_Wilson_vector<N, radix> h(v_in[X], dir, -sign);
            vtemp[-dir][X] = gauge[dir][X].adjoint() * h;
        }
        onsites(opp_parity(par)) {
            half_Wilson_vector<N, radix> h(v_in[X], dir, sign);
            vtemp[dir][X] = h;
        }

        vtemp[dir].start_fetch(dir, par);
        vtemp[-dir].start_fetch(-dir, par);
    }

    // Calculate the derivatives. This
    foralldir(dir) {
        onsites(par) {
            v_out[X] = v_out[X] -
                       (kappa * gauge[dir][X] * vtemp[dir][X + dir]).expand(dir, sign) -
                       (kappa * vtemp[dir][X + dir]).expand(dir, -sign);
        }
    }
}

/// Apply the hopping term to v_in and overwrite v_out
template <int N, typename radix, typename matrix>
inline void Dirac_Wilson_hop_set(const Field<matrix> *gauge, const double kappa,
                                 const Field<Wilson_vector<N, radix>> &v_in,
                                 Field<Wilson_vector<N, radix>> &v_out, parity par,
                                 int sign) {
    Field<half_Wilson_vector<N, radix>>(&vtemp)[2 * NDIM] =
        wilson_dirac_temp_vector<N, radix>;
    for (int dir = 0; dir < 2 * NDIM; dir++) {
        vtemp[dir].copy_boundary_condition(v_in);
    }

    // Run neighbour fetches and multiplications
    foralldir(dir) {
        // First multiply the by conjugate before communicating
        onsites(opp_parity(par)) {
            half_Wilson_vector<N, radix> h(v_in[X], dir, -sign);
            vtemp[-dir][X] = gauge[dir][X].adjoint() * h;
        }
        onsites(opp_parity(par)) {
            half_Wilson_vector<N, radix> h(v_in[X], dir, sign);
            vtemp[dir][X] = h;
        }

        vtemp[dir].start_fetch(dir, par);
        vtemp[-dir].start_fetch(-dir, par);
    }
    // Set on first direction
    direction dir = direction(0);
    onsites(par) {
        v_out[X] = -(kappa * gauge[dir][X] * vtemp[dir][X + dir]).expand(dir, sign) -
                   (kappa * vtemp[dir][X + dir]).expand(dir, -sign);
    }
    // Add for all other directions
    for (int d = 1; d < NDIM; d++) {
        direction dir = direction(d);
        onsites(par) {
            half_Wilson_vector<N, radix> h1(v_in[X + dir], dir, sign);
            v_out[X] = v_out[X] -
                       (kappa * gauge[dir][X] * vtemp[dir][X + dir]).expand(dir, sign) -
                       (kappa * vtemp[dir][X + dir]).expand(dir, -sign);
        }
    }
}

/// The diagonal part of the operator. Without clover this is just the identity
template <int N, typename radix>
inline void Dirac_Wilson_diag(const Field<Wilson_vector<N, radix>> &v_in,
                              Field<Wilson_vector<N, radix>> &v_out, parity par) {
    v_out[par] = v_in[X];
}

/// Inverse of the diagonal part. Without clover this does nothing.
template <int N, typename radix>
inline void Dirac_Wilson_diag_inverse(Field<Wilson_vector<N, radix>> &v, parity par) {}

/// Calculate derivative  d/dA_x,mu (chi D psi)
/// Necessary for the HMC force calculation.
template <int N, typename radix, typename gaugetype, typename momtype>
inline void Dirac_Wilson_calc_force(const Field<gaugetype> *gauge, const double kappa,
                                    const Field<Wilson_vector<N, radix>> &chi,
                                    const Field<Wilson_vector<N, radix>> &psi,
                                    Field<momtype> (&out)[NDIM], parity par, int sign) {
    Field<half_Wilson_vector<N, radix>>(&vtemp)[2 * NDIM] =
        wilson_dirac_temp_vector<N, radix>;
    vtemp[0].copy_boundary_condition(chi);
    vtemp[1].copy_boundary_condition(chi);

    foralldir(dir) {
        onsites(opp_parity(par)) {
            half_Wilson_vector<N, radix> hw(chi[X], dir, -sign);
            vtemp[0][X] = hw;
        }
        onsites(par) {
            half_Wilson_vector<N, radix> hw(psi[X], dir, sign);
            vtemp[1][X] = hw;
        }

        out[dir][ALL] = 0;
        out[dir][par] =
            -kappa * ((vtemp[0][X + dir].expand(dir, -sign)).outer_product(psi[X]));
        out[dir][opp_parity(par)] =
            out[dir][X] -
            kappa * ((vtemp[1][X + dir].expand(dir, sign)).outer_product(chi[X]));
    }
}

/// An operator class that applies the Wilson Dirac operator
/// D.apply(in, out) aplies the operator
/// D.dagger(int out) aplies the conjugate of the operator
///
/// This is useful for defining inverters as composite
/// operators. For example the conjugate gradient inverter
/// is CG<Dirac_Wilson>.
template <typename matrix> class Dirac_Wilson {
  private:
    /// A reference to the gauge links used in the dirac operator
    Field<matrix> (&gauge)[NDIM];

  public:
    /// The hopping parameter, kappa = 1/(8-2m)
    double kappa;
    /// Size of the gauge matrix and color dimension of the Wilson vector
    static constexpr int N = matrix::size;

    using radix = number_type<matrix>;
    /// The wilson vector type
    using vector_type = Wilson_vector<N, radix>;
    /// The matrix type
    using matrix_type = matrix;

    /// Single precision type in case the base type is double precision.
    /// This is used to precondition the inversion of this operator
    using type_flt = Dirac_Wilson<typename gauge_field_base<matrix>::gauge_type_flt>;

    /// The parity this operator applies to
    parity par = ALL;

    /// Constructor: initialize mass and gauge
    Dirac_Wilson(Dirac_Wilson &d) : gauge(d.gauge), kappa(d.kappa) {}
    /// Constructor: initialize mass and gauge
    Dirac_Wilson(double k, Field<matrix> (&U)[NDIM]) : gauge(U), kappa(k) {}
    /// Constructor: initialize mass and gauge
    Dirac_Wilson(double k, gauge_field_base<matrix> &g) : gauge(g.gauge), kappa(k) {}

    /// Construct from another Dirac_Wilson operator of a different type.
    template <typename M>
    Dirac_Wilson(Dirac_Wilson<M> &d, gauge_field_base<matrix> &g)
        : gauge(g.gauge), kappa(d.kappa) {}

    /// Applies the operator to in
    inline void apply(const Field<vector_type> &in, Field<vector_type> &out) {
        Dirac_Wilson_diag(in, out, ALL);
        Dirac_Wilson_hop(gauge, kappa, in, out, ALL, 1);
    }

    /// Applies the conjugate of the operator
    inline void dagger(const Field<vector_type> &in, Field<vector_type> &out) {
        Dirac_Wilson_diag(in, out, ALL);
        Dirac_Wilson_hop(gauge, kappa, in, out, ALL, -1);
    }

    /// Applies the derivative of the Dirac operator with respect
    /// to the gauge Field
    template <typename momtype>
    inline void force(const Field<vector_type> &chi, const Field<vector_type> &psi,
                      Field<momtype> (&force)[NDIM], int sign = 1) {
        Dirac_Wilson_calc_force(gauge, kappa, chi, psi, force, ALL, sign);
    }
};

/// Multiplying from the left applies the standard Dirac operator
template <int N, typename radix, typename matrix>
Field<Wilson_vector<N, radix>> operator*(Dirac_Wilson<matrix> D,
                                         const Field<Wilson_vector<N, radix>> &in) {
    Field<Wilson_vector<N, radix>> out;
    D.apply(in, out);
    return out;
}

/// Multiplying from the right applies the conjugate
template <int N, typename radix, typename matrix>
Field<Wilson_vector<N, radix>> operator*(const Field<Wilson_vector<N, radix>> &in,
                                         Dirac_Wilson<matrix> D) {
    Field<Wilson_vector<N, radix>> out;
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
template <typename matrix> class Dirac_Wilson_evenodd {
  private:
    /// A reference to the gauge links used in the dirac operator
    Field<matrix> (&gauge)[NDIM];

  public:
    /// The hopping parameter, kappa = 1/(8-2m)
    double kappa;
    /// Size of the gauge matrix and color dimension of the Wilson vector
    static constexpr int N = matrix::size;
    /// The wilson vector type
    using radix = number_type<matrix>;
    using vector_type = Wilson_vector<N, radix>;
    /// The matrix type
    using matrix_type = matrix;

    /// Single precision type in case the base type is double precision.
    /// This is used to precondition the inversion of this operator
    using type_flt =
        Dirac_Wilson_evenodd<typename gauge_field_base<matrix>::gauge_type_flt>;

    /// The parity this operator applies to
    parity par = EVEN;

    /// Constructor: initialize mass and gauge
    Dirac_Wilson_evenodd(Dirac_Wilson_evenodd &d) : gauge(d.gauge), kappa(d.kappa) {}
    /// Constructor: initialize mass and gauge
    Dirac_Wilson_evenodd(double k, Field<matrix> (&U)[NDIM]) : gauge(U), kappa(k) {}
    /// Constructor: initialize mass and gauge
    Dirac_Wilson_evenodd(double k, gauge_field_base<matrix> &g)
        : gauge(g.gauge), kappa(k) {}

    /// Construct from another Dirac_Wilson operator of a different type.
    template <typename M>
    Dirac_Wilson_evenodd(Dirac_Wilson_evenodd<M> &d, gauge_field_base<matrix> &g)
        : gauge(g.gauge), kappa(d.kappa) {}

    /// Applies the operator to in
    inline void apply(const Field<vector_type> &in, Field<vector_type> &out) {
        Dirac_Wilson_diag(in, out, EVEN);

        Dirac_Wilson_hop_set(gauge, kappa, in, out, ODD, 1);
        Dirac_Wilson_diag_inverse(out, ODD);
        Dirac_Wilson_hop(gauge, -kappa, out, out, EVEN, 1);
        out[ODD] = 0;
    }

    /// Applies the conjugate of the operator
    inline void dagger(const Field<vector_type> &in, Field<vector_type> &out) {
        Dirac_Wilson_diag(in, out, EVEN);

        Dirac_Wilson_hop_set(gauge, kappa, in, out, ODD, -1);
        Dirac_Wilson_diag_inverse(out, ODD);
        Dirac_Wilson_hop(gauge, -kappa, out, out, EVEN, -1);
        out[ODD] = 0;
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
        Dirac_Wilson_hop_set(gauge, kappa, chi, tmp, ODD, -sign);
        Dirac_Wilson_diag_inverse(tmp, ODD);
        Dirac_Wilson_calc_force(gauge, -kappa, tmp, psi, force, EVEN, sign);

        tmp[ALL] = 0;
        Dirac_Wilson_hop_set(gauge, kappa, psi, tmp, ODD, sign);
        Dirac_Wilson_diag_inverse(tmp, ODD);
        Dirac_Wilson_calc_force(gauge, -kappa, chi, tmp, force2, ODD, sign);

        foralldir(dir) force[dir][ALL] = force[dir][X] + force2[dir][X];
    }
};

#endif