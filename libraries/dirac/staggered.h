#ifndef __DIRAC_STAGGERED_H__
#define __DIRAC_STAGGERED_H__

#include "../plumbing/defs.h"
#include "../datatypes/cmplx.h"
#include "../datatypes/matrix.h"
#include "../datatypes/sun.h"
#include "../plumbing/field.h"
#include "../../libraries/hmc/gauge_field.h"

template<typename vector>
field<vector> staggered_dirac_temp[NDIM];


/// Initialize the staggered eta field
inline void init_staggered_eta(field<double> (&staggered_eta)[NDIM]){
  // Initialize the staggered eta field
  foralldir(d){
    onsites(ALL){
      element<coordinate_vector> l = X.coordinates();
      element<int> sumcoord = 0;
      for(int d2=e_x;d2<d;d2++){
        sumcoord += l[d2];
      }
      // +1 if sumcoord divisible by 2, -1 otherwise
      // If statements not yet implemented for vectors
      staggered_eta[d][X] = (sumcoord%2)*2-1; 
    }
  }
}



/// Apply the mass term v_out = m*v_in
template<typename vtype>
void dirac_staggered_diag(
  const double mass,
  const field<vtype> &v_in,
  field<vtype> &v_out,
  parity par)
{
  v_out[par] = v_out[X] + mass*v_in[X];
}

/// Apply the inverse of the diagonal part,
/// v_out = 1/m * v_in
template<typename vtype>
void dirac_staggered_diag_inverse(
  const double mass,
  field<vtype> &v_out,
  parity par)
{
  v_out[par] = (1.0/mass) * v_out[X];
}


/// Apply the derivative part
template<typename mtype, typename vtype>
void dirac_staggered_hop(
  const field<mtype> *gauge,
  const field<vtype> &v_in,
  field<vtype> &v_out,
  field<double> (&staggered_eta)[NDIM],
  parity par, int sign)
{
  field<vtype> (&vtemp)[NDIM] = staggered_dirac_temp<vtype>;
  foralldir(dir){
    vtemp[dir].copy_boundary_condition(v_in);
    v_in.start_fetch(dir, par);
  }

  // First multiply the by conjugate before communicating the vector
  foralldir(dir){
    vtemp[dir][opp_parity(par)] = gauge[dir][X].adjoint()*v_in[X];
    vtemp[dir].start_fetch(-dir, par);
  }

  // Run neighbour fetches and multiplications
  foralldir(dir){
    v_out[par] = v_out[X] + 0.5 * sign * staggered_eta[dir][X] *
      ( gauge[dir][X]*v_in[X+dir] -  vtemp[dir][X-dir]);
  }
}


/// Calculate derivative  d/dA_x,mu (chi D psi)
/// Necessary for the HMC force calculation.
template<typename gaugetype, typename momtype, typename vtype>
void dirac_staggered_calc_force(
  const field<gaugetype> *gauge,
  const field<vtype> &chi,
  const field<vtype> &psi,
  field<momtype> (&out)[NDIM],
  field<double> (&staggered_eta)[NDIM],
  int sign, parity par)
{
  foralldir(dir){
    out[dir][ALL] = 0;
    out[dir][par] = -sign * 0.5 *
      staggered_eta[dir][X] * chi[X+dir].outer_product(psi[X]);
    out[dir][opp_parity(par)] = out[dir][X] + sign*0.5 *
      staggered_eta[dir][X+dir] * psi[X+dir].outer_product(chi[X]);
  }
}


/// An operator class that applies the staggered Dirac operator
/// D.apply(in, out) aplies the operator
/// D.dagger(int out) aplies the conjugate of the operator
///
/// This is useful for defining inverters as composite
/// operators. For example the conjugate gradient inverter 
/// is CG<dirac_staggered>.
template<typename matrix>
class dirac_staggered {
  private:
    /// The eta field in the staggered operator, eta_x,\nu -1^(sum_mu<nu x_\mu)
    field<double> staggered_eta[NDIM];

  public:
    /// the fermion mass
    double mass;
    /// The SU(N) vector type
    using vector_type = SU_vector<matrix::size, typename matrix::base_type>;
    /// The matrix type
    using matrix_type = matrix;
    /// A reference to the gauge links used in the dirac operator
    field<matrix> (&gauge)[NDIM];

    /// Single precision type in case the base type is double precision.
    /// This is used to precondition the inversion of this operator
    using type_flt = dirac_staggered<typename gauge_field_base<matrix>::gauge_type_flt>;

    /// The parity this operator applies to
    parity par = ALL;

    // Constructor: initialize mass, gauge and eta
    dirac_staggered(dirac_staggered &d) : gauge(d.gauge), mass(d.mass) {
      // Initialize the eta field (Share this?)
      init_staggered_eta(staggered_eta);
    }
    // Constructor: initialize mass, gauge and eta
    dirac_staggered(double m, field<matrix> (&g)[NDIM]) : gauge(g), mass(m)  {
      init_staggered_eta(staggered_eta);
    }
    // Constructor: initialize mass, gauge and eta
    dirac_staggered(double m, gauge_field_base<matrix> &g) : gauge(g.gauge), mass(m) {
      init_staggered_eta(staggered_eta);
    }

    /// Construct from another Dirac_Wilson operator of a different type.
    template<typename M>
    dirac_staggered(dirac_staggered<M> &d, gauge_field_base<matrix> &g) : gauge(g.gauge), mass(d.mass) 
    {
      init_staggered_eta(staggered_eta);
    }


    /// Applies the operator to in
    void apply( const field<vector_type> & in, field<vector_type> & out){
      out[ALL] = 0;
      dirac_staggered_diag(mass, in, out, ALL);
      dirac_staggered_hop(gauge, in, out, staggered_eta, ALL, 1);
    }

    /// Applies the conjugate of the operator
    void dagger( const field<vector_type> & in, field<vector_type> & out){
      out[ALL] = 0;
      dirac_staggered_diag(mass, in, out, ALL);
      dirac_staggered_hop(gauge, in, out, staggered_eta, ALL, -1);
    }

    /// Applies the derivative of the Dirac operator with respect
    /// to the gauge field
    template<typename momtype>
    void force( const field<vector_type> & chi, const field<vector_type> & psi, field<momtype> (&force)[NDIM], int sign=1){
      dirac_staggered_calc_force(gauge, chi, psi, force, staggered_eta, sign, ALL);
    }
};



/// Multiplying from the left applies the standard Dirac operator
template<typename vector, typename matrix>
field<vector> operator* (dirac_staggered<matrix> D, const field<vector> & in) {
  field<vector> out;
  out.copy_boundary_condition(in);
  D.apply(in, out);
  return out;
}

/// Multiplying from the right applies the conjugate
template<typename vector, typename matrix>
field<vector> operator* (const field<vector> & in, dirac_staggered<vector> D) {
  field<vector> out;
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
/// As a side effect, the output field becomes 
/// out = D_{diag}^{-1} D_{odd to even} in
///
template<typename matrix>
class dirac_staggered_evenodd {
  private:
    /// The eta field in the staggered operator, eta_x,\nu -1^(sum_mu<nu x_\mu)
    field<double> staggered_eta[NDIM];

    /// A reference to the gauge links used in the dirac operator
    field<matrix> (&gauge)[NDIM];
  public:
    /// the fermion mass
    double mass;
    /// The SU(N) vector type
    using vector_type = SU_vector<matrix::size, typename matrix::base_type>;
    /// The matrix type
    using matrix_type = matrix;

    /// Single precision type in case the base type is double precision.
    /// This is used to precondition the inversion of this operator
    using type_flt = dirac_staggered_evenodd<typename gauge_field_base<matrix>::gauge_type_flt>;

    /// The parity this operator applies to
    parity par = EVEN;

    /// Constructor: initialize mass, gauge and eta
    dirac_staggered_evenodd(dirac_staggered_evenodd &d) : gauge(d.gauge), mass(d.mass) {
      init_staggered_eta(staggered_eta);
    }
    /// Constructor: initialize mass, gauge and eta
    dirac_staggered_evenodd(double m, field<matrix> (&U)[NDIM]) : gauge(U), mass(m) {
      init_staggered_eta(staggered_eta);
    }
    /// Constructor: initialize mass, gauge and eta
    dirac_staggered_evenodd(double m, gauge_field_base<matrix> &g) : gauge(g.gauge), mass(m) {
      init_staggered_eta(staggered_eta);
    }

    /// Construct from another Dirac_Wilson operator of a different type.
    template<typename M>
    dirac_staggered_evenodd(dirac_staggered_evenodd<M> &d, gauge_field_base<matrix> &g) : gauge(g.gauge), mass(d.mass) {
      init_staggered_eta(staggered_eta);
    }



    /// Applies the operator to in
    inline void apply(  field<vector_type> & in, field<vector_type> & out){
      out[ALL] = 0;
      dirac_staggered_diag(mass, in, out, EVEN);

      dirac_staggered_hop(gauge, in, out, staggered_eta, ODD, 1);
      dirac_staggered_diag_inverse(mass, out, ODD);
      dirac_staggered_hop(gauge, out, out, staggered_eta, EVEN, 1);
    }

    /// Applies the conjugate of the operator
    inline void dagger( field<vector_type> & in, field<vector_type> & out){
      out[ALL] = 0;
      dirac_staggered_diag(mass, in, out, EVEN);

      dirac_staggered_hop(gauge, in, out, staggered_eta, ODD, -1);
      dirac_staggered_diag_inverse(mass, out, ODD);
      dirac_staggered_hop(gauge, out, out, staggered_eta, EVEN, -1);
    }

    /// Applies the derivative of the Dirac operator with respect
    /// to the gauge field
    template<typename momtype>
    inline void force(const field<vector_type> & chi, const field<vector_type> & psi, field<momtype> (&force)[NDIM], int sign){
      field<momtype> force2[NDIM];
      field<vector_type> tmp;
      tmp.copy_boundary_condition(chi);

      tmp[ALL] = 0;
      dirac_staggered_hop(gauge, chi, tmp, staggered_eta, ODD, -sign);
      dirac_staggered_diag_inverse(mass, tmp, ODD);
      dirac_staggered_calc_force(gauge, tmp, psi, force, staggered_eta, sign, EVEN);

      tmp[ALL] = 0;
      dirac_staggered_hop(gauge, psi, tmp, staggered_eta, ODD, sign);
      dirac_staggered_diag_inverse(mass, tmp, ODD);
      dirac_staggered_calc_force(gauge, chi, tmp, force2, staggered_eta, sign, ODD);

      foralldir(dir){
        force[dir][ALL] = force[dir][X] + force2[dir][X];
      }
    }
};



#endif