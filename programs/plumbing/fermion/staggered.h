#ifndef __DIRAC_STAGGERED_H__
#define __DIRAC_STAGGERED_H__

#include "../../plumbing/defs.h"
#include "../../datatypes/cmplx.h"
#include "../../datatypes/general_matrix.h"
#include "../../datatypes/wilson_vector.h"
#include "../../plumbing/field.h"


template<typename vector>
field<vector> staggered_dirac_temp[NDIM];


inline void init_staggered_eta(field<double> (&staggered_eta)[NDIM]){
  // Initialize the staggered eta field
  foralldir(d){
    onsites(ALL){
      element<coordinate_vector> l = X.coordinates();
      element<int> sumcoord = 0;
      for(int d2=XUP;d2<d;d2++){
        sumcoord += l[d2];
      }
      // +1 if sumcoord divisible by 2, -1 otherwise
      // If statements not yet implemented for vectors
      staggered_eta[d][X] = (sumcoord%2)*2-1; 
    }
  }
}



// Apply the mass diagonally
template<typename vtype>
void dirac_staggered_diag(
  const double mass,
  const field<vtype> &v_in,
  field<vtype> &v_out,
  parity par)
{
  v_out[par] = v_out[X] + mass*v_in[X];
}

// Apply the mass diagonally
template<typename vtype>
void dirac_staggered_diag_inverse(
  const double mass,
  field<vtype> &v_out,
  parity par)
{
  v_out[par] = (1.0/mass) * v_out[X];
}


template<typename mtype, typename vtype>
void dirac_staggered_hop(
  const field<mtype> *gauge,
  const field<vtype> &v_in,
  field<vtype> &v_out,
  field<double> (&staggered_eta)[NDIM],
  parity par, int sign)
{
  field<vtype> (&vtemp)[NDIM] = staggered_dirac_temp<vtype>;
  foralldir(dir)
    vtemp[dir].copy_boundary_condition(v_in);

  // First multiply the by conjugate before communicating the vector
  foralldir(dir){
    direction odir = opp_dir( (direction)dir );
    vtemp[dir][opp_parity(par)] = gauge[dir][X].conjugate()*v_in[X];
    //vtemp[dir].start_get(odir);
  }

  // Run neighbour fetches and multiplications
  foralldir(dir){
    direction odir = opp_dir( (direction)dir );
    v_out[par] = v_out[X] + 0.5 * sign * staggered_eta[dir][X] *
      ( gauge[dir][X]*v_in[X+dir] -  vtemp[dir][X-dir]);

  }
}



template<typename mtype, typename vtype>
void dirac_staggered_calc_force(
  const field<mtype> *gauge,
  const field<vtype> &chi,
  const field<vtype> &psi,
  field<mtype> (&out)[NDIM],
  field<double> (&staggered_eta)[NDIM],
  int sign, parity par)
{
  foralldir(dir){
    out[dir][ALL] = 0;
    out[dir][par] = -sign * 0.5 *
      staggered_eta[dir][X] * chi[X+dir].outer_product(psi[X]);
    out[dir][opp_parity(par)] = out[dir][X] + sign*0.5 *
      staggered_eta[dir][X+dir] * psi[X+dir].outer_product(chi[X]);

    out[dir][ALL] = gauge[dir][X]*out[dir][X];
  }
}



template<typename vector, typename matrix>
class dirac_staggered {
  private:
    double mass;
    field<double> staggered_eta[NDIM];

    // Note array of fields, changes with the field
    field<matrix> (&gauge)[NDIM];
  public:

    using vector_type = vector;
    parity par = ALL;

    // Constructor: initialize mass, gauge and eta
    dirac_staggered(dirac_staggered &d) : gauge(d.gauge), mass(d.mass) {
      // Initialize the eta field (Share this?)
      init_staggered_eta(staggered_eta);
    }
    dirac_staggered(double m, field<matrix> (&g)[NDIM]) : gauge(g), mass(m)  {
      init_staggered_eta(staggered_eta);
    }

    // Applies the operator to in
    void apply( const field<vector> & in, field<vector> & out){
      out[ALL] = 0;
      dirac_staggered_diag(mass, in, out, ALL);
      dirac_staggered_hop(gauge, in, out, staggered_eta, ALL, 1);
    }

    // Applies the conjugate of the operator
    void dagger( const field<vector> & in, field<vector> & out){
      out[ALL] = 0;
      dirac_staggered_diag(mass, in, out, ALL);
      dirac_staggered_hop(gauge, in, out, staggered_eta, ALL, -1);
    }

    // Applies the derivative of the Dirac operator with respect
    // to the gauge field
    void force( const field<vector> & chi, const field<vector> & psi, field<matrix> (&force)[NDIM], int sign=1){
      dirac_staggered_calc_force(gauge, chi, psi, force, staggered_eta, sign, ALL);
    }
};


// Multiplying from the left applies the standard Dirac operator
template<typename vector, typename matrix>
field<vector> operator* (dirac_staggered<field<vector>, field<matrix>> D, const field<vector> & in) {
  field<vector> out;
  out.copy_boundary_condition(in);
  D.apply(in, out);
  return out;
}

// Multiplying from the right applies the conjugate
template<typename vector, typename matrix>
field<vector> operator* (const field<vector> & in, dirac_staggered<field<vector>, field<matrix>> D) {
  field<vector> out;
  out.copy_boundary_condition(in);
  D.dagger(in, out);
  return out;
}




template<typename vector, typename matrix>
class dirac_staggered_evenodd {
  private:
    double mass;
    field<double> staggered_eta[NDIM];

    field<matrix> (&gauge)[NDIM];
  public:

    using vector_type = vector;
    using matrix_type = matrix;

    parity par = EVEN;

    dirac_staggered_evenodd(dirac_staggered_evenodd &d) : gauge(d.gauge), mass(d.mass) {
      init_staggered_eta(staggered_eta);
    }
    dirac_staggered_evenodd(double m, field<matrix> (&U)[NDIM]) : gauge(U), mass(m) {
      init_staggered_eta(staggered_eta);
    }


    // Applies the operator to in
    inline void apply(  field<vector_type> & in, field<vector_type> & out){
      out[ALL] = 0;
      dirac_staggered_diag(mass, in, out, EVEN);

      dirac_staggered_hop(gauge, in, out, staggered_eta, ODD, 1);
      dirac_staggered_diag_inverse(mass, out, ODD);
      dirac_staggered_hop(gauge, out, out, staggered_eta, EVEN, 1);
    }

    // Applies the conjugate of the operator
    inline void dagger( field<vector_type> & in, field<vector_type> & out){
      out[ALL] = 0;
      dirac_staggered_diag(mass, in, out, EVEN);

      dirac_staggered_hop(gauge, in, out, staggered_eta, ODD, -1);
      dirac_staggered_diag_inverse(mass, out, ODD);
      dirac_staggered_hop(gauge, out, out, staggered_eta, EVEN, -1);
    }

    // Applies the derivative of the Dirac operator with respect
    // to the gauge field
    inline void force(const field<vector_type> & chi, const field<vector_type> & psi, field<matrix_type> (&force)[NDIM], int sign){
      field<matrix_type> force2[NDIM];
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