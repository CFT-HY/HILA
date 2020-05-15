#ifndef __DIRAC_WILSON_H__
#define __DIRAC_WILSON_H__

#include "../../plumbing/defs.h"
#include "../../datatypes/cmplx.h"
#include "../../datatypes/general_matrix.h"
#include "../../datatypes/wilson_vector.h"
#include "../../plumbing/field.h"



template<typename vector>
field<half_Wilson_vector<vector>> wilson_dirac_temp_vector[NDIM];



template<typename vector, typename matrix>
void dirac_wilson_hop(
  const field<matrix> *gauge, const double kappa,
  const field<Wilson_vector<vector>> &v_in,
  field<Wilson_vector<vector>> &v_out,
  parity par, int sign)
{
  field<half_Wilson_vector<vector>> (&vtemp)[NDIM] = wilson_dirac_temp_vector<vector>;
  
  // Run neighbour fetches and multiplications
  foralldir(dir){
    direction odir = opp_dir( (direction)dir );
    // First multiply the by conjugate before communicating
    onsites(opp_parity(par)){
      half_Wilson_vector<vector> h(v_in[X], dir, -sign);
      vtemp[dir][X] = gauge[dir][X].conjugate()*h;
    }
    vtemp[dir].start_get(odir);
  }
  foralldir(dir){
    direction odir = opp_dir( (direction)dir );
    onsites(par){
      half_Wilson_vector<vector> h1(v_in[X+dir], dir, sign);
      v_out[X] = v_out[X] 
               - (kappa*gauge[dir][X]*h1).expand(dir, sign)
               - (kappa*vtemp[dir][X-dir]).expand(dir, -sign);
    }
  }
}



template<typename vector>
void dirac_wilson_diag(
  const field<Wilson_vector<vector>> &v_in,
  field<Wilson_vector<vector>> &v_out,
  parity par)
{
  v_out[par] += v_in[X];
}


template<typename vector>
void dirac_wilson_diag_inverse(
  field<Wilson_vector<vector>> &v,
  parity par)
{}



template<typename vector, typename matrix>
void dirac_wilson_calc_force(
  const field<matrix> *gauge,
  const double kappa,
  const field<Wilson_vector<vector>> &psi,
  const field<Wilson_vector<vector>> &chi,
  field<matrix> (&out)[NDIM],
  parity par,
  int sign)
{
  field<half_Wilson_vector<vector>> (&vtemp)[NDIM] = wilson_dirac_temp_vector<vector>;

  foralldir(dir){
    direction odir = opp_dir( (direction)dir );
    onsites(par){
      vtemp[0][X] = half_Wilson_vector<vector>(psi[X], dir, -sign);
      vtemp[1][X] = half_Wilson_vector<vector>(chi[X], dir, sign);
    }

    out[dir][par] = -kappa * (
        ( vtemp[0][X+dir].expand(dir,-sign) ).outer_product(chi[X])
      + ( vtemp[1][X+dir].expand(dir, sign) ).outer_product(psi[X])
    );
    out[dir][par] = gauge[dir][X]*out[dir][X];
  }
}



template<typename vector, typename matrix>
class dirac_wilson {
  private:
    double kappa;

    // Note array of fields, changes with the field
    field<matrix> (&gauge)[NDIM];
  
  public:

    using vector_type = Wilson_vector<vector>;
    using matrix_type = matrix;

    // Constructor: initialize mass and gauge
    dirac_wilson(dirac_wilson &d) : gauge(d.gauge), kappa(d.kappa) {}
  
    // Constructor: initialize mass and gauge
    dirac_wilson(double k, field<matrix> (&U)[NDIM]) : gauge(U), kappa(k) {}


    void diag(const field<vector_type> & in, field<vector_type> & out, parity par){
      dirac_wilson_diag(in, out, par);
    }

    void diag_inverse(field<vector_type> & vec, parity par){
      dirac_wilson_diag_inverse(vec, par);
    }

    void hop(const field<vector_type> & in, field<vector_type> & out, parity par, int sign){
      dirac_wilson_hop(gauge, kappa, in, out, par, sign);
    }

    void hop_deriv( const field<vector_type> & psi, const field<vector_type> & chi, field<matrix> (&force)[NDIM], parity par, int sign){
      dirac_wilson_calc_force(gauge, kappa, psi, chi, force, par, sign);
    }

    // Applies the operator to in
    void apply( const field<vector_type> & in, field<vector_type> & out){
      out[ALL] = 0;
      dirac_wilson_diag(in, out, ALL);
      dirac_wilson_hop(gauge, kappa, in, out, ALL, 1);
    }

    // Applies the conjugate of the operator
    void dagger( const field<vector_type> & in, field<vector_type> & out){
      out[ALL] = 0;
      dirac_wilson_diag(in, out, ALL);
      dirac_wilson_hop(gauge, kappa, in, out, ALL, -1);
    }

    // Applies the derivative of the Dirac operator with respect
    // to the gauge field
    void force( const field<vector_type> & psi, const field<vector_type> & chi, field<matrix> (&force)[NDIM], int sign=1){
      dirac_wilson_calc_force(gauge, kappa, psi, chi, force, ALL, sign);
    }
};


// Multiplying from the left applies the standard Dirac operator
template<typename vector, typename matrix>
field<Wilson_vector<vector>> operator* (dirac_wilson<vector, matrix> D, const field<Wilson_vector<vector>> & in) {
  field<Wilson_vector<vector>> out;
  D.apply(in, out);
  return out;
}

// Multiplying from the right applies the conjugate
template<typename vector, typename matrix>
field<Wilson_vector<vector>> operator* (const field<Wilson_vector<vector>> & in, dirac_wilson<vector, matrix> D) {
  field<Wilson_vector<vector>> out;
  D.dagger(in, out);
  return out;
}





template<typename dirac_op>
class precondition_evenodd {
  private:
    dirac_op D;
  
  public:

    using vector_type = typename dirac_op::vector_type;
    using matrix_type = typename dirac_op::matrix_type;

    precondition_evenodd(precondition_evenodd &d) : D(d.D) {}
    precondition_evenodd(dirac_op d) : D(d) {}

    // Applies the operator to in
    void apply( const field<vector_type> & in, field<vector_type> & out){
      out[ALL] = 0;
      D.diag(in, out, EVEN);

      D.hop(in, out, ODD, 1);
      D.diag_inverse(out, ODD);
      D.hop(out, out, EVEN, 1);
      out[ODD] = 0;
    }

    // Applies the conjugate of the operator
    void dagger( const field<vector_type> & in, field<vector_type> & out){
      out[ALL] = 0;
      D.diag(in, out, EVEN);

      D.hop(in, out, ODD, -1);
      D.diag_inverse(out, ODD);
      D.hop(out, out, EVEN, -1);
      out[ODD] = 0;
    }

    // Applies the derivative of the Dirac operator with respect
    // to the gauge field
    void force( const field<vector_type> & chi, const field<vector_type> & psi, field<matrix_type> (&force)[NDIM], int sign){
      field<vector_type> tmp;
      tmp[ALL] = 0;
      D.hop(psi, tmp, ODD, sign);
      D.diag_inverse(tmp, ODD);
      D.hop_deriv(tmp, chi, force, ALL, sign);

      tmp[ALL] = 0;
      D.hop(chi, tmp, ODD, -sign);
      D.diag_inverse(tmp, ODD);
      D.hop_deriv(psi, tmp, force, ALL, sign);
    }
};


#endif