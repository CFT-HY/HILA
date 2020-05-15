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
void dirac_wilson_diag_set(
  const field<Wilson_vector<vector>> &v_in,
  field<Wilson_vector<vector>> &v_out,
  parity par)
{
  v_out[par] = v_in[X];
}





template<typename vector, typename matrix>
void dirac_wilson_calc_force(
  const field<matrix> *gauge,
  const double kappa,
  const field<Wilson_vector<vector>> &psi,
  const field<Wilson_vector<vector>> &chi,
  field<matrix> (&out)[NDIM],
  int sign)
{
  field<half_Wilson_vector<vector>> (&vtemp)[NDIM] = wilson_dirac_temp_vector<vector>;

  foralldir(dir){
    direction odir = opp_dir( (direction)dir );
    onsites(ALL){
      vtemp[0][X] = half_Wilson_vector<vector>(psi[X], dir, -sign);
      vtemp[1][X] = half_Wilson_vector<vector>(chi[X], dir, sign);
    }

    out[dir][ALL] = -kappa * (
        ( vtemp[0][X+dir].expand(dir,-sign) ).outer_product(chi[X])
      + ( vtemp[1][X+dir].expand(dir, sign) ).outer_product(psi[X])
    );
    out[dir][ALL] = gauge[dir][X]*out[dir][X];
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

    // Constructor: initialize mass, gauge and eta
    dirac_wilson(dirac_wilson &d) : gauge(d.gauge), kappa(d.kappa) {}
  
    // Constructor: initialize mass, gauge and eta
    dirac_wilson(double k, field<matrix> (&U)[NDIM]) : gauge(U), kappa(k) {}

    // Applies the operator to in
    void apply( const field<Wilson_vector<vector>> & in, field<Wilson_vector<vector>> & out){
      dirac_wilson_diag_set(in, out, ALL);
      dirac_wilson_hop(gauge, kappa, in, out, ALL, 1);
    }

    // Applies the conjugate of the operator
    void dagger( const field<Wilson_vector<vector>> & in, field<Wilson_vector<vector>> & out){
      dirac_wilson_diag_set(in, out, ALL);
      dirac_wilson_hop(gauge, kappa, in, out, ALL, -1);
    }

    // Applies the derivative of the Dirac operator with respect
    // to the gauge field
    void force( const field<Wilson_vector<vector>> & psi, const field<Wilson_vector<vector>> & chi, field<matrix> (&force)[NDIM], int sign=1){
      dirac_wilson_calc_force(gauge, kappa, psi, chi, force, sign);
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


#endif