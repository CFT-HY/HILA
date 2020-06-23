#ifndef __DIRAC_WILSON_H__
#define __DIRAC_WILSON_H__

#include "../plumbing/defs.h"
#include "../datatypes/cmplx.h"
#include "../datatypes/matrix.h"
#include "../datatypes/wilson_vector.h"
#include "../plumbing/field.h"
#include "../hmc/gauge_field.h"


template<int N, typename radix>
field<half_Wilson_vector<N, radix>> wilson_dirac_temp_vector[NDIM];



template<int N, typename radix, typename matrix>
inline void Dirac_Wilson_hop(
  const field<matrix> *gauge, const double kappa,
  const field<Wilson_vector<N, radix>> &v_in,
  field<Wilson_vector<N, radix>> &v_out,
  parity par, int sign)
{
  field<half_Wilson_vector<N, radix>> (&vtemp)[NDIM] = wilson_dirac_temp_vector<N, radix>;
  foralldir(dir)
    vtemp[dir].copy_boundary_condition(v_in);

  // Run neighbour fetches and multiplications
  foralldir(dir){
    direction odir = opp_dir( (direction)dir );
    // First multiply the by conjugate before communicating
    onsites(opp_parity(par)){
      half_Wilson_vector<N, radix> h(v_in[X], dir, -sign);
      vtemp[dir][X] = gauge[dir][X].conjugate()*h;
    }

    vtemp[dir].set_boundary_condition(dir, v_in.get_boundary_condition(dir));
    //vtemp[dir].start_get(odir);
  }
  foralldir(dir){
    direction odir = opp_dir( (direction)dir );
    onsites(par){
      half_Wilson_vector<N, radix> h1(v_in[X+dir], dir, sign);
      v_out[X] = v_out[X] 
               - (kappa*gauge[dir][X]*h1).expand(dir, sign)
               - (kappa*vtemp[dir][X-dir]).expand(dir, -sign);
    }
  }
}



template<int N, typename radix>
inline void Dirac_Wilson_diag(
  const field<Wilson_vector<N, radix>> &v_in,
  field<Wilson_vector<N, radix>> &v_out,
  parity par)
{
  v_out[par] = v_in[X];
}


template<int N, typename radix>
inline void Dirac_Wilson_diag_inverse(
  field<Wilson_vector<N, radix>> &v,
  parity par)
{}



template<int N, typename radix, typename gaugetype, typename momtype>
inline void Dirac_Wilson_calc_force(
  const field<gaugetype> *gauge,
  const double kappa,
  const field<Wilson_vector<N, radix>> &chi,
  const field<Wilson_vector<N, radix>> &psi,
  field<momtype> (&out)[NDIM],
  parity par,
  int sign)
{
  field<half_Wilson_vector<N, radix>> (&vtemp)[NDIM] = wilson_dirac_temp_vector<N, radix>;
  vtemp[0].copy_boundary_condition(chi);
  vtemp[1].copy_boundary_condition(chi);
  
  foralldir(dir){
    onsites(opp_parity(par)){
      half_Wilson_vector<N, radix> hw(chi[X], dir, -sign);
      vtemp[0][X] = hw;
    }
    onsites(par){
      half_Wilson_vector<N, radix> hw(psi[X], dir, sign);
      vtemp[1][X] = hw;
    }

    out[dir][ALL] = 0;
    out[dir][par] = -kappa * (
      ( vtemp[0][X+dir].expand(dir,-sign) ).outer_product(psi[X])
    );
    out[dir][opp_parity(par)] = out[dir][X] - kappa * (
      ( vtemp[1][X+dir].expand(dir, sign) ).outer_product(chi[X])
    );

    out[dir][ALL] = gauge[dir][X]*out[dir][X];
  }
}



template<int N, typename radix, typename matrix>
class Dirac_Wilson {
  private:
    double kappa;

    // Note array of fields, changes with the field
    field<matrix> (&gauge)[NDIM];
  
  public:

    using vector_type = Wilson_vector<N, radix>;
    using matrix_type = matrix;

    parity par = ALL;

    // Constructor: initialize mass and gauge
    Dirac_Wilson(Dirac_Wilson &d) : gauge(d.gauge), kappa(d.kappa) {}
    Dirac_Wilson(double k, field<matrix> (&U)[NDIM]) : gauge(U), kappa(k) {}
    Dirac_Wilson(double k, gauge_field<matrix> &g) : gauge(g.gauge), kappa(k) {}
    Dirac_Wilson(double k, represented_gauge_field<matrix> &g) : gauge(g.gauge), kappa(k) {}

    // Applies the operator to in
    inline void apply( const field<vector_type> & in, field<vector_type> & out){
      Dirac_Wilson_diag(in, out, ALL);
      Dirac_Wilson_hop(gauge, kappa, in, out, ALL, 1);
    }

    // Applies the conjugate of the operator
    inline void dagger( const field<vector_type> & in, field<vector_type> & out){
      Dirac_Wilson_diag(in, out, ALL);
      Dirac_Wilson_hop(gauge, kappa, in, out, ALL, -1);
    }

    // Applies the derivative of the Dirac operator with respect
    // to the gauge field
    template<typename momtype>
    inline void force(const field<vector_type> & chi,  const field<vector_type> & psi, field<momtype> (&force)[NDIM], int sign=1){
      Dirac_Wilson_calc_force(gauge, kappa, chi, psi, force, ALL, sign);
    }
};


// Multiplying from the left applies the standard Dirac operator
template<int N, typename radix, typename matrix>
field<Wilson_vector<N, radix>> operator* (Dirac_Wilson<N, radix, matrix> D, const field<Wilson_vector<N, radix>> & in) {
  field<Wilson_vector<N, radix>> out;
  D.apply(in, out);
  return out;
}

// Multiplying from the right applies the conjugate
template<int N, typename radix, typename matrix>
field<Wilson_vector<N, radix>> operator* (const field<Wilson_vector<N, radix>> & in, Dirac_Wilson<N, radix, matrix> D) {
  field<Wilson_vector<N, radix>> out;
  D.dagger(in, out);
  return out;
}







template<int N, typename radix, typename matrix>
class Dirac_Wilson_evenodd {
  private:
    double kappa;
    field<matrix> (&gauge)[NDIM];
  public:

    using vector_type = Wilson_vector<N, radix>;
    using matrix_type = matrix;

    // The parity 
    parity par = EVEN;

    Dirac_Wilson_evenodd(Dirac_Wilson_evenodd &d) : gauge(d.gauge), kappa(d.kappa) {}
    Dirac_Wilson_evenodd(double k, field<matrix> (&U)[NDIM]) : gauge(U), kappa(k) {}
    Dirac_Wilson_evenodd(double k, gauge_field<matrix> &g) : gauge(g.gauge), kappa(k) {}
    Dirac_Wilson_evenodd(double k, represented_gauge_field<matrix> &g) : gauge(g.gauge), kappa(k) {}


    // Applies the operator to in
    inline void apply( field<vector_type> & in, field<vector_type> & out){
      out[ALL] = 0;
      Dirac_Wilson_diag(in, out, EVEN);

      Dirac_Wilson_hop(gauge, kappa, in, out, ODD, 1);
      Dirac_Wilson_diag_inverse(out, ODD);
      Dirac_Wilson_hop(gauge, -kappa, out, out, EVEN, 1);
      out[ODD] = 0;
    }

    // Applies the conjugate of the operator
    inline void dagger( field<vector_type> & in, field<vector_type> & out){
      out[ALL] = 0;
      Dirac_Wilson_diag(in, out, EVEN);

      Dirac_Wilson_hop(gauge, kappa, in, out, ODD, -1);
      Dirac_Wilson_diag_inverse(out, ODD);
      Dirac_Wilson_hop(gauge, -kappa, out, out, EVEN, -1);
      out[ODD] = 0;
    }

    // Applies the derivative of the Dirac operator with respect
    // to the gauge field
    template<typename momtype>
    inline void force(const field<vector_type> & chi, const field<vector_type> & psi, field<momtype> (&force)[NDIM], int sign){
      field<momtype> force2[NDIM];
      field<vector_type> tmp;
      tmp.copy_boundary_condition(chi);

      tmp[ALL] = 0;
      Dirac_Wilson_hop(gauge, kappa, chi, tmp, ODD, -sign);
      Dirac_Wilson_diag_inverse(tmp, ODD);
      Dirac_Wilson_calc_force(gauge, -kappa, tmp, psi, force, EVEN, sign);

      tmp[ALL] = 0;
      Dirac_Wilson_hop(gauge, kappa, psi, tmp, ODD, sign);
      Dirac_Wilson_diag_inverse(tmp, ODD);
      Dirac_Wilson_calc_force(gauge, -kappa, chi, tmp, force2, ODD, sign);

      foralldir(dir)
        force[dir][ALL] = force[dir][X] + force2[dir][X];
    }
};


#endif