#ifndef __DIRAC_WILSON_H__
#define __DIRAC_WILSON_H__

#include "../../plumbing/defs.h"
#include "../../datatypes/cmplx.h"
#include "../../datatypes/general_matrix.h"
#include "../../datatypes/wilson_vector.h"
#include "../../plumbing/field.h"




template<int n, typename radix=double>
void dirac_wilson_apply(
  const field<SU<n, radix>> *gauge,
  const double kappa,
  const field<Wilson_vector<n, radix>> &v_in,
  field<Wilson_vector<n, radix>> &v_out,
  field<half_Wilson_vector<n, radix>> (&vtemp)[NDIM])
{
  // Start getting neighbours
  foralldir(dir){
    v_in.start_get(dir);
  }

  // The diagonal term
  v_out[ALL] = v_in[X];

  // Run neighbour fetches and multiplications
  foralldir(dir){
    direction odir = opp_dir( (direction)dir );
    // First multiply the by conjugate before communicating
    onsites(ALL){
      half_Wilson_vector<n, radix> h(v_in[X], odir);
      vtemp[dir][X] = gauge[dir][X].conjugate()*h;
    }
    vtemp[dir].start_get(odir);
  }
  foralldir(dir){
    direction odir = opp_dir( (direction)dir );
    onsites(ALL){
      half_Wilson_vector<n, radix> h1(v_in[X+dir], dir);
      v_out[X] = v_out[X] 
               - (kappa*gauge[dir][X]*h1).expand(dir)
               - (kappa*vtemp[dir][X-dir]).expand(odir);
    }
  }
}


template<int n, typename radix=double>
void dirac_wilson_dagger(
  const field<SU<n, radix>> *gauge,
  const double kappa,
  const field<Wilson_vector<n, radix>> &v_in,
  field<Wilson_vector<n, radix>> &v_out,
  field<half_Wilson_vector<n, radix>> (&vtemp)[NDIM])
{
  // Start getting neighbours
  foralldir(dir){
    v_in.start_get(dir);
  }

  // The diagonal term
  v_out[ALL] = v_in[X];

  // Run neighbour fetches and multiplications
  foralldir(dir){
    direction odir = opp_dir( (direction)dir );
    // First multiply the by conjugate before communicating
    onsites(ALL){
      half_Wilson_vector<n, radix> h(v_in[X], dir);
      vtemp[dir][X] = gauge[dir][X].conjugate() * h;
    }
    vtemp[dir].start_get(odir);
  }
  foralldir(dir){
    direction odir = opp_dir( (direction)dir );
    onsites(ALL){
      half_Wilson_vector<n, radix> h1(v_in[X+dir], odir);
      v_out[X] = v_out[X] 
               - (kappa*gauge[dir][X]*h1).expand(odir)
               - (kappa*vtemp[dir][X-dir]).expand(dir);
    }
  }
}




template<int n, typename radix=double>
void dirac_wilson_calc_force(
  const field<SU<n, radix>> *gauge,
  const double kappa,
  const field<Wilson_vector<n, radix>> &psi,
  const field<Wilson_vector<n, radix>> &chi,
  field<SU<n, radix>> (&out)[NDIM],
  field<half_Wilson_vector<n, radix>> (&vtemp)[NDIM])
{
  foralldir(dir){
    direction odir = opp_dir( (direction)dir );
    onsites(ALL){
      vtemp[0][X] = half_Wilson_vector<n, radix>(psi[X], odir);
      vtemp[1][X] = half_Wilson_vector<n, radix>(chi[X], dir);
    }
    
    out[dir][ALL] = -kappa * (
        ( vtemp[0][X+dir].expand(odir) ).outer_product(chi[X])
      + ( vtemp[1][X+dir].expand(dir)  ).outer_product(psi[X])
    );
    
    out[dir][ALL] = gauge[dir][X]*out[dir][X];
  }
}



template<int n, typename radix=double>
class dirac_wilson {
  private:
    double kappa;
    field<half_Wilson_vector<n, radix>> vtemp[NDIM];

    // Note array of fields, changes with the field
    field<SU<n, radix>> *gauge;
  
  public:

    using vector_type = Wilson_vector<n, radix>;

    // Constructor: initialize mass, gauge and eta
    dirac_wilson(dirac_wilson &d) {
      kappa = d.kappa;
      gauge = d.gauge;
    }
  
    // Constructor: initialize mass, gauge and eta
    dirac_wilson(double k, field<SU<n, radix>> *U) {
      // Set mass and gauge field
      kappa = k;
      gauge = (field<SU<n, radix>>*) U;
    }


    // Applies the operator to in
    void apply( const field<Wilson_vector<n, radix>> & in, field<Wilson_vector<n, radix>> & out){
      dirac_wilson_apply(gauge, kappa, in, out, vtemp);
    }

    // Applies the conjugate of the operator
    void dagger( const field<Wilson_vector<n, radix>> & in, field<Wilson_vector<n, radix>> & out){
      dirac_wilson_dagger(gauge, kappa, in, out, vtemp);
    }

    // Applies the derivative of the Dirac operator with respect
    // to the gauge field
    void force( const field<Wilson_vector<n, radix>> & psi, const field<Wilson_vector<n, radix>> & chi, field<SU<n, radix>> (&force)[NDIM]){
      dirac_wilson_calc_force(gauge, kappa, psi, chi, force, vtemp);
    }
};


// Multiplying from the left applies the standard Dirac operator
template<int n, typename radix>
field<Wilson_vector<n, radix>> operator* (dirac_wilson<n, radix> D, const field<Wilson_vector<n, radix>> & in) {
  field<Wilson_vector<n, radix>> out;
  D.apply(in, out);
  return out;
}

// Multiplying from the right applies the conjugate
template<int n, typename radix>
field<Wilson_vector<n, radix>> operator* (const field<Wilson_vector<n, radix>> & in, dirac_wilson<n, radix> D) {
  field<Wilson_vector<n, radix>> out;
  D.dagger(in, out);
  return out;
}


#endif