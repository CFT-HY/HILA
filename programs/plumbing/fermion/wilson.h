#ifndef __DIRAC_WILSON_H__
#define __DIRAC_WILSON_H__

#include "../../plumbing/defs.h"
#include "../../datatypes/cmplx.h"
#include "../../datatypes/general_matrix.h"
#include "../../datatypes/wilson_vector.h"
#include "../../plumbing/field.h"




template<typename vector, typename matrix>
void dirac_wilson_apply(
  const field<matrix> *gauge,
  const double kappa,
  const field<Wilson_vector<vector>> &v_in,
  field<Wilson_vector<vector>> &v_out,
  field<half_Wilson_vector<vector>> (&vtemp)[NDIM])
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
      half_Wilson_vector<vector> h(v_in[X], odir);
      vtemp[dir][X] = gauge[dir][X].conjugate()*h;
    }
    vtemp[dir].start_get(odir);
  }
  foralldir(dir){
    direction odir = opp_dir( (direction)dir );
    onsites(ALL){
      half_Wilson_vector<vector> h1(v_in[X+dir], dir);
      v_out[X] = v_out[X] 
               - (kappa*gauge[dir][X]*h1).expand(dir)
               - (kappa*vtemp[dir][X-dir]).expand(odir);
    }
  }
}


template<typename vector, typename matrix>
void dirac_wilson_dagger(
  const field<matrix> *gauge,
  const double kappa,
  const field<Wilson_vector<vector>> &v_in,
  field<Wilson_vector<vector>> &v_out,
  field<half_Wilson_vector<vector>> (&vtemp)[NDIM])
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
      half_Wilson_vector<vector> h(v_in[X], dir);
      vtemp[dir][X] = gauge[dir][X].conjugate() * h;
    }
    vtemp[dir].start_get(odir);
  }
  foralldir(dir){
    direction odir = opp_dir( (direction)dir );
    onsites(ALL){
      half_Wilson_vector<vector> h1(v_in[X+dir], odir);
      v_out[X] = v_out[X] 
               - (kappa*gauge[dir][X]*h1).expand(odir)
               - (kappa*vtemp[dir][X-dir]).expand(dir);
    }
  }
}




template<typename vector, typename matrix>
void dirac_wilson_calc_force(
  const field<matrix> *gauge,
  const double kappa,
  const field<Wilson_vector<vector>> &psi,
  const field<Wilson_vector<vector>> &chi,
  field<matrix> (&out)[NDIM],
  field<half_Wilson_vector<vector>> (&vtemp)[NDIM],
  int sign)
{
  foralldir(dir){
    if(sign == 1){
      direction odir = opp_dir( (direction)dir );
      onsites(ALL){
        vtemp[0][X] = half_Wilson_vector<vector>(psi[X], odir);
        vtemp[1][X] = half_Wilson_vector<vector>(chi[X], dir);
      }

      out[dir][ALL] = -kappa * (
          ( vtemp[0][X+dir].expand(odir) ).outer_product(chi[X])
        + ( vtemp[1][X+dir].expand(dir)  ).outer_product(psi[X])
      );
    } else {
      direction odir = opp_dir( (direction)dir );
      onsites(ALL){
        vtemp[0][X] = half_Wilson_vector<vector>(psi[X], dir);
        vtemp[1][X] = half_Wilson_vector<vector>(chi[X], odir);
      }

      out[dir][ALL] = -kappa * (
          ( vtemp[0][X+dir].expand(dir) ).outer_product(chi[X])
        + ( vtemp[1][X+dir].expand(odir)  ).outer_product(psi[X])
      );
    }
    out[dir][ALL] = gauge[dir][X]*out[dir][X];
  }
}


template<typename vector, typename matrix>
class dirac_wilson {
  private:
    double kappa;
    field<half_Wilson_vector<vector>> vtemp[NDIM];

    // Note array of fields, changes with the field
    field<matrix> *gauge;
  
  public:

    using vector_type = Wilson_vector<vector>;

    // Constructor: initialize mass, gauge and eta
    dirac_wilson(dirac_wilson &d) {
      kappa = d.kappa;
      gauge = d.gauge;
    }
  
    // Constructor: initialize mass, gauge and eta
    dirac_wilson(double k, field<matrix> *U) {
      // Set mass and gauge field
      kappa = k;
      gauge = (field<matrix>*) U;
    }


    // Applies the operator to in
    void apply( const field<Wilson_vector<vector>> & in, field<Wilson_vector<vector>> & out){
      dirac_wilson_apply(gauge, kappa, in, out, vtemp);
    }

    // Applies the conjugate of the operator
    void dagger( const field<Wilson_vector<vector>> & in, field<Wilson_vector<vector>> & out){
      dirac_wilson_dagger(gauge, kappa, in, out, vtemp);
    }

    // Applies the derivative of the Dirac operator with respect
    // to the gauge field
    void force( const field<Wilson_vector<vector>> & psi, const field<Wilson_vector<vector>> & chi, field<matrix> (&force)[NDIM], int sign=1){
      dirac_wilson_calc_force(gauge, kappa, psi, chi, force, vtemp, sign);
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