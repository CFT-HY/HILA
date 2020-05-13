#ifndef __DIRAC_STAGGERED_NJL_H__
#define __DIRAC_STAGGERED_NJL_H__

#include "../plumbing/defs.h"
#include "../datatypes/cmplx.h"
#include "../datatypes/general_matrix.h"
#include "../datatypes/wilson_vector.h"
#include "../plumbing/field.h"

#include "../plumbing/dirac.h"




template<int N>
class dirac_staggered_NJL {
  private:
    double mass;
    field<SU_vector<N>> vtemp[NDIM];
    field<double> staggered_eta[NDIM];

    // Note array of fields, changes with the field
    field<SU<N>> *gauge;
    field<double> &sigma;
    field<double> &pi;
  
  public:

    using vector_type = SU_vector<N>;

    // Constructor: initialize mass, gauge and eta
    dirac_staggered_NJL(dirac_staggered_NJL &d) : sigma(d.sigma), pi(d.pi) {
      mass = d.mass;
      gauge = d.gauge;

      // Initialize the eta field (Share this?)
      init_staggered_eta(staggered_eta);
    }
  
    // Constructor: initialize mass, gauge and eta
    dirac_staggered_NJL(double m, field<SU<N>> *U, field<double> &s, field<double> &p) : sigma(s), pi(p) {
      // Set mass and gauge field
      mass = m;
      gauge = (field<SU<N>>*) U;

      // Initialize the eta field
      init_staggered_eta(staggered_eta);
    }

    // Update mass
    void set_mass(double m){
      mass = m;
    }



    // Applies the operator to in
    void apply( const field<SU_vector<N>> & in, field<SU_vector<N>> & out){
      dirac_staggered_apply(gauge, mass, in, out, staggered_eta, vtemp);
    }

    // Applies the conjugate of the operator
    void dagger( const field<SU_vector<N>> & in, field<SU_vector<N>> & out){
      dirac_staggered_dagger(gauge, mass, in, out, staggered_eta, vtemp);
    }

    // Applies the derivative of the Dirac operator with respect
    // to the gauge field
    void force( const field<SU_vector<N>> & psi, const field<SU_vector<N>> & chi, field<SU<N>> (&force)[NDIM]){
      dirac_staggered_calc_force(gauge, mass, psi, chi, force, staggered_eta, vtemp);
    }
};




#endif