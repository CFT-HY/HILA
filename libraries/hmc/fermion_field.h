#ifndef FERMION_FIELD_H
#define FERMION_FIELD_H

#include "gauge_field.h"
#include "dirac/Hasenbusch.h"
#include "dirac/conjugate_gradient.h"
#include "MRE_guess.h"
#include <cmath>



template<typename gauge_field, typename DIRAC_OP>
class fermion_action : public action_base{
  public:
    using vector_type = typename DIRAC_OP::vector_type;
    using momtype = squarematrix<gauge_field::N, cmplx<typename gauge_field::basetype>>;
    gauge_field &gauge;
    DIRAC_OP &D;
    field<vector_type> chi;

    // We save a few previous invertions to build an initial guess.
    // old_chi contains a list of these
    int MRE_size = 0;
    std::vector<field<vector_type>> old_chi_inv;

    void setup(int mre_guess_size){
      #if NDIM > 3
      //chi.set_boundary_condition(TUP, boundary_condition_t::ANTIPERIODIC);
      //chi.set_boundary_condition(TDOWN, boundary_condition_t::ANTIPERIODIC);
      #endif
      MRE_size = mre_guess_size;
      old_chi_inv.resize(MRE_size);
      for(int i=0; i<MRE_size; i++){
        old_chi_inv[i][ALL] = 0;
      }
    }

    fermion_action(DIRAC_OP &d, gauge_field &g)
    : D(d), gauge(g){ 
      chi = 0.0;  // Allocates chi and sets it to zero
      setup(0);
    }

    fermion_action(DIRAC_OP &d, gauge_field &g, int mre_guess_size)
    : D(d), gauge(g){ 
      chi = 0.0;  // Allocates chi and sets it to zero
      setup(mre_guess_size);
    }


    fermion_action(fermion_action &fa)
    : gauge(fa.gauge), D(fa.D){
      chi = fa.chi;  // Copies the field
      setup(fa.MRE_size);
    }

      }
    }

    // Return the value of the action with the current
    // field configuration
    double action(){ 
      field<vector_type> psi;
      psi.copy_boundary_condition(chi);
      CG<DIRAC_OP> inverse(D);
      double action = 0;
    
      gauge.refresh();

      psi=0;
      inverse.apply(chi,psi);
      onsites(D.par){
        action += chi[X].rdot(psi[X]);
      }
      hila::output << "fermion action " << action << "\n";
      return action;
    }

    /// Generate a pseudofermion field with a distribution given
    /// by the action chi 1/(D_dagger D) chi
    void draw_gaussian_fields(){
      field<vector_type> psi;
      psi.copy_boundary_condition(chi);
      gauge.refresh();

      onsites(D.par){
        if(disable_avx[X]==0){};
        psi[X].gaussian();
      }
      D.dagger(psi,chi);
    }


    // Build an initial guess for the fermion matrix inversion
    // by inverting first in the limited space of a few previous
    // solutions. These are saved in old_chi.
    void initial_guess(field<vector_type> & chi, field<vector_type> & psi){
      if(MRE_size > 0){
        MRE_guess(psi, chi, D, old_chi_inv);
      } else {
        psi[ALL]=0;
      }
      // If the gauge type is double precision, solve first in single precision
      if constexpr (std::is_same<double, typename gauge_field::basetype>::value){
        output0 << "Starting with single precision inversion\n";

        auto single_precision = gauge.get_single_precision();
        typename DIRAC_OP::type_flt D_flt(D, single_precision);
        field<typename DIRAC_OP::type_flt::vector_type> c, p, t1, t2;
        c[ALL] = chi[X];
        p[ALL] = psi[X];
        CG inverse(D_flt, 1e-6);
        inverse.apply(c, p);

        D_flt.apply(p, t1);
        D_flt.dagger(t1, t2);
        psi[ALL] = p[X];
      }
    }

    // Add new solution to the list
    void save_new_solution(field<vector_type> & psi){
      if(MRE_size > 0){
        for(int i=1; i<MRE_size; i++){
          old_chi_inv[i] = old_chi_inv[i-1];
        }
        old_chi_inv[0] = psi;
      }
    }


    // Update the momentum with the derivative of the fermion
    // action
    void force_step(double eps){
      field<vector_type> psi, Mpsi;
      psi.copy_boundary_condition(chi);
      Mpsi.copy_boundary_condition(chi);
      field<momtype> force[NDIM], force2[NDIM];

      CG<DIRAC_OP> inverse(D);
      gauge.refresh();

      output0 << "base force\n";
      initial_guess(chi, psi);
      inverse.apply(chi, psi);
      save_new_solution(psi);

      D.apply(psi, Mpsi);

      D.force(Mpsi, psi, force, 1);
      D.force(psi, Mpsi, force2, -1);

      foralldir(dir){
        force[dir][ALL] = -eps*(force[dir][X] + force2[dir][X]);
      }
      gauge.add_momentum(force);
    }
};











/* The Hasenbusch method for updating fermion fields:
 * Split the Dirac determinant into two parts, 
 * D_h1 = D + mh and
 * D_h2 = D * 1 / (D + mh)^dagger
 */


/* The first action term, with D_h1 = D + mh.
 * Since the only real difference here is an addition
 * to the original operator, we can use fermion_action
 */
template<typename gauge_field, typename DIRAC_OP>
class Hasenbusch_action_1 : public action_base {
  public:
    Hasenbusch_operator<DIRAC_OP> D_h;
    fermion_action<gauge_field, Hasenbusch_operator<DIRAC_OP>> base_action;
    double _mh;

    Hasenbusch_action_1(DIRAC_OP &d, gauge_field &g, double mh)
     : _mh(mh), D_h(d, mh), base_action(D_h, g) {}
    Hasenbusch_action_1(DIRAC_OP &d, gauge_field &g, double mh, int mre_guess_size)
     : _mh(mh), D_h(d, mh), base_action(D_h, g, mre_guess_size) {}

    Hasenbusch_action_1(Hasenbusch_action_1 &fa) : _mh(fa._mh), D_h(fa.D_h), base_action(fa.base_action){}

    double action(){return(base_action.action());}
    void draw_gaussian_fields(){base_action.draw_gaussian_fields();}
    void force_step(double eps){base_action.force_step(eps);}

};


/* The second Hasenbusch action term, D_h2 = D/(D^dagger + mh)
 */
template<typename gauge_field, typename DIRAC_OP>
class Hasenbusch_action_2 : public action_base {
  public:
    using vector_type = typename DIRAC_OP::vector_type;
    using momtype = squarematrix<gauge_field::N, cmplx<typename gauge_field::basetype>>;
    gauge_field &gauge;
    DIRAC_OP D;
    Hasenbusch_operator<DIRAC_OP> D_h;
    double mh;
    field<vector_type> chi;
    
    // We save a few previous invertions to build an initial guess.
    // old_chi contains a list of these
    int MRE_size = 0;
    std::vector<field<vector_type>> old_chi_inv;

    void setup(int mre_guess_size){
      #if NDIM > 3
      //chi.set_boundary_condition(TUP, boundary_condition_t::ANTIPERIODIC);
      //chi.set_boundary_condition(TDOWN, boundary_condition_t::ANTIPERIODIC);
      #endif
      MRE_size = mre_guess_size;
      old_chi_inv.resize(MRE_size);
      for(int i=0; i<MRE_size; i++){
        old_chi_inv[i][ALL] = 0;
      }
    }

    Hasenbusch_action_2(DIRAC_OP &d, gauge_field &g, double _mh) 
     : mh(_mh), D(d), D_h(d, _mh), gauge(g) {
      chi = 0.0;  // Allocates chi and sets it to zero
      setup(0);
    }
    Hasenbusch_action_2(DIRAC_OP &d, gauge_field &g, double _mh, int mre_guess_size) 
     : mh(_mh), D(d), D_h(d, _mh), gauge(g) {
      chi = 0.0;  // Allocates chi and sets it to zero
      setup(mre_guess_size);
    }

    Hasenbusch_action_2(Hasenbusch_action_2 &fa) 
     : mh(fa.mh), D(fa.D), D_h(fa.D_h), gauge(fa.gauge){
      chi = fa.chi;  // Copies the field
      setup(fa.MRE_size);
    }


    // Return the value of the action with the current
    // field configuration
    double action(){
      field<vector_type> psi;
      field<vector_type> v;
      psi.copy_boundary_condition(chi);
      v.copy_boundary_condition(chi);
      double action = 0;

      gauge.refresh();

      CG<DIRAC_OP> inverse(D);

      v[ALL] = 0;
      D_h.dagger(chi, psi);
      inverse.apply(psi, v);
      D_h.apply(v, psi);
      onsites(EVEN){
        if(disable_avx[X]==0){};
        action += chi[X].rdot(psi[X]);
      }
      
      hila::output << "h2 " << action << "\n";
      return action;
    }

    
    /// Generate a pseudofermion field with a distribution given
    /// by the action chi 1/(D_dagger D) chi
    void draw_gaussian_fields(){
      field<vector_type> psi;
      field<vector_type> v;
      psi.copy_boundary_condition(chi);
      v.copy_boundary_condition(chi);
      CG inverse_h(D_h, 1e-12); // Applies 1/(D_h^dagger D_h)
      gauge.refresh();

      onsites(D.par){
        if(disable_avx[X]==0){};
        psi[X].gaussian();
      }
      v=0;
      D_h.dagger(psi, v);
      inverse_h.apply(v, psi); // 1/D_h * psi
      D.apply(psi, chi);
    }


    // Build an initial guess for the fermion matrix inversion
    // by inverting first in the limited space of a few previous
    // solutions. These are saved in old_chi.
    void initial_guess(field<vector_type> & chi, field<vector_type> & psi){
      if(MRE_size > 0){
        MRE_guess(psi, chi, D, old_chi_inv);
      } else {
        psi[ALL]=0;
      }
      // If the gauge type is double precision, solve first in single precision
      if constexpr (std::is_same<double, typename gauge_field::basetype>::value){
        output0 << "Starting with single precision inversion\n";

        auto single_precision = gauge.get_single_precision();
        typename DIRAC_OP::type_flt D_flt(D, single_precision);
        field<typename DIRAC_OP::type_flt::vector_type> c, p, t1, t2;
        c[ALL] = chi[X];
        p[ALL] = psi[X];
        CG inverse(D_flt, 1e-6);
        inverse.apply(c, p);

        D_flt.apply(p, t1);
        D_flt.dagger(t1, t2);
        psi[ALL] = p[X];
      }
    }

    // Add new solution to the list
    void save_new_solution(field<vector_type> & psi){
      if(MRE_size > 0){
        for(int i=1; i<MRE_size; i++){
          old_chi_inv[i] = old_chi_inv[i-1];
        }
        old_chi_inv[0] = psi;
      }
    }


    // Update the momentum with the derivative of the fermion
    // action
    void force_step(double eps){
      field<vector_type> psi, Mpsi;
      field<vector_type> Dhchi;
      psi.copy_boundary_condition(chi);
      Mpsi.copy_boundary_condition(chi);
      Dhchi.copy_boundary_condition(chi);
      field<momtype> force[NDIM], force2[NDIM];
      CG<DIRAC_OP> inverse(D);

      gauge.refresh();

      D_h.dagger(chi, Dhchi); 

      output0 << "light hasenbusch force\n";
      initial_guess(Dhchi, psi);
      inverse.apply(Dhchi, psi);
      save_new_solution(psi);

      D.apply(psi, Mpsi);

      Mpsi[D.par] = Mpsi[X] - chi[X];

      D.force(Mpsi, psi, force, 1);
      D.force(psi, Mpsi, force2, -1);

      foralldir(dir){
        force[dir][ALL] = -eps*(force[dir][X] + force2[dir][X]);
      }
      gauge.add_momentum(force);
    }

};





#endif