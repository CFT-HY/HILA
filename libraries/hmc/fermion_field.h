#ifndef FERMION_FIELD_H
#define FERMION_FIELD_H


#include "../../libraries/dirac/Hasenbusch.h"
#include "../dirac/conjugate_gradient.h"
#include <cmath>




template<typename gauge_field, typename DIRAC_OP>
class fermion_action : action_base{
  public:
    using vector_type = typename DIRAC_OP::vector_type;
    using momtype = squarematrix<gauge_field::N, cmplx<typename gauge_field::basetype>>;
    gauge_field &gauge;
    DIRAC_OP &D;
    field<vector_type> chi;

    fermion_action(DIRAC_OP &d, gauge_field &g)
    : D(d), gauge(g){ 
      chi = 0.0;
      #if NDIM > 3
      chi.set_boundary_condition(TUP, boundary_condition_t::ANTIPERIODIC);
      chi.set_boundary_condition(TDOWN, boundary_condition_t::ANTIPERIODIC);
      #endif
    }

    fermion_action(fermion_action &fa)
    : gauge(fa.gauge), D(fa.D)  {
      chi = fa.chi;
      #if NDIM > 3
      chi.set_boundary_condition(TUP, boundary_condition_t::ANTIPERIODIC);
      chi.set_boundary_condition(TDOWN, boundary_condition_t::ANTIPERIODIC);
      #endif
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

    // Update the momentum with the derivative of the fermion
    // action
    void force_step(double eps){
      field<vector_type> psi, Mpsi;
      psi.copy_boundary_condition(chi);
      Mpsi.copy_boundary_condition(chi);
      field<momtype> force[NDIM], force2[NDIM];
      CG<DIRAC_OP> inverse(D);

      gauge.refresh();

      psi[ALL]=0;
      inverse.apply(chi, psi);

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

    Hasenbusch_action_1(DIRAC_OP &d, gauge_field &g, double mh) : _mh(mh), D_h(d, mh), base_action(D_h, g) {}

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

    Hasenbusch_action_2(DIRAC_OP &d, gauge_field &g, double _mh) : mh(_mh), D(d), D_h(d, _mh), gauge(g) {
      chi = 0.0;
      #if NDIM > 3
      chi.set_boundary_condition(TUP, boundary_condition_t::ANTIPERIODIC);
      chi.set_boundary_condition(TDOWN, boundary_condition_t::ANTIPERIODIC);
      #endif
    }

    Hasenbusch_action_2(Hasenbusch_action_2 &fa) : mh(fa.mh), D(fa.D), D_h(fa.D_h), gauge(fa.gauge){
      chi = fa.chi;
      #if NDIM > 3
      chi.set_boundary_condition(TUP, boundary_condition_t::ANTIPERIODIC);
      chi.set_boundary_condition(TDOWN, boundary_condition_t::ANTIPERIODIC);
      #endif
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
      CG inverse_h(D_h); // Applies 1/(D_h^dagger D_h)
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

      psi[ALL]=0;
      inverse.apply(Dhchi, psi);

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