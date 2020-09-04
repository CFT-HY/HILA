#ifndef FERMION_FIELD_H
#define FERMION_FIELD_H


#include "../../libraries/dirac/Hasenbusch.h"
#include "../dirac/conjugate_gradient.h"
#include <cmath>




template<typename gauge_field, typename DIRAC_OP>
class fermion_action{
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

    // Make a copy of fields updated in a trajectory
    void backup_fields(){}

    // Restore the previous backup
    void restore_backup(){}

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


// Sum operator for creating an action_sum object
template<typename matrix, typename DIRAC_OP, typename action2>
action_sum<action2, fermion_action<matrix, DIRAC_OP>> operator+(action2 a1, fermion_action<matrix, DIRAC_OP> a2){
  action_sum<action2, fermion_action<matrix, DIRAC_OP>> sum(a1, a2);
  return sum;
}









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
class Hasenbusch_action_1 {
  public:
    Hasenbusch_operator<DIRAC_OP> D_h;
    fermion_action<gauge_field, Hasenbusch_operator<DIRAC_OP>> base_action;
    double _mh;

    Hasenbusch_action_1(DIRAC_OP &d, gauge_field &g, double mh) : _mh(mh), D_h(d, mh), base_action(D_h, g) {}

    Hasenbusch_action_1(Hasenbusch_action_1 &fa) : _mh(fa._mh), D_h(fa.D_h), base_action(fa.base_action){}

    double action(){return(base_action.action());}
    void backup_fields(){}
    void restore_backup(){}
    void draw_gaussian_fields(){base_action.draw_gaussian_fields();}
    void force_step(double eps){base_action.force_step(eps);}

};


template<typename gauge_field, typename DIRAC_OP>
class Hasenbusch_action_2 {
  public:
    Hasenbusch_operator<DIRAC_OP> D_h;
    fermion_action<gauge_field, DIRAC_OP> base_action;
    double _mh;

    Hasenbusch_action_2(DIRAC_OP &d, gauge_field &g, double mh) : _mh(mh), D_h(d, mh), base_action(d, g) {}

    Hasenbusch_action_2(Hasenbusch_action_2 &fa) : _mh(fa._mh), D_h(fa.D_h), base_action(fa.base_action){}

    double action(){return(0);}
    void backup_fields(){}
    void restore_backup(){}
    void draw_gaussian_fields(){base_action.draw_gaussian_fields();}
    void force_step(double eps){}

};





#endif