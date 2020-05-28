#ifndef FERMION_FIELD_H
#define FERMION_FIELD_H


#include "../../plumbing/algorithms/conjugate_gradients.h"
#include <cmath>




template<typename matrix, typename DIRAC_OP>
class fermion_action{
  public:
    field<matrix> (&gauge)[NDIM];
    field<matrix> (&momentum)[NDIM];
    DIRAC_OP &D;
    using vector_type = typename DIRAC_OP::vector_type;
    field<vector_type> chi;


    fermion_action(DIRAC_OP &d, field<matrix> (&g)[NDIM], field<matrix> (&m)[NDIM], parity p=ALL)
    : D(d), gauge(g), momentum(m){ 
      chi = 0.0;
      #if NDIM > 3
      chi.set_boundary_condition(TUP, boundary_condition_t::ANTIPERIODIC);
      chi.set_boundary_condition(TDOWN, boundary_condition_t::ANTIPERIODIC);
      #endif
    }

    fermion_action(fermion_action &fa)
    : gauge(fa.gauge), momentum(fa.momentum), D(fa.D)  {
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
    
      psi=0;
      inverse.apply(chi,psi);
      onsites(D.par){
        action += chi[X].rdot(psi[X]);
        
      }
      return action;
    }

    // Make a copy of fields updated in a trajectory
    void back_up_fields(){}

    // Restore the previous backup
    void restore_backup(){}

    /// Generate a pseudofermion field with a distribution given
    /// by the action chi 1/(D_dagger D) chi
    void draw_gaussian_fields(){
      field<vector_type> psi;
      psi.copy_boundary_condition(chi);
      onsites(D.par){
        if(disable_avx[X]==0){};
        psi[X].gaussian();
      }
      psi.mark_changed(ALL);
      D.dagger(psi,chi);
    }

    // Update the momentum with the derivative of the fermion
    // action
    void force_step(double eps){
      field<vector_type> psi, Mpsi;
      psi.copy_boundary_condition(chi);
      Mpsi.copy_boundary_condition(chi);
      field<matrix> force[NDIM], force2[NDIM];
      CG<DIRAC_OP> inverse(D);

      psi[ALL]=0;
      inverse.apply(chi, psi);

      D.apply(psi, Mpsi);

      D.force(Mpsi, psi, force, 1);
      D.force(psi, Mpsi, force2, -1);

      foralldir(dir){
        onsites(ALL){
          force[dir][X] = force[dir][X] + force2[dir][X];
          project_antihermitean(force[dir][X]);
          momentum[dir][X] = momentum[dir][X] - eps*force[dir][X];
        }
      }
    }

};

// Sum operator for creating an action_sum object
template<typename matrix, typename DIRAC_OP, typename action2>
action_sum<fermion_action<matrix, DIRAC_OP>, action2> operator+(fermion_action<matrix, DIRAC_OP> a1, action2 a2){
  action_sum<fermion_action<matrix, DIRAC_OP>, action2> sum(a1, a2);
  return sum;
}



#endif