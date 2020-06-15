#ifndef FERMION_FIELD_H
#define FERMION_FIELD_H


#include "../dirac/conjugate_gradient.h"
#include <cmath>




template<typename momtype, typename DIRAC_OP>
class fermion_action{
  public:
    field<momtype> (&momentum)[NDIM];
    DIRAC_OP &D;
    using vector_type = typename DIRAC_OP::vector_type;
    field<vector_type> chi;


    fermion_action(DIRAC_OP &d, field<momtype> (&m)[NDIM])
    : D(d), momentum(m){ 
      chi = 0.0;
      #if NDIM > 3
      chi.set_boundary_condition(TUP, boundary_condition_t::ANTIPERIODIC);
      chi.set_boundary_condition(TDOWN, boundary_condition_t::ANTIPERIODIC);
      #endif
    }

    fermion_action(fermion_action &fa)
    : momentum(fa.momentum), D(fa.D)  {
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
    void backup_fields(){}

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
action_sum<action2, fermion_action<matrix, DIRAC_OP>> operator+(action2 a1, fermion_action<matrix, DIRAC_OP> a2){
  action_sum<action2, fermion_action<matrix, DIRAC_OP>> sum(a1, a2);
  return sum;
}


template<typename sun, typename representation, typename DIRAC_OP>
class high_representation_fermion_action {
  public:
    using momtype = squarematrix<representation::size, cmplx<typename representation::base_type>>;
    field<momtype> represented_force[NDIM];
    field<representation> (&represented_gauge)[NDIM];
    field<sun> (&momentum)[NDIM];
    field<sun> (&gauge)[NDIM];
    fermion_action<momtype, DIRAC_OP> base_action;

    high_representation_fermion_action(DIRAC_OP &d, field<sun> (&m)[NDIM], field<sun> (&g)[NDIM], field<representation> (&rg)[NDIM]) : momentum(m), gauge(g), represented_gauge(rg), 
    base_action(d, represented_force){};

    high_representation_fermion_action(high_representation_fermion_action &hrfa) : momentum(hrfa.momentum), gauge(hrfa.gauge), represented_gauge(hrfa.represented_gauge), base_action(hrfa.base_action) {}

    // Return the value of the action with the current
    // field configuration
    double action(){ 
      foralldir(dir){
        onsites(ALL){
          if(disable_avx[X]==0){};
          represented_gauge[dir][X].represent(gauge[dir][X]);
        }
      }
      return base_action.action();
    }

    // Make a copy of fields updated in a trajectory
    void backup_fields(){}

    // Restore the previous backup
    void restore_backup(){}

    /// Generate a pseudofermion field with a distribution given
    /// by the action chi 1/(D_dagger D) chi
    void draw_gaussian_fields(){
      foralldir(dir){
        onsites(ALL){
          if(disable_avx[X]==0){}; 
          represented_gauge[dir][X].represent(gauge[dir][X]);
        }
      }
      base_action.draw_gaussian_fields();
    }

    void force_step(double eps){
      foralldir(dir){
        onsites(ALL){
          if(disable_avx[X]==0){}; 
          represented_gauge[dir][X].represent(gauge[dir][X]);
          represented_force[dir][X] = 0;
        }
      }

      field<typename DIRAC_OP::vector_type> psi, Mpsi;
      psi.copy_boundary_condition(base_action.chi);
      Mpsi.copy_boundary_condition(base_action.chi);
      field<momtype> force[NDIM], force2[NDIM];
      CG<DIRAC_OP> inverse(base_action.D);

      psi[ALL]=0;
      inverse.apply(base_action.chi, psi);

      base_action.D.apply(psi, Mpsi);

      base_action.D.force(Mpsi, psi, force, 1);
      base_action.D.force(psi, Mpsi, force2, -1);

      foralldir(dir){
        onsites(ALL){
          if(disable_avx[X]==0){}; 
          // Disable for now. Return type of project_force needs to
          // be fixed
          force[dir][X] = force[dir][X] + force2[dir][X];
          element<sun> fforce;
          fforce = representation::project_force(force[dir][X]);
          momentum[dir][X] = momentum[dir][X] - eps*fforce;
        }
      }
    }
};

// Sum operator for creating an action_sum object
template<typename sun, typename matrix, typename DIRAC_OP, typename action2>
action_sum<action2, high_representation_fermion_action<sun, matrix, DIRAC_OP>> operator+(action2 a1, high_representation_fermion_action<sun, matrix, DIRAC_OP> a2){
  action_sum<action2, high_representation_fermion_action<sun, matrix, DIRAC_OP>> sum(a1, a2);
  return sum;
}


#endif