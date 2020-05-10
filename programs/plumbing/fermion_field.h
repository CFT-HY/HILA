#ifndef FERMION_FIELD_H
#define FERMION_FIELD_H


#include "../plumbing/algorithms/conjugate_gradients.h"
#include <cmath>


/// Generate a pseudofermion field with a distribution given
/// by the action chi 1/(D_dagger D) chi
template<typename VECTOR, typename DIRAC_OP>
void generate_pseudofermion(field<VECTOR> &chi, DIRAC_OP D){
  field<VECTOR> psi, tmp;
  onsites(ALL){
    psi[X].gaussian();
  }
  D.dagger(psi,chi);
}


/// Calculate the action of a fermion term
template<typename VECTOR, typename DIRAC_OP>
double pseudofermion_action(field<VECTOR> &chi, DIRAC_OP D){
  field<VECTOR> psi, tmp;
  CG<field<VECTOR>, DIRAC_OP> inverse(D);
  double action = 0;

  psi=0;
  inverse.apply(chi,psi);
  onsites(ALL){
    action += (chi[X]*psi[X]).re;
  }
  return action;
}



/// Apply the force of the gauge field on the momentum field 
template<typename SUN, typename VECTOR, typename DIRAC_OP>
void fermion_force(field<VECTOR> &chi, field<SUN> (&momentum)[NDIM], DIRAC_OP &D, double eps){
  field<VECTOR> psi, Mpsi;
  field<SUN> force[NDIM], force2[NDIM];
  CG<field<VECTOR>, DIRAC_OP> inverse(D);
  
  psi=0;
  inverse.apply(chi, psi);
  
  D.apply(psi, Mpsi);

  D.force(Mpsi, psi, force);
  D.force(psi, Mpsi, force2);

  foralldir(dir){
    onsites(ALL){
      force[dir][X] = force[dir][X] - force2[dir][X];
      project_antihermitean(force[dir][X]);
      momentum[dir][X] = momentum[dir][X] - eps*force[dir][X];
    }
  }
}





template<typename SUN, typename DIRAC_OP>
class fermion_action{
  public:
    field<SUN> (&gauge)[NDIM];
    field<SUN> (&momentum)[NDIM];
    DIRAC_OP &D;
    field<typename DIRAC_OP::vector_type> chi;


    fermion_action(DIRAC_OP &d, field<SUN> (&g)[NDIM], field<SUN> (&m)[NDIM])
    : D(d), gauge(g), momentum(m){ chi = 0.0; }

    fermion_action(fermion_action &fa)
    : gauge(fa.gauge), momentum(fa.momentum), D(fa.D)  {
      chi = fa.chi;
    }

    // Return the value of the action with the current
    // field configuration
    double action(){ 
      return pseudofermion_action(chi, D);
    }

    // Make a copy of fields updated in a trajectory
    void back_up_fields(){}

    // Restore the previous backup
    void restore_backup(){}

    /// Gaussian random momentum for each element
    void draw_gaussian_fields(){
      generate_pseudofermion(chi, D);
    }

    // Update the momentum with the derivative of the fermion
    // action
    void force_step(double eps){
      fermion_force( chi, momentum, D, eps );
    }

};

// Sum operator for creating an action_sum object
template<typename SUN, typename DIRAC_OP, typename action2>
action_sum<fermion_action<SUN, DIRAC_OP>, action2> operator+(fermion_action<SUN, DIRAC_OP> a1, action2 a2){
  action_sum sum(a1, a2);
  return sum;
}



#endif