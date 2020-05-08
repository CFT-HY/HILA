#ifndef FERMION_FIELD_H
#define FERMION_FIELD_H



/// Generate a pseudofermion field with a distribution given
/// by the action chi 1/(D_dagger D) chi
template<typename VECTOR, typename DIRAC_OP>
void generate_pseudofermion(field<VECTOR> &chi, DIRAC_OP D){
  field<VECTOR> psi, tmp;
  onsites(ALL){
    psi[X].random();
  }
  D.apply(psi,tmp);
  D.dagger(tmp,chi);
}




/// Apply the force of the gauge field on the momentum field 
template<typename gauge_action_type, typename VECTOR, typename DIRAC_OP>
void fermion_force(field<VECTOR> &chi, gauge_action_type gauge, double eps){
  
}





template<typename gauge_action_type, typename VECTOR, typename DIRAC_OP>
class fermion_action{
  public:
    gauge_action_type &ga;
    DIRAC_OP &D;
    field<SUN> (&gauge)[NDIM];
    field<VECTOR> chi;

    fermion_action(gauge_action_type &g, DIRAC_OP &d) : ga(g), D(d),gauge(g.gauge) {
      chi = 0.0;
    }

    // Return the value of the action with the current
    // field configuration
    double action(){
      return ga.action();
    }

    // Make a copy of fields updated in a trajectory
    void back_up_fields(){
      ga.back_up_fields();
    }

    // Restore the previous backup
    void restore_backup(){
      ga.restore_backup();
    }

    /// Gaussian random momentum for each element
    void generate_momentum(){
      ga.generate_momentum();
      generate_pseudofermion(chi, D);
    }

    // Update the momentum with the derivative of the fermion
    // action
    void force_step(double eps){
      
    }

    // Update the gauge field with momentum
    void momentum_step(double eps){
      ga.step(eps);
    }

};





#endif