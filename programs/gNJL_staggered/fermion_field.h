

template<typename SUN, typename VECTOR, typename DIRAC_OP>
void fermion_force_gNJL(field<VECTOR> &chi, field<SUN> (&momentum)[NDIM], field<double> &sigma_momentum, field<double> &pi_momentum, DIRAC_OP &D, double eps){
  field<VECTOR> psi, Mpsi;
  psi.copy_boundary_condition(chi);
  Mpsi.copy_boundary_condition(chi);
  field<SUN> force[NDIM], force2[NDIM];
  field<double> smom, smom2, pmom, pmom2;
  CG<DIRAC_OP> inverse(D);
  
  psi=0;
  inverse.apply(chi, psi);
  
  D.apply(psi, Mpsi);

  D.force(Mpsi, psi, force, smom, pmom, 1);
  D.force(psi, Mpsi, force2, smom2, pmom2, -1);

  foralldir(dir){
    onsites(ALL){
      force[dir][X] = force[dir][X] + force2[dir][X];
      project_antihermitean(force[dir][X]);
      momentum[dir][X] = momentum[dir][X] - eps*force[dir][X];
    }
  }

  onsites(ALL){
    smom[X] = smom[X] + smom2[X];
    pmom[X] = pmom[X] + pmom2[X];
    sigma_momentum[X] = sigma_momentum[X] - eps * smom[X];
    pi_momentum[X] = pi_momentum[X] - eps * pmom[X];
  }
}


template<typename vector, typename matrix>
class gNJL_fermion_action{
  public:
    field<matrix> (&momentum)[NDIM];
    field<double> &sigma_momentum, &pi_momentum;
    dirac_staggered_gNJL<vector, matrix> &D;
    field<vector> chi;


    gNJL_fermion_action(dirac_staggered_gNJL<vector, matrix> &d, field<matrix> (&m)[NDIM], field<double> &sm, field<double> &pm)
    : D(d), momentum(m), sigma_momentum(sm), pi_momentum(pm){ chi = 0.0; }

    gNJL_fermion_action(gNJL_fermion_action &fa)
    : momentum(fa.momentum), D(fa.D), sigma_momentum(fa.sigma_momentum), pi_momentum(fa.pi_momentum)  {
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
      fermion_force_gNJL( chi, momentum, sigma_momentum, pi_momentum, D, eps );
    }

};