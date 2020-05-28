

template<typename SUN, typename VECTOR, typename DIRAC_OP>
void fermion_force_gNJL(field<VECTOR> &chi, field<SUN> (&momentum)[NDIM], field<double> &sigma_momentum, field<double> &pi_momentum, DIRAC_OP &D, double eps){
  
}


template<typename matrix, typename DIRAC_OP>
class gNJL_fermion_action : public fermion_action<matrix,DIRAC_OP> {
  public:
    field<double> &sigma_momentum, &pi_momentum;

    using vector_type = typename DIRAC_OP::vector_type;
    using fa = fermion_action<matrix,DIRAC_OP>;

    gNJL_fermion_action(DIRAC_OP &d, field<matrix> (&m)[NDIM], field<double> &sm, field<double> &pm)
    : fermion_action<matrix,DIRAC_OP>(d, m), sigma_momentum(sm), pi_momentum(pm){}

    gNJL_fermion_action(gNJL_fermion_action &fa)
    : fermion_action<matrix,DIRAC_OP>(fa.D, fa.momentum), sigma_momentum(fa.sigma_momentum), pi_momentum(fa.pi_momentum){}


    // Update the gauge momentum and the auxiliary field momenta
    void force_step(double eps){
      field<vector_type> psi, Mpsi;
      psi.copy_boundary_condition(fa::chi);
      Mpsi.copy_boundary_condition(fa::chi);
      field<matrix> force[NDIM], force2[NDIM];
      field<double> smom, smom2, pmom, pmom2;
      CG<DIRAC_OP> inverse(fa::D);

      psi=0;
      inverse.apply(fa::chi, psi);

      fa::D.apply(psi, Mpsi);

      fa::D.force(Mpsi, psi, force, smom, pmom, 1);
      fa::D.force(psi, Mpsi, force2, smom2, pmom2, -1);

      foralldir(dir){
        onsites(ALL){
          force[dir][X] = force[dir][X] + force2[dir][X];
          project_antihermitean(force[dir][X]);
          fa::momentum[dir][X] = fa::momentum[dir][X] - eps*force[dir][X];
        }
      }

      onsites(ALL){
        smom[X] = smom[X] + smom2[X];
        pmom[X] = pmom[X] + pmom2[X];
        sigma_momentum[X] = sigma_momentum[X] - eps * smom[X];
        pi_momentum[X] = pi_momentum[X] - eps * pmom[X];
      }
    }

};



