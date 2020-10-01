#ifndef FERMION_FIELD_H
#define FERMION_FIELD_H

#include "gauge_field.h"
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
    bool use_MRE_guess = true;

    // We save a few previous invertions to build an initial guess.
    // old_chi contains a list of these
    int MRE_size = 4;
    std::vector<field<vector_type>> old_chi_inv;

    fermion_action(DIRAC_OP &d, gauge_field &g)
    : D(d), gauge(g){ 
      chi = 0.0;
      #if NDIM > 3
      chi.set_boundary_condition(TUP, boundary_condition_t::ANTIPERIODIC);
      chi.set_boundary_condition(TDOWN, boundary_condition_t::ANTIPERIODIC);
      #endif
      old_chi_inv.resize(MRE_size);
      for(int i=0; i<MRE_size; i++){
        old_chi_inv[i][ALL] = 0;
      }
    }

    fermion_action(fermion_action &fa)
    : gauge(fa.gauge), D(fa.D)  {
      chi = fa.chi;
      #if NDIM > 3
      chi.set_boundary_condition(TUP, boundary_condition_t::ANTIPERIODIC);
      chi.set_boundary_condition(TDOWN, boundary_condition_t::ANTIPERIODIC);
      #endif
      old_chi_inv.resize(MRE_size);
      for(int i=0; i<MRE_size; i++){
        old_chi_inv[i][ALL] = 0;
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
    void MRE_guess(field<vector_type> & psi, field<vector_type> & ){
      double M[MRE_size][MRE_size];
      double v[MRE_size];
      field<vector_type> basis[MRE_size];
      field<vector_type> tmp;

      // Build an orthogonal basis from the previous solutions
      for(int i=0; i<MRE_size; i++) {
        // Start with the original solution vector
        basis[i][ALL] = old_chi_inv[i][X];
        // Remove the projected components of all previous vectors
        for(int j=0; j<i; j++){
          double vdot = 0, norm = 0;
          onsites(D.par) {
            norm += basis[i][X].rdot(basis[i][X]);
            vdot += basis[j][X].rdot(basis[i][X]);
          }
          if(norm*norm > 1e-32){
            basis[i][D.par] -= vdot/norm*basis[j][X];
          }
        }
      }

      // Build the projected matrix, M[i][j] = v[i].v[j]
      for(int i=0; i<MRE_size; i++) {
        field<vector_type> Dchi, DDchi;
        D.apply(basis[i], Dchi);
        D.dagger(Dchi, DDchi);
        for(int j=0; j<MRE_size; j++) {
          double sum = 0;
          onsites(D.par) {
            sum += basis[j][X].rdot(DDchi[X]);
          }
          M[j][i] = sum;
        }
      }
      // And the projected vector
      for(int i=0; i<MRE_size; i++){
        double sum = 0;
        onsites(D.par){
          sum += basis[i][X].rdot(chi[X]);
        }
        v[i] = sum;
      }

      // Now invert the small matrix M (Gaussian elimination)
      for(int i=0; i<MRE_size; i++){
        // Check that the diagonal element is nonzero
        if(M[i][i]*M[i][i] > 1e-32){ 
          // Normalize the i:th row
          double diag_inv = 1.0/M[i][i];
          for(int j=0; j<MRE_size; j++) {
            M[i][j] *= diag_inv;
          }
          v[i] *= diag_inv;
          // Subtract from all other rows
          for(int k=0; k<MRE_size; k++) if(k!=i){
            double weight = M[k][i];
            for(int j=0; j<MRE_size; j++) {
              M[k][j] -= weight*M[i][j];
            }
            v[k] -= weight*v[i];
          }
        } else {
          // In the matrix element is too small, just remove it from the basis
          v[i] = 0;
          for(int j=0; j<MRE_size; j++) { 
            M[j][i] = 0;
          }
        }
      }

      // Construct the solution in the original basis
      psi[ALL] = 0;
      for(int i=0; i<MRE_size; i++) if(!isnan(v[i])){
        psi[D.par] += v[i]*basis[i][X];
      }
    }

    // Add new solution to the list
    void save_new_solution(field<vector_type> & psi){
      for(int i=1; i<MRE_size; i++){
        old_chi_inv[i] = old_chi_inv[i-1];
      }
      old_chi_inv[0] = psi;
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

      if(MRE_size > 0){
        MRE_guess(psi, chi);
      } else {
        psi[ALL]=0;
      }
      inverse.apply(chi, psi);
      if(MRE_size > 0){
        save_new_solution(psi);
      }

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