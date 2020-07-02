#include "test.h"
#include "datatypes/vector.h"
#include "datatypes/sun.h"
#include "datatypes/representations.h"
#include "plumbing/field.h"
#include "hmc/update.h"
#include "hmc/gauge_field.h"
#include "hmc/smearing.h"
#include "dirac/wilson.h"
#include "dirac/staggered.h"
#include "hmc/fermion_field.h"

constexpr int N = 2;




template<typename dirac, typename gauge_field_type>
void check_forces(dirac &D, gauge_field_type &gauge){
  using sun = typename gauge_field_type::fund_type;
  using forcetype = squarematrix<gauge_field_type::N, cmplx<double>>;
  field<forcetype> force[NDIM];
  double eps = 1e-5;

  fermion_action fa(D, gauge);

  for(int ng = 0; ng < 1 ; ng++ ){ //matrix::generator_count(); ng++){
    gauge.set_unity();
    fa.draw_gaussian_fields();
    gauge.zero_momentum();
    foralldir(dir){
      force[dir][ALL] = 0;
    }

    sun g1 = gauge.get_gauge(0).get_value_at(50);
    sun h = sun(1) + cmplx(0, eps) * sun::generator(ng);
    sun g12 = h*g1;

    field<typename dirac::vector_type> psi, chi, tmp, tmp2;
    #if NDIM > 3
    psi.set_boundary_condition(TUP, boundary_condition_t::ANTIPERIODIC);
    psi.set_boundary_condition(TDOWN, boundary_condition_t::ANTIPERIODIC);
    #endif
    chi.copy_boundary_condition(psi);
    tmp.copy_boundary_condition(psi);
    tmp2.copy_boundary_condition(psi);

    onsites(D.par){
      if(disable_avx[X]==0){};
      psi[X].gaussian();
      chi[X].gaussian();
    }
    
    double s1 = 0;
    D.apply(psi,tmp);
    onsites(D.par){
      s1 += chi[X].rdot(tmp[X]);
    }

    if(mynode()==0){
      gauge.get_gauge(0).set_value_at(g12,50);
    }
    gauge.get_gauge(0).mark_changed(ALL);
    gauge.refresh();

    double s2 = 0;
    D.apply(psi,tmp);
    onsites(D.par){
      s2 += chi[X].rdot(tmp[X]);
    }

    if(mynode()==0)
      gauge.get_gauge(0).set_value_at(g1, 50);
    gauge.get_gauge(0).mark_changed(ALL);
    gauge.refresh();

    D.force(chi, psi, force, 1);

    gauge.add_momentum(force);
    sun f = gauge.get_momentum(0).get_value_at(50);
    double f1 = (s2-s1)/eps;
    double f2 = (f*cmplx(0,1)*sun::generator(ng)).trace().re;
    double diff = f2-f1;

    if(mynode()==0) {
      hila::output << "Action 1 " << s1 << "\n";
      hila::output << "Action 2 " << s2 << "\n";
      hila::output << "Calculated deriv " << f2 << "\n";
      hila::output << "Actual deriv " << (s2-s1)/eps << "\n";
      hila::output << "Fermion deriv " << ng << " diff " << diff << "\n";
      assert( diff*diff < eps*10 && "Fermion deriv" );
    }



    onsites(D.par){
      if(disable_avx[X]==0){};
      psi[X].gaussian();
      chi[X].gaussian();
    }
    gauge.zero_momentum();

    s1 = 0;
    D.dagger(psi,tmp);
    onsites(D.par){
      s1 += chi[X].rdot(tmp[X]);
    }

    if(mynode()==0){
      gauge.get_gauge(0).set_value_at(g12,50);
    }
    gauge.get_gauge(0).mark_changed(ALL);
    gauge.refresh();

    s2 = 0;
    D.dagger(psi,tmp);
    onsites(D.par){
      s2 += chi[X].rdot(tmp[X]);
    }

    if(mynode()==0)
      gauge.get_gauge(0).set_value_at(g1, 50);
    gauge.get_gauge(0).mark_changed(ALL);
    gauge.refresh();

    D.force(chi, psi, force, -1);
    gauge.add_momentum(force);
    f = gauge.get_momentum(0).get_value_at(50);
    f1 = (s2-s1)/eps;
    f2 = (f*cmplx(0,1)*sun::generator(ng)).trace().re;
    diff = f2-f1;

    if(mynode()==0) {
      //hila::output << "Action 1 " << s1 << "\n";
      //hila::output << "Action 2 " << s2 << "\n";
      //hila::output << "Calculated deriv " << f2 << "\n";
      //hila::output << "Actual deriv " << (s2-s1)/eps << "\n";
      //hila::output << "Fermion dg deriv " << ng << " diff " << diff << "\n";
      assert( diff*diff < eps*10 && "Fermion dg deriv" );
    }


    gauge.zero_momentum();

    if(mynode()==0){
      gauge.get_gauge(0).set_value_at(g12,50);
    }
    gauge.get_gauge(0).mark_changed(ALL);
    gauge.refresh();
    s2 = fa.action();

    if(mynode()==0)
      gauge.get_gauge(0).set_value_at(g1, 50);
    gauge.get_gauge(0).mark_changed(ALL);
    gauge.refresh();
    s1 = fa.action();

    fa.force_step(1.0);
    f = gauge.get_momentum(0).get_value_at(50);
    f1 = (s2-s1)/eps;
    f2 = (f*cmplx(0,1)*sun::generator(ng)).trace().re;
    diff = f2-f1;

    if(mynode()==0) {
      //hila::output << "Action 1 " << s1 << "\n";
      //hila::output << "Action 2 " << s2 << "\n";
      //hila::output << "Calculated force " << f2 << "\n";
      //hila::output << "Actual force " << (s2-s1)/eps << "\n";
      //hila::output << "Fermion force " << ng << " diff " << diff << "\n";
      assert( diff*diff < eps*10 && "Fermion force" );
    }
  }
}







int main(int argc, char **argv){

  /* Use a smaller lattice size since
   * the inversion takes a while */
  #if NDIM==1
  lattice->setup( 64, argc, argv );
  #elif NDIM==2
  lattice->setup( 32, 16, argc, argv );
  #elif NDIM==3
  lattice->setup( 32, 8, 8, argc, argv );
  #elif NDIM==4
  lattice->setup( 16, 8, 8, 8, argc, argv );
  #endif
  seed_random(2);

  // Test the force calculation by varying one gauge link
  // (Needs to be moved to tests)
  gauge_field<SU<N,double>> gauge;
  gauge.random();
  double eps = 1e-5;

  gauge_action ga(gauge, 1.0);
  
  for(int ng = 0; ng < SU<N>::generator_count(); ng++){
    foralldir(dir){
      onsites(ALL){
        gauge.momentum[dir][X] = 0;
      }
    }
    SU<N> g1 = gauge.gauge[0].get_value_at(50);
    SU<N> h = SU<N>(1) + cmplx(0,eps)* SU<N>::generator(ng);
    SU<N> g12 = h*g1;

    double s1 = ga.action();

    if(mynode()==0)
      gauge.gauge[0].set_value_at(g12,50);
    gauge.gauge[0].mark_changed(ALL);

    double s2 = ga.action();

    if(mynode()==0)
      gauge.gauge[0].set_value_at(g1, 50);
    gauge.gauge[0].mark_changed(ALL);

    ga.force_step(1.0);
    SU<N> f = gauge.momentum[0].get_value_at(50);
    double diff = (f*cmplx(0,1)*SU<N>::generator(ng)).trace().re - (s2-s1)/eps;

    if(mynode()==0) {
      //hila::output << "Force 1 " << (f*SU<N>::generator(ng)).trace().re << "\n";
      //hila::output << "Force 2 " << (s2-s1)/eps << "\n";
      //hila::output << "Force " << ng << " diff " << diff << "\n";
      h = cmplx(0,1)*SU<N>::generator(ng);
      assert( diff*diff < eps*10 && "Gauge force" );
    }
  }


  using SUN=SU<N, double>;
  using adj=adjoint<N, double>;
  using sym=symmetric<N, double>;
  using asym=antisymmetric<N, double>;

  // Define represented gauge fields for the Dirac force tests
  represented_gauge_field<adj> adj_gauge(gauge);
  represented_gauge_field<sym> sym_gauge(gauge);
  represented_gauge_field<asym> asym_gauge(gauge);
  stout_smeared_field<SUN> stout_gauge(gauge, 0.01);


  // Staggered Forces
  output0 << "Checking staggered forces:\n";
  dirac_staggered D_Stg(1.5, gauge);
  check_forces(D_Stg, gauge);

  output0 << "Checking evenodd preconditioned staggered forces:\n";
  dirac_staggered_evenodd D_Stg_eo(1.5, gauge);
  check_forces(D_Stg_eo, gauge);

  output0 << "Checking stout smeared forces:\n";
  dirac_staggered_evenodd D_stg_stout(1.5, stout_gauge);
  check_forces(D_stg_stout, stout_gauge);

  output0 << "Checking adjoint staggered forces:\n";
  dirac_staggered_evenodd D_stg_adj(1.5, adj_gauge);
  check_forces(D_stg_adj, adj_gauge);

  output0 << "Checking symmetric staggered forces:\n";
  dirac_staggered_evenodd D_stg_sym(1.5, sym_gauge);
  check_forces(D_stg_sym, sym_gauge);

  output0 << "Checking antisymmetric staggered forces:\n";
  dirac_staggered_evenodd D_stg_asym(1.5, asym_gauge);
  check_forces(D_stg_asym, asym_gauge);
  

  // Wilson forces
  output0 << "Checking evenodd preconditioned Wilson forces:\n";
  Dirac_Wilson D_W(0.05, gauge);
  check_forces(D_W, gauge);

  output0 << "Checking evenodd preconditioned Wilson forces:\n";
  Dirac_Wilson_evenodd D_W_eo(0.05, gauge);
  check_forces(D_W_eo, gauge);

  output0 << "Checking adjoint Wilson forces:\n";
  Dirac_Wilson_evenodd D_W_adj(0.05, adj_gauge);
  check_forces(D_W_adj, adj_gauge);

  output0 << "Checking symmetric Wilson forces:\n";
  Dirac_Wilson_evenodd D_W_sym(0.05, sym_gauge);
  check_forces(D_W_sym, sym_gauge);

  output0 << "Checking antisymmetric Wilson forces:\n";
  Dirac_Wilson_evenodd D_W_asym(0.05, asym_gauge);
  check_forces(D_W_asym, asym_gauge);


  // Check also the momentum action and derivative
  for(int ng = 0; ng < SU<N>::generator_count(); ng++){
    ga.draw_gaussian_fields();

    double s1 = ga.action();
    SU<N> h = gauge.momentum[0].get_value_at(0);
    h += eps * cmplx(0,1)*SU<N>::generator(ng);
    if(mynode()==0)
      gauge.momentum[0].set_value_at(h, 0);
    double s2 = ga.action();

    double diff = (h*cmplx(0,1)*SU<N>::generator(ng)).trace().re + (s2-s1)/eps;
    if(mynode()==0) {
      //hila::output << "Mom 1 " << (h*SU<N>::generator(ng)).trace().re << "\n";
      //hila::output << "Mom 2 " << (s2-s1)/eps << "\n";
      //hila::output << "Mom " << ng << " diff " << diff << "\n";
      h = cmplx(0,1)*SU<N>::generator(ng);
      assert( diff*diff < eps*10 && "Momentum derivative" );
    }
  }

  finishrun();
}

