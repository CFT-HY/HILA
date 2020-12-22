// Representations introduces a datatype and must be before field.h
#include "datatypes/representations.h"
// Test.h includes field.h
#include "test.h"
#include "hmc/hmc.h"
#include "hmc/gauge_field.h"
#include "hmc/smearing.h"
#include "dirac/wilson.h"
#include "dirac/staggered.h"
#include "hmc/fermion_field.h"

constexpr int N = 2;




template<typename fermion_action, typename dirac, typename gauge_field_type>
void check_forces(fermion_action &fa, dirac &D, gauge_field_type &gauge){
  using sun = typename gauge_field_type::fund_type;
  using forcetype = SquareMatrix<gauge_field_type::N, cmplx<double>>;
  field<forcetype> force[NDIM];
  double eps = 1e-5;

  for(int ng = 0; ng < 1 ; ng++ ){ //matrix::generator_count(); ng++){
    //gauge.set_unity();
    fa.draw_gaussian_fields();
    gauge.zero_momentum();
    foralldir(dir){
      force[dir][ALL] = 0;
    }

    coordinate_vector coord(0);

    sun g1 = gauge.get_gauge(0).get_element(coord);
    sun h = sun(1) + cmplx<double>(0, eps) * sun::generator(ng);
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

    gauge.get_gauge(0).set_element(g12, coord);
    gauge.refresh();

    double s2 = 0;
    D.apply(psi,tmp);
    onsites(D.par){
      s2 += chi[X].rdot(tmp[X]);
    }

    gauge.get_gauge(0).set_element(g1, coord);
    gauge.refresh();

    D.force(chi, psi, force, 1);

    gauge.add_momentum(force);
    sun f = gauge.get_momentum(0).get_element(coord);
    double f1 = (s2-s1)/eps;
    double f2 = (f*cmplx<double>(0,1)*sun::generator(ng)).trace().re;
    double diff = f2-f1;

    if(mynode()==0) {
      //output0 << "Action 1 " << s1 << "\n";
      //output0 << "Action 2 " << s2 << "\n";
      //output0 << "Calculated deriv " << f2 << "\n";
      //output0 << "Actual deriv " << (s2-s1)/eps << "\n";
      //output0 << "Fermion deriv " << ng << " diff " << diff << "\n";
      assert( diff*diff < eps && "Fermion deriv" );
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

    gauge.get_gauge(0).set_element(g12, coord);
    gauge.refresh();

    s2 = 0;
    D.dagger(psi,tmp);
    onsites(D.par){
      s2 += chi[X].rdot(tmp[X]);
    }

    gauge.get_gauge(0).set_element(g1, coord);
    gauge.refresh();

    D.force(chi, psi, force, -1);
    gauge.add_momentum(force);
    f = gauge.get_momentum(0).get_element(coord);
    f1 = (s2-s1)/eps;
    f2 = (f*cmplx<double>(0,1)*sun::generator(ng)).trace().re;
    diff = f2-f1;

    if(mynode()==0) {
      //output0 << "Action 1 " << s1 << "\n";
      //output0 << "Action 2 " << s2 << "\n";
      //output0 << "Calculated deriv " << f2 << "\n";
      //output0 << "Actual deriv " << (s2-s1)/eps << "\n";
      //output0 << "Fermion dg deriv " << ng << " diff " << diff << "\n";
      assert( diff*diff < eps && "Fermion dg deriv" );
    }


    gauge.zero_momentum();
    field<double> sf1, sf2;
    sf1[ALL]=0; sf2[ALL]=0;

    gauge.get_gauge(0).set_element(g12, coord);
    gauge.refresh();
    fa.action(sf2);

    gauge.get_gauge(0).set_element(g1, coord);
    gauge.refresh();
    fa.action(sf1);

    double dS = 0;
    s1 = s2 = 0;
    onsites(ALL){
      dS += sf2[X]-sf1[X];
      s1 += sf1[X];
      s2 += sf2[X];
    }


    fa.force_step(1.0);
    f = gauge.get_momentum(0).get_element(coord);
    f1 = dS/eps;
    f2 = (f*cmplx<double>(0,1)*sun::generator(ng)).trace().re;
    diff = f2-f1;

    if(mynode()==0) {
      //output0 << "Action 1 " << s1 << "\n";
      //output0 << "Action 2 " << s2 << " " << dS << " " << dS/s2 << "\n";
      //output0 << "Calculated force " << f2 << "\n";
      //output0 << "Actual force " << f1 << "\n";
      //output0 << "Fermion force " << ng << " diff " << diff << "\n";
      assert( diff*diff < eps && "Fermion force" );
    }
  }
}







int main(int argc, char **argv){

  /* Use a smaller lattice size since
   * the inversion takes a while */
  #if NDIM==1
  const int nd[NDIM] = {64};
  #elif NDIM==2
  const int nd[NDIM] = {32,8};
  #elif NDIM==3
  const int nd[NDIM] = {16,8,8};
  #elif NDIM==4
  const int nd[NDIM] = {16,8,8,8};
  #endif
  hila::initialize(argc,argv);
  lattice->setup(nd);
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
    SU<N> h = SU<N>(1) + cmplx<double>(0,eps)* SU<N>::generator(ng);
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
    double diff = (f*cmplx<double>(0,1)*SU<N>::generator(ng)).trace().re - (s2-s1)/eps;

    if(mynode()==0) {
      //hila::output << "Force " << (f*cmplx<double>(0,1)*SU<N>::generator(ng)).trace().re << "\n";
      //hila::output << "Force " << (f*SU<N>::generator(ng)).trace().re << "\n";
      //hila::output << "Deriv " << (s2-s1)/eps << "\n";
      //hila::output << "Force " << ng << " diff " << diff << "\n";
      h = cmplx<double>(0,1)*SU<N>::generator(ng);
      assert( diff*diff < eps*10 && "Gauge force" );
    }
  }


  using SUN=SU<N, double>;
  using adj=adjointRep<N, double>;
  using sym=symmetric<N, double>;
  using asym=antisymmetric<N, double>;

  // Define represented gauge fields for the Dirac force tests
  gauge.random();
  represented_gauge_field<adj> adj_gauge(gauge);
  represented_gauge_field<sym> sym_gauge(gauge);
  represented_gauge_field<asym> asym_gauge(gauge);
  stout_smeared_field<SUN> stout_gauge(gauge, 0.1, 4, 10);
  HEX_smeared_field<SUN> hex_gauge(gauge, 10);

  // Staggered Forces
  output0 << "Checking staggered forces:\n";
  {
    dirac_staggered D(5.0, gauge);
    fermion_action fa(D, gauge);
    check_forces(fa, D, gauge);
  }

  {
    output0 << "Checking evenodd preconditioned staggered forces:\n";
    dirac_staggered_evenodd D(5.0, gauge);
    fermion_action fa(D, gauge);
    check_forces(fa, D, gauge);
  }

  {
    output0 << "Checking stout smeared forces:\n";
    dirac_staggered_evenodd D(5.0, stout_gauge);
    fermion_action fa(D, stout_gauge);
    check_forces(fa, D, stout_gauge);
  }

  {
    output0 << "Checking HEX smeared forces:\n";
    dirac_staggered_evenodd D(5.0, hex_gauge);
    fermion_action fa(D, hex_gauge);
    check_forces(fa, D, hex_gauge);
  }

  {
    output0 << "Checking adjoint staggered forces:\n";
    dirac_staggered_evenodd D(5.0, adj_gauge);
    fermion_action fa(D, adj_gauge);
    check_forces(fa, D, adj_gauge);
  }

  {
    output0 << "Checking symmetric staggered forces:\n";
    dirac_staggered_evenodd D(5.0, sym_gauge);
    fermion_action fa(D, sym_gauge);
    check_forces(fa, D, sym_gauge);
  }

  {
    output0 << "Checking antisymmetric staggered forces:\n";
    dirac_staggered_evenodd D(5.0, asym_gauge);
    fermion_action fa(D, asym_gauge);
    check_forces(fa, D, asym_gauge);
  }

  // Wilson forces
  {
    output0 << "Checking Wilson forces:\n";
    Dirac_Wilson D(0.05, gauge);
    fermion_action fa(D, gauge);
    check_forces(fa, D, gauge);
  }

  {
    output0 << "Checking evenodd Wilson forces:\n";
    Dirac_Wilson_evenodd D(0.12, gauge);
    fermion_action fa(D, gauge);
    check_forces(fa, D, gauge);

    output0 << "Checking hasenbusch 1:\n";
    Hasenbusch_action_1 fa1(D, gauge, 0);
    check_forces(fa1, D, gauge);
 
    output0 << "Checking hasenbusch 2:\n";
    Hasenbusch_action_2 fa2(D, gauge, 0);
    check_forces(fa2, D, gauge);
  }

  {
    output0 << "Checking adjoint Wilson forces:\n";
    Dirac_Wilson_evenodd D(0.05, adj_gauge);
    fermion_action fa(D, adj_gauge);
    check_forces(fa, D, adj_gauge);
  }

  {
    output0 << "Checking symmetric Wilson forces:\n";
    Dirac_Wilson_evenodd D(0.05, sym_gauge);
    fermion_action fa(D, sym_gauge);
    check_forces(fa, D, sym_gauge);
  }

  {
    output0 << "Checking antisymmetric Wilson forces:\n";
    Dirac_Wilson_evenodd D(0.05, asym_gauge);
    fermion_action fa(D, asym_gauge);
    check_forces(fa, D, asym_gauge);
  }

  // Check also the momentum action and derivative
  for(int ng = 0; ng < SU<N>::generator_count(); ng++){
    ga.draw_gaussian_fields();

    double s1 = ga.action();
    SU<N> h = gauge.momentum[0].get_value_at(0);
    h += eps * cmplx<double>(0,1)*SU<N>::generator(ng);
    if(mynode()==0)
      gauge.momentum[0].set_value_at(h, 0);
    double s2 = ga.action();

    double diff = (h*SU<N>::generator(ng)).trace().re + (s2-s1)/eps;
    if(mynode()==0) {
      //hila::output << "Mom 1 " << (h*SU<N>::generator(ng)).trace().re << "\n";
      //hila::output << "Mom 2 " << (s2-s1)/eps << "\n";
      //hila::output << "Mom " << ng << " diff " << diff << "\n";
      h = cmplx<double>(0,1)*SU<N>::generator(ng);
      assert( diff*diff < eps*10 && "Momentum derivative" );
    }
  }

  hila::finishrun();
}

