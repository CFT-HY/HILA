/*********************************************
 * Simulates Wilson fermions with a gauge *
 * interaction                               *
 *********************************************/


#include "Wilson.h"


void test_gamma_matrices(){
  Wilson_vector<N> w1, w2, w3;
  half_Wilson_vector<N> h1;
  SU<N> U; U.random();
  w1.gaussian();

  w2 = w1-gamma5*(gamma5*w1);
  assert(w2.norm_sq() < 0.0001 && "g5*g5 = 1");

  foralldir(d){
    w2 = w1-gamma_matrix[d]*(gamma_matrix[d]*w1);
    assert(w2.norm_sq() < 0.0001 && "gamma_d*gamma_d = 1");

    w2 = w1 + gamma_matrix[d]*w1;
    h1 = half_Wilson_vector(w1,d);
    double diff = w2.norm_sq() - h1.norm_sq();
    assert(diff*diff < 0.0001 && "half_Wilson_vector projection norm to direction XUP");

    w3 = (U*h1).expand(d) - U*w2;
    assert(w3.norm_sq() < 0.0001 && "half_wilson_vector expand");


    w2 = w1 - gamma_matrix[d]*w1;
    h1 = half_Wilson_vector(w1,opp_dir(d));
    diff = w2.norm_sq() - h1.norm_sq();
    assert(diff*diff < 0.0001 && "half_Wilson_vector projection norm to direction XUP");

    w3 = (U*h1).expand(opp_dir(d)) - U*w2;
    assert(w3.norm_sq() < 0.0001 && "half_wilson_vector expand");

  }

}



int main(int argc, char **argv){

  input parameters = input();
  parameters.import("parameters");
  double beta = parameters.get("beta");
  double kappa = parameters.get("kappa");
  int seed = parameters.get("seed");
	double hmc_steps = parameters.get("hmc_steps");
	double traj_length = parameters.get("traj_length");

  lattice->setup( nd[0], nd[1], nd[2], nd[3], argc, argv );
  seed_random(seed);

  test_gamma_matrices();

  // Define gauge field and momentum field
  field<SUN> gauge[NDIM];
  field<SUN> momentum[NDIM];

  // Define gauge and momentum action terms
  gauge_momentum_action ma(gauge, momentum);
  gauge_action ga(gauge, momentum, beta);

  ga.set_unity();

  // Define a Dirac operator
  dirac_wilson<N> D(kappa, gauge);
  fermion_action fa(D, gauge, momentum);

  // Check conjugate of the Dirac operator
  field<Wilson_vector<N>> a, b, Db, Ddaggera, DdaggerDb;
  onsites(ALL){
    a[X].gaussian();
    b[X].gaussian();
  }
  double diffre = 0, diffim = 0;
  D.apply(b, Db);
  D.dagger(a, Ddaggera);
  onsites(ALL){
    diffre += a[X].dot(Db[X]).re - Ddaggera[X].dot(b[X]).re;
    diffim += a[X].dot(Db[X]).im - Ddaggera[X].dot(b[X]).im;
  }

  assert(diffre*diffre < 1e-16 && "test Wilson dirac conjugate");
  assert(diffim*diffim < 1e-16 && "test Wilson dirac conjugate");

    double eps = 1e-6;
    for(int ng = 0; ng < ga.n_generators(); ng++){
      fa.draw_gaussian_fields();
      foralldir(dir){
        onsites(ALL){
          gauge[dir][X].random();
          momentum[dir][X] = 0;
        }
      }

      SU<N> g1 = gauge[0].get_value_at(50);
      SU<N> h = 1;
      h += eps * ga.generator(ng);
      SU<N> g12 = h*g1;

      static field<Wilson_vector<N>> psi, chi, tmp, tmp2;
      onsites(ALL){
        psi[X].gaussian();
        chi[X].gaussian();
      }
      double s1 = 0;
      D.apply(psi,tmp);
      onsites(ALL){
        s1 += chi[X].rdot(tmp[X]);
      }

      if(mynode()==0){
        gauge[0].set_value_at(g12,50);
      }
      gauge[0].mark_changed(ALL);
      double s2 = 0;
      D.apply(psi,tmp);
      onsites(ALL){
        s2 += chi[X].rdot(tmp[X]);
      }

      if(mynode()==0)
        gauge[0].set_value_at(g1, 50);
      gauge[0].mark_changed(ALL);

      D.force(chi, psi, momentum);
      SU<N> f = momentum[0].get_value_at(50);
      double diff = (f*ga.generator(ng)).trace().re - (s2-s1)/eps;

      if(mynode()==0) {
        hila::output << "Action 1 " << s1 << "\n";
        hila::output << "Action 2 " << s2 << "\n";
        hila::output << "Calculated deriv " << (f*ga.generator(ng)).trace().re << "\n";
        hila::output << "Actual deriv " << (s2-s1)/eps << "\n";
        hila::output << "deriv " << ng << " diff " << diff << "\n";
        assert( diff*diff < eps*eps*1000 && "Fermion deriv" );
      }
    }




  

  finishrun();

  return 0;
}
