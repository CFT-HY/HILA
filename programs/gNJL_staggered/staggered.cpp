/*********************************************
 * Simulates staggered fermions with a gauge *
 * interaction                               *
 *********************************************/


#include "staggered.h"





void measure(double mass, double gamma, field<SU<N>> (&gauge)[NDIM], field<double> &sigma, field<double> &pi){
  static int iter = 0;

  output0 << " Measure_start " << iter << "\n";
  double plaq = plaquette(gauge);
  output0 << "Plaq: " << plaq << "\n";

  double sigmasq = 0, sigma_ave = 0;
  double pisq = 0, pi_ave = 0;
  onsites(ALL){
    sigma_ave += sigma[X];
    sigmasq += sigma[X]*sigma[X];
    pi_ave += pi[X];
    pisq += pi[X]*pi[X];
  }
  sigma_ave /= lattice->volume();
  sigmasq /= lattice->volume();
  pi_ave /= lattice->volume();
  pisq /= lattice->volume();
  output0 << "Sigma: " << sigma_ave << "\n";
  output0 << "Sigmasq: " << sigmasq << "\n";
  output0 << "Pi: " << pi_ave << "\n";
  output0 << "Pisq: " << pisq << "\n";

  output0 << "AUXSQ: " << sigmasq+pisq << "\n";


  double susc = -16*sigma_ave*sigma_ave;
  onsites(ALL){
    susc += sigma_ave*sigma[X];
  }
  double norm = 16*2*gamma*gamma*gamma*gamma;
  output0 << "SUSC: " << susc/norm << "\n";


  // Susceptibility: inefficient, but clearly correct method:
  coordinate_vector x,y;
  bool is_same = true;
  foralldir(dir){
    x[dir] = (int)(lattice->size(dir)*hila_random());
    y[dir] = (int)(lattice->size(dir)*hila_random());
    if(x[dir] != y[dir]) is_same = false;
  }



  dirac_staggered_gNJL<VEC, SUN> D(mass, gauge, sigma, pi);
  field<VEC> src, prop, tmp;
  CG inverse(D);

  VEC v = 0; v.c[0].re = 1;
  src[ALL] = 0;
  src.set_element(v, x);

  prop[ALL] = 0;
  D.dagger(src,tmp);
  inverse.apply(tmp, prop);

  double psibarpsi_x = prop.get_element(x).c[0].re;
  output0 << "PSIBARPSI: " << psibarpsi_x << "\n";

  if(is_same){
    output0 << "PSISUSC: 0\n";
    output0 << "PSISUSC_disc: 0\n";
    output0 << "PSISUSC_conn: 0\n";
  } else {

    D.dagger(prop,tmp);
    inverse.apply(tmp, prop);
    double susc_conn = 2*prop.get_element(x).c[0].re;

    v = 0; v.c[0].re = 1;
    src[ALL] = 0;
    src.set_element(v, y);

    prop[ALL] = 0;
    D.dagger(src,tmp);
    inverse.apply(tmp, prop);
    double psibarpsi_y = prop.get_element(y).c[0].re;
    prop[ALL] = 0;
    pi[ALL] = -pi[X];
    D.dagger(src,tmp);
    inverse.apply(tmp, prop);
    pi[ALL] = -pi[X];
    psibarpsi_y += prop.get_element(y).c[0].re;

    double susc_disc = lattice->volume()*psibarpsi_x*psibarpsi_y;

    output0 << "PSISUSC: " << susc_conn + susc_disc << "\n";
    output0 << "PSISUSC_disc: " << susc_disc << "\n";
    output0 << "PSISUSC_conn: " << susc_conn << "\n";

  }

  output0 << " Measure_end " << iter << "\n";
  iter++;
}




int main(int argc, char **argv){

  // Read parameters
  input parameters = input();
  parameters.import("parameters");
  double beta = parameters.get("beta");
  double gamma = parameters.get("gamma");
  double mass = parameters.get("mass");
  int seed = parameters.get("seed");
	double hmc_steps = parameters.get("hmc_steps");
	double traj_length = parameters.get("traj_length");
	double alpha = parameters.get("alpha");
	std::string configfile = parameters.get("configuration_file");

  lattice->setup( nd[0], nd[1], nd[2], nd[3], argc, argv );
  seed_random(seed);


  // Define gauge field and momentum field
  field<SUN> gauge[NDIM];
  field<SUN> momentum[NDIM];

  // Define gauge and momentum action terms
  gauge_momentum_action gma(gauge, momentum);
  gauge_action ga(gauge, momentum, beta);

  // Initialize the gauge field
  ga.set_unity();
  
  // Define auxiliary fields and momentum
  field<double> sigma, pi, sigma_mom, pi_mom;

  // And action
  auxiliary_momentum_action ama(sigma, pi, sigma_mom, pi_mom, alpha);
  auxiliary_action aa(sigma, pi, sigma_mom, pi_mom, gamma);

  // Initialize
  sigma[ALL] = 0.5;
  pi[ALL] = 0;

  // Define a Dirac operator (2 flavors)
  dirac_staggered_gNJL<VEC, SUN> D(mass, gauge, sigma, pi);
  gNJL_fermion_action fa(D, momentum, sigma_mom, pi_mom);

  // Build two integrator levels. Gauge is on the lowest level and
  // the fermions are on higher level
  integrator integrator_level_1(ga+aa, gma+ama);
  integrator integrator_level_2(fa, integrator_level_1);


  int file_found = (bool)std::ifstream(configfile);
  broadcast(file_found);
  if( file_found )
  {
    output0 << "Found configuration file, reading\n";
    read_fields(configfile, gauge[0], gauge[1], gauge[2], gauge[3], sigma, pi);
  } else {
    output0 << "No config file " << configfile << ", starting new run\n";
  }

  // Run HMC using the integrator
  for(int step = 0; step < 1000; step ++){
    // Run update
    update_hmc(integrator_level_2, hmc_steps, traj_length);
    measure(mass, gamma, gauge, sigma, pi);
  
    write_fields(configfile, gauge[0], gauge[1], gauge[2], gauge[3], sigma, pi);
  }



  finishrun();

  return 0;
}
