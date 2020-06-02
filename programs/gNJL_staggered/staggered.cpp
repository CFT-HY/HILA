/*********************************************
 * Simulates staggered fermions with a gauge *
 * interaction                               *
 *********************************************/


#include "staggered.h"



double polyakov(direction dir, field<SU<N>> (&gauge)[NDIM]){
  coordinate_vector vol = lattice->size();
  field<SU<N>> polyakov; polyakov[ALL] = 1;
  for(int t=0; t<vol[dir]; t++){
    onsites(ALL){
      if(X.coordinates()[dir]==(t+1)%vol[dir]){
        polyakov[X] = polyakov[X] * gauge[dir][X-dir];
      }
    }
  }

  double poly=0;
  onsites(ALL){
    if(X.coordinates()[dir]==0){
      poly += polyakov[X].trace().re;
    }
  }

  double v3 = lattice->volume()/vol[dir];
  return poly/(N*v3);
}



void measure(double mass, double gamma, field<SU<N>> (&gauge)[NDIM], field<double> &sigma, field<double> &pi){
  static int iter = 0;
  coordinate_vector vol = lattice->size();

  output0 << " Measure_start " << iter << "\n";
  double plaq = plaquette(gauge);
  output0 << "Plaq: " << plaq << "\n";

  double poly = polyakov(TUP, gauge);
  output0 << "POLYAKOV_T: " << poly << "\n";
  poly = polyakov(XUP, gauge);
  output0 << "POLYAKOV_X: " << poly << "\n";
  poly = polyakov(YUP, gauge);
  output0 << "POLYAKOV_Y: " << poly << "\n";
  poly = polyakov(ZUP, gauge);
  output0 << "POLYAKOV_Z: " << poly << "\n";
  

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
    x[dir] = (int)(vol[dir]*hila_random());
    y[dir] = (int)(vol[dir]*hila_random());
    if(x[dir] != y[dir]) is_same = false;
  }

  std::vector<double> sigma_prop(vol[TUP]);
  std::fill(sigma_prop.begin(), sigma_prop.end(), 0);
  coordinate_vector c(0);
  for(int t=0; t<vol[TUP]; t++){
    std::vector<double> sigma_prop_t(vol[TUP]);
    std::fill(sigma_prop_t.begin(), sigma_prop_t.end(), 0);
    c[TUP] = t;
    double sigma_0 = sigma.get_element(c);
    broadcast(sigma_0);
    onsites(ALL){
      coordinate_vector l = X.coordinates();
      int t = l[TUP];

      sigma_prop_t[t] += sigma_0 * sigma[X];
    }
    for(int t2=0; t2<vol[TUP]; t2++){
      sigma_prop[t2] += sigma_prop_t[t2];
    }
  }
  for(int t=0; t<vol[TUP]; t+=2){
    output0 << "SIGMA_PROP: " << t << " " << sigma_prop[t] / lattice->volume() << "\n";
  }

  std::vector<double> pi_prop(vol[TUP]);
  std::fill(pi_prop.begin(), pi_prop.end(), 0);
  c = 0;
  for(int t=0; t<vol[TUP]; t++){
    std::vector<double> pi_prop_t(vol[TUP]);
    std::fill(pi_prop_t.begin(), pi_prop_t.end(), 0);
    c[TUP] = t;
    double sigma_0 = sigma.get_element(c);
    broadcast(sigma_0);
    onsites(ALL){
      coordinate_vector l = X.coordinates();
      int t = l[TUP];

      pi_prop_t[t] += sigma_0 * sigma[X];
    }
    for(int t2=0; t2<vol[TUP]; t2++){
      pi_prop[t2] += pi_prop_t[t2];
    }
  }
  for(int t=0; t<vol[TUP]; t+=2){
    output0 << "PI_PROP: " << t << " " << pi_prop[t] / lattice->volume() << "\n";
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
	int n_trajectories = parameters.get("n_trajectories");
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
  for(int step = 0; step < n_trajectories; step ++){
    // Run update
    update_hmc(integrator_level_2, hmc_steps, traj_length);
    measure(mass, gamma, gauge, sigma, pi);
  
    write_fields(configfile, gauge[0], gauge[1], gauge[2], gauge[3], sigma, pi);
  }



  finishrun();

  return 0;
}
