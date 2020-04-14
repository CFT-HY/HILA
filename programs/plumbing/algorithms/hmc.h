#ifndef INTEGRATOR_H
#define INTEGRATOR_H




template<class action_term>
void leapfrog_step(action_term gauge, double eps){
  gauge.momentum_step(0.5*eps);
  gauge.force_step(eps);
  gauge.momentum_step(0.5*eps);
}


template<class action_term>
void O2_step(action_term gauge, double eps){
  double zeta = eps*0.1931833275037836;
  double middlestep = eps-2*zeta;
  gauge.force_step(zeta);
  gauge.momentum_step(0.5*eps);
  gauge.force_step(middlestep);
  gauge.momentum_step(0.5*eps);
  gauge.force_step(zeta);
}


template<class action_term>
void update_hmc(action_term gt, int steps, double traj_length){
  static int accepted=0, trajectory=1;

  // Draw the momentum
  gt.gaussian_momentum();

  // Make a copy of the gauge field in case the update is rejected
  field<SUN> gauge_copy[NDIM];
  foralldir(dir) gauge_copy[dir] = gt.gauge[dir];

  // Calculate the starting action and print
  double S_gauge = gt.gauge_action();
  double S_mom = gt.momentum_action();
  double start_action = S_gauge + S_mom;
  output0 << "Begin HMC Trajectory " << trajectory << ": Action " 
          << S_gauge << " " << S_mom << " " << start_action << "\n";



  // Run the integator
  for(int step=0; step < steps; step++){
    O2_step(gt, traj_length/steps);
  }


  // Recalculate the action
  double S2_gauge = gt.gauge_action();
  double S2_mom = gt.momentum_action();
  double edS = exp(-(S2_gauge+S2_mom - start_action));

  output0 << "End HMC: Action " << S2_gauge << " " << S2_mom << " "
        << S2_gauge + S2_mom  << " " << S2_gauge+S2_mom - start_action
        << " exp(-dS) " << edS << "\n";


  // Accept or reject
  if(hila_random() < edS){
    output0 << "Accepted!\n";
    accepted++;
  } else {
    output0 << "Rejected!\n";
  }

  output0 << "Acceptance " << accepted << "/" << trajectory 
          << " " << (double)accepted/(double)trajectory << "\n";
  trajectory++;
}



#endif
