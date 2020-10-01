#ifndef HMC_H
#define HMC_H


#include <sys/time.h>
#include <ctime>
#include "integrator.h"



/// The Hybrid Montecarlo algorithm.
// Consists of an integration step following equations of
// motion implemented in the integrator class gt
// and an accept-reject step using the action
//
// The integrator class must implement at least two functions, 
// action() an integrator_step(double eps) 
template<class integrator_type>
void update_hmc(integrator_type &integrator, int steps, double traj_length){
  
  static int accepted=0, trajectory=1;
  struct timeval start, end;
  double timing;

  // Draw the momentum
  integrator.draw_gaussian_fields();

  // Make a copy of the gauge field in case the update is rejected
  integrator.backup_fields();
  
  gettimeofday(&start, NULL);

  // Calculate the starting action and print
  double start_action = integrator.action();
  output0 << "Begin HMC Trajectory " << trajectory << ": Action " 
          << start_action << "\n";

  // Run the integator
  for(int step=0; step < steps; step++){
    integrator.step(traj_length/steps);
  }

  // Recalculate the action
  double end_action = integrator.action();
  double edS = exp(-(end_action - start_action));

  output0 << "End HMC: Action " << end_action << " "
        << end_action - start_action
        << " exp(-dS) " << edS << "\n";

  // Accept or reject
  if(hila_random() < edS){
    output0 << "Accepted!\n";
    accepted++;
  } else {
    output0 << "Rejected!\n";
    integrator.restore_backup();
  }


  output0 << "Acceptance " << accepted << "/" << trajectory 
          << " " << (double)accepted/(double)trajectory << "\n";

  gettimeofday(&end, NULL);
  timing = (double)(end.tv_sec - start.tv_sec) + 1e-6*(end.tv_usec - start.tv_usec);

  output0 << "HMC done in " << timing << " seconds \n";
  trajectory++;
}



#endif
