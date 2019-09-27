#ifndef TIMERS_H
#define TIMERS_H

#include "globals.h"
#include "clock.h"

class timer {
private:
  double t_start, t_total, t_initial;
  unsigned long long count;   // need more than 32 bits

public:  
  // initialize timer to this timepoint
  timer() { 
    reset();
  }
  ~timer() {}
  
  void reset();
  double start();
  double end();
  void report(const char * label);
  
};

#endif  /* timers */
