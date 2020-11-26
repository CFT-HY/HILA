#ifndef TIMING_H
#define TIMING_H

#include "plumbing/defs.h"


////////////////////////////////////////////////////////////////
// Timer class
////////////////////////////////////////////////////////////////

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
  void report(const char * label, int print_header = -1);
  
};

//////////////////////////////////////////////////////////////////
// Prototypes
//////////////////////////////////////////////////////////////////

double gettime();
void inittime();

bool time_to_exit();
void setup_timelimit(long seconds);

void timestamp(const char *msg);


#endif  /* timing */
