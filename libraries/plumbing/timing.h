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
  std::string label;

public:  
  // initialize timer to this timepoint
  timer() {}
  timer(const char * tag) {
    init(tag);
  }

  ~timer() { remove(); }
  
  void init(const char * tag);
  void remove();
  void reset();
  double start();
  double end();
  void report();
  
};

void report_timers();

//////////////////////////////////////////////////////////////////
// Prototypes
//////////////////////////////////////////////////////////////////

double gettime();
void inittime();

bool time_to_exit();
void setup_timelimit(long seconds);

void timestamp(const char *msg);


#endif  /* timing */
