#ifndef TIMING_H
#define TIMING_H

#include "plumbing/defs.h"

namespace hila {

////////////////////////////////////////////////////////////////
/// This file defines timer class and other timing related utilities

struct timer_value {
    double time;   // time accumulated in this timer (in sec)
    int64_t count; // how many times this timer has been used (start - stop interval)
};

// clang-format off
////////////////////////////////////////////////////////////////
///
/// Timers are used to time recurring events.  Usage:
///
///         static timer loop_timer("Loop");
///
///         loop_timer.start();
///            < timed section >
///         loop_timer.stop();
///
/// All timer values are automatically reported on program exit (hila::finishrun calls
/// report_timers())
///
/// Timer can be reset with
///       loop_timer.reset();
///
/// Timer value can be inquired with
///       timer_value tval = loop_timer.value();
/// where timer_value is a struct defined below.
///
///
/// Other time related functions:
///       double gettime();                    - returns the time in seconds from program
///       start void timestamp("message");     - prints msg + current date+time +
///                                              time from start
///
///       void setup_timelimit(long seconds);  - sets the cpu time limit (see
///                                              time_to_exit()) 
///       bool time_to_exit();                 - to be called periodically
///                                              at a suitable spot in
///                                              the program; returns true if the program
///                                              should exit now.
///
////////////////////////////////////////////////////////////////
// clang-format on


class timer {
  private:
    double t_start, t_total;
    int64_t count; // need more than 32 bits
    std::string label;
    bool is_on, is_error;

  public:
    // initialize timer to this timepoint
    timer() { init(nullptr); }
    timer(const char *tag) { init(tag); }

    ~timer() { remove(); }

    void init(const char *tag);
    void remove();
    void reset();
    double start();
    double stop();
    void error();
    void report(bool print_not_timed = false);

    timer_value value();
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


} // namespace hila

#endif /* timing */
