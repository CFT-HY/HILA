
#include <time.h>
#include <chrono>

#include "timing.h"

//////////////////////////////////////////////////////////////////
/// Time related routines (runtime - timing - timelimit)
/// Check timing.h for details
//////////////////////////////////////////////////////////////////

#ifdef USE_MPI
#include "com_mpi.h"
#endif


// This stores the start time of the program
static double start_time = -1.0;


/////////////////////////////////////////////////////////////////
/// Timer routines - for high-resolution event timing
/////////////////////////////////////////////////////////////////

// store all timers in use
std::vector<timer *> timer_list = {};

// initialize timer to this timepoint
void timer::init(const char * tag) {
  if (tag != nullptr) label = tag;
  reset();
  // has this been inited already?  Check
  bool found = false;
  for (auto p : timer_list) if (p == this) {
    found = true;
    break;
  }
  // if not, push it to the list
  if (!found) timer_list.push_back(this);
}

// remove the timer also from the list
void timer::remove() {
  for (auto it = timer_list.begin(); it != timer_list.end(); ++it) {
    if (*it == this) {
      timer_list.erase(it);
      return;
    }
  }
}

void timer::reset() {
  t_start = t_total = 0.0;
  count = 0;
}

double timer::start() {
  t_start = gettime();
  return t_start;
}
  
double timer::stop() {
  double e = gettime();
  t_total += (e - t_start);
  count++;
  return e;
}

timer_value timer::value() {
  timer_value r;
  r.time = t_total;
  r.count = count;
  return r;
}

void timer::report(bool print_not_timed ) {
  if (hila::myrank() == 0) {
    char line[202];

    // time used during the counter activity
    double ttime = gettime();
    if (count > 0) {
      if (t_total/count < 0.01) {
        std::snprintf(line,200,"%-20s: %14.3f %14llu %10.3f us %8.4f\n",
                      label.c_str(), t_total, count, 1e6 * t_total/count, t_total/ttime );
      } else {
        std::snprintf(line,200,"%-20s: %14.3f %14llu %10.3f s  %8.4f\n",
                      label.c_str(), t_total, count, t_total/count, t_total/ttime );
      }
      hila::output << line;
    } else if (print_not_timed) {
      std::snprintf(line,200,"%-20s: no timed calls made\n",label.c_str());
      hila::output << line;      
    }
  }
}

void report_timers() {
  if (hila::myrank() == 0) {
    if (timer_list.size() > 0) {
      hila::output << "TIMER REPORT:             total(sec)          calls     time/call  fraction\n";
      hila::output << "---------------------------------------------------------------------------\n";

      for (auto tp : timer_list) tp->report();

      hila::output << "---------------------------------------------------------------------------\n";
    } else {
      hila::output << "No timers defined\n";
    }
  }
}


/////////////////////////////////////////////////////////////////
/// Use clock_gettime() to get the accurate time
/// (alternative: use gettimeofday()  or MPI_Wtime())
/// gettime returns the time in secs since program start


double gettime() {
  struct timespec tp;

  if (start_time == -1.0) inittime();

  clock_gettime(CLOCK_MONOTONIC, &tp);
  return( ((double)tp.tv_sec - start_time) + 1.0e-9*(double)tp.tv_nsec );
}

void inittime() {
  if (start_time == -1.0) {
    start_time = 0.0;
    start_time = gettime();
  }
}
 
//////////////////////////////////////////////////////////////////
/// Routines for checking remaining cpu-time
/// time_to_exit() is called periodically on a point where exit can be done.
/// It uses the max of the time intervals for the estimate for one further round.
/// If not enough time returns true, else false
/// 


static double timelimit = 0;

void setup_timelimit(long seconds) {
  timelimit = (double) seconds;
}

bool time_to_exit() {
  static double max_interval = 0.0;
  static double previous_time = 0.0;
  bool finish;

  // if no time limit set
  if (timelimit == 0.0) return false;

  if (hila::myrank() == 0) {
    double this_time = gettime();
    if (this_time - previous_time > max_interval) 
      max_interval = this_time - previous_time;
    previous_time = this_time;

    // Give 5 min margin for the exit - perhaps needed for writing etc.
    if (timelimit - this_time < max_interval + 5*60.0) finish = true;
    else finish = false;

    hila::output << "TIMECHECK: " << this_time << "s used, " 
                 << timelimit - this_time << "s remaining\n";
    if (finish)
      hila::output << "CPU TIME LIMIT, EXITING THE PROGRAM\n";
  }
  broadcast(finish);
  return finish;
}


/*****************************************************
 * Time stamp
 */

void timestamp(const char *msg)
{
  if (hila::myrank() == 0) {
    std::time_t ct = std::time(NULL);
    if (msg != NULL) hila::output << msg;
    std::string d = ctime(&ct);
    d.resize(d.size()-1);   // take away \n at the end
    hila::output << " -- date " << d << "  run time " << gettime() << "s\n";
  }
}

