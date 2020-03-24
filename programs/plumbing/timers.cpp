
#include "globals.h"
#include "timers.h"
#include "clock.h"

  // initialize timer to this timepoint
void timer::reset() {
  t_start = t_total = 0.0;
  t_initial = gettime();
  count = 0;
}

double timer::start() {
  if (mynode() == 0) {
    t_start = gettime();
    return t_start;
  } else return 0.0;
}
  
double timer::end() {
  if (mynode() == 0) {
    double e = gettime();
    t_total += (e - t_start);
    count++;
    return e;
  } else return 0.0;
}

void timer::report(const char * label, int print_header) {
  if (mynode() == 0) {
    static bool first = true;
    char line[200];

    if (print_header > 0 || (first && print_header == -1)) {
      first = false;
      hila::output << "TIMER REPORT:         total(sec)          calls   usec/call  fraction\n";
    }
    // time used during the counter activity
    t_initial = gettime() - t_initial;
    if (count > 0) {
      std::snprintf(line,200," %16s: %14.3f %14llu %11.3f %8.4f\n",
                    label, t_total, count, 1e6 * t_total/count, t_total/t_initial );
    } else {
      std::snprintf(line,200," %16s: no timed calls made\n",label);
    }
    hila::output << line;      
  }
}
  
