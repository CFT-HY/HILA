#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <assert.h>
#include <sys/time.h>
#include <ctime>

#include "plumbing/defs.h"
#include "datatypes/matrix.h"
#include "datatypes/sun.h"
//#include "datatypes/wilson_vector.h"
#include "plumbing/field.h"
//#include "dirac/staggered.h"
//#include "dirac/wilson.h"

// Minimum time to run each benchmark
// in microseconds
constexpr double mintime = 1000;

// Direct output to stdout
// std::ostream &hila::output = std::cout;

// Calculate time difference in milliseconds
static inline double timediff(timeval start, timeval end) {
    long long t1 = (long long)(start.tv_usec) + 1000000 * (long long)(start).tv_sec;
    long long t2 = (long long)(end.tv_usec) + 1000000 * (long long)(end).tv_sec;
    return 1e-3 * (double)(t2 - t1);
}
