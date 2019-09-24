#ifndef CLOCK_H
#define CLOCK_H
//**************************************************
//  Fast clock independent of MPI
//
// #include <sys/time.h>
#include <time.h>

#if !defined(BGL)

// clock_gettime() is the clock for modern Linux - TODO: check architectures!

inline double gettime()
{
  struct timespec tp;

  clock_gettime(CLOCK_MONOTONIC, &tp);
  return( (double)tp.tv_sec + 1.0e-9*(double)resource.tv_nsec );
}

#else

/* For BLUE GENE/L, timer is rts_get_timebase(),
 * which returns the cycles since boot
 *
 * This does not work in BG/P! Use gettimeofday() above
 *
 */

// #include <rts.h>
// 
// inline double time()
// {
//   static double clockspeed=1.0e-6/700.0;
// 
//   return ( rts_get_timebase() * clockspeed );
// }


#endif

#endif // CLOCK_H