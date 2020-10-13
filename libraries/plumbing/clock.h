#ifndef CLOCK_H
#define CLOCK_H
//**************************************************
//  Fast clock independent of MPI
//
// #include <sys/time.h>
#include <time.h>

#if !defined(USE_MPI)

// clock_gettime() is the clock for modern Linux - TODO: check architectures!

inline double gettime() {
  struct timespec tp;

  clock_gettime(CLOCK_MONOTONIC, &tp);
  return( (double)tp.tv_sec + 1.0e-9*(double)tp.tv_nsec );
}

#else

// Use MPI_Wtime() when MPI is available,

#include <mpi.h>

inline double gettime() { 
  return MPI_Wtime();
}

#endif

#endif // CLOCK_H