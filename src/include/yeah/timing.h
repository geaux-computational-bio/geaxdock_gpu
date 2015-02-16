//#ifndef _TIMING_H_
//#define _TIMING_H_


#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>




typedef struct
{
  double start; // start time stamp
  double span; // accumulated time span
  double avrg; // time span average
  int n; // number of records
} Timing;




double
HostTimeNow ()
{
  struct timeval t;
  gettimeofday (&t, NULL);
  double mytime_second = (double) t.tv_sec + (double) t.tv_usec * 1e-6;
  return mytime_second;
}


// fail ?
/*
double
HostTimeNow2 ()
{
  struct timespec tt; 
  clock_gettime (CLOCK_REALTIME, &tt);
  double mytime_second = (double) tt.tv_sec + (double) tt.tv_nsec * 1e-9;
  return mytime_second;
}
*/




void
HostTimingInit (Timing * t, int n)
{
  for (int i = 0; i < n; i++) {
    t[i].start = 0.0;
    t[i].span = 0.0;
    t[i].n = 0;
  }
}

void
HostTiming (Timing * t, int idx, int mode)
{
  switch (mode) {
  case 0: // marking start time stamp
    t[idx].start = HostTimeNow ();
    break;
  case 1: // marking end time stamp, and accumulated time span
    t[idx].span += t[idx].start - HostTimeNow ();
    t[idx].n += 1;
    break;
  default: // culculate average. this action can be triggered only once
    t[idx].avrg = t[idx].span / t[idx].n;
  }
}







//#endif

