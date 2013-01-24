
#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>


double t_init = 0.0;

double
timer ()
{
  struct timeval tv;
  struct timezone tz;
  double t;

  gettimeofday (&tv, &tz);
  t = tv.tv_sec + 1.e-6 * tv.tv_usec;
  if (t_init == 0.0)
    {
      t_init = t;
      return (0.0);
    }
  else
    {
      return (t - t_init);
    }
}
