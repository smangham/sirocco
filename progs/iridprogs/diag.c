




#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

int eplinit = 0;

int
open_diagfile ()
{

  if (eplinit == 0)
    {
      epltptr = fopen ("python.ext", "w");
      eplinit = 1;
    }
//  diag_on_off = 0;            // 0=off everything else is on
  return (0);
}
