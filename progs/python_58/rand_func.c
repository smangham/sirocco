


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  
	double gen_rand (func, xmin, xmax) returns a raandom number
	between xmin and xmax

  Description:	

  Arguments:  


  Returns:

  Notes:

 

  History:
	06oct	ksl	Coded to replace or supplement various routines
			for generating random numbers when the probability
			distribution can be expressed as a function

 ************************************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


#include "log.h"

#define LINELENGTH 132
#define MAX	1000
#define MAXRAND 2147486748.

double gen_rand_value;
double gen_rand_xmin;
double gen_rand_xmax;
double (*gen_rand_func) ();

double
gen_rand (func, xmin, xmax)
     double (*func) (double);
     double xmin, xmax;
{
  double total;
  double x;
  double qromb (), zbrent (), xgen_rand_func ();

  total = qromb (func, xmin, xmax, 1.e-6);
  gen_rand_value = rand () / MAXRAND * total;


  gen_rand_func = func;
  gen_rand_xmin = xmin;
  gen_rand_xmax = xmax;


  x = zbrent (xgen_rand_func, xmin, xmax, 1.e-4);
//  printf("foo %f %f %f %f %f\n",total,gen_rand_value,xmin,xmax,x);
  return (x);


}


double
xgen_rand_func (x)
     double x;

{
  double z, total;
  double qromb ();
  if (gen_rand_xmin == x)
    return (-gen_rand_value);
  total = qromb (gen_rand_func, gen_rand_xmin, x, 1.e-4);
  z = total - gen_rand_value;
//  printf("%f %f\n",x,z);
  return (z);
}
