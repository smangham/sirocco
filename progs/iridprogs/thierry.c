
/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:

  Description:	These routines are intended to read the ?? star models, which
	Ivan and Thierry gave me.

  Arguments:		

  Returns:

  Notes:
	r proceeds from the outside of the star in

  History:
	01aug	ksl	Coded
	

 ************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"

#define LINELENGTH 160

#include "python.h"

//Internal structure for use with this stellar gridded model
#define NDEPTHS 100
struct xxstar
{
  int nd;			// Number of depth points
  double r[NDEPTHS];		//distance to the center of star (outside in_order)
  double depth[NDEPTHS];	//distance from outer edge of model
  double v[NDEPTHS];
  double v_offset[NDEPTHS];	//kluge for numerical recipes
  double v_spline2[NDEPTHS];	//calculated spline coefficients
  double rho[NDEPTHS];
  double rho_offset[NDEPTHS];	//kluge for numerical recipes
  double rho_spline2[NDEPTHS];	//calculated spline coefficients
  double t[NDEPTHS];
}
xstar;



/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	get_thierry_wind_params gets input data which is necessary for Ivan and Thierry's
	wind models
Arguments:		

Returns:
 
Description:	
Notes:
History:
	01sept	ksl	Coded.
**************************************************************/


int
get_thierry_params ()
{
  char infile[LINELENGTH];
  int read_thierry ();

  Log ("Creating a wind model from Thierry's code\n");

  geo.wind_thetamin = 0.0;
  geo.wind_rmax = geo.rmax;
  geo.wind_thetamax = 90. / RADIAN;
  geo.wind_rmin = geo.rstar;

  geo.diskrad = 0.0;
  geo.disk_type = 0;

  strcpy (infile, "models/RVT30");

  rdstr ("Filename_for_thierry_wind", infile);

  read_thierry (infile);

  geo.wind_rmin = geo.rstar = xstar.r[xstar.nd - 1];
  geo.wind_rmax = geo.rmax = xstar.r[0];
  geo.wind_rho_min = 0;
  geo.wind_rho_max = geo.rmax = xstar.r[0];
  geo.wind_thetamin = 0.0;
  geo.wind_thetamax = 90. / RADIAN;

  geo.xlog_scale = 0.3 * geo.rstar;
  geo.zlog_scale = 0.3 * geo.rstar;

/* Moved the wind_cone definition back to main in 01dec ksl */

  return (0);
}

int
read_thierry (infile)
     char infile[];
{
  char line[200], firstword[80];
  double array[NDEPTHS];
  int i, mchar;
  FILE *fopen (), *fptr;

/* Initialize the xstar structure */
  for (i = 0; i < NDEPTHS; i++)
    {
      xstar.r[i] = 0.0;
      xstar.v[i] = 0.0;
      xstar.rho[i] = 0.0;
      xstar.t[i] = 0.0;
    }


/* Open the input file.  Exit if it is not opened successfully */
  if ((fptr = fopen (infile, "r")) == NULL)
    {
      printf ("Failed to open file %s\n", infile);
      exit (0);
    }

/* Read and print the input file */
  while (fgets (line, 160, fptr) != NULL)
    {
      sscanf (line, "%s", firstword);
      mchar = strlen (line);
      printf ("%s totol chars %d\n", firstword, mchar);

      if (strncmp (firstword, "ND:", 3) == 0)
	{
	  printf ("Got ND points\n");
	  sscanf (line, "%*s %d", &xstar.nd);

	}
      if (strncmp (firstword, "Radius", 3) == 0)
	{
	  printf ("Got radius\n");
	  thierry_read_set (fptr, array, xstar.nd);
	  for (i = 0; i < xstar.nd; i++)
	    xstar.r[i] = 1e10 * array[i];	// Convert to cm
	}

      if (strncmp (firstword, "Velocity", 3) == 0)
	{
	  printf ("Got vel\n");
	  thierry_read_set (fptr, array, xstar.nd);
	  for (i = 0; i < xstar.nd; i++)
	    xstar.v[i] = 1e5 * array[i];	// Convert to cm/s
	}
      if (strncmp (firstword, "Mass", 3) == 0)
	{
	  printf ("Got rho\n");
	  thierry_read_set (fptr, array, xstar.nd);
	  for (i = 0; i < xstar.nd; i++)
	    xstar.rho[i] = array[i];
	}
      if (strncmp (firstword, "Temperature", 3) == 0)
	{
	  printf ("Got t\n");
	  thierry_read_set (fptr, array, xstar.nd);
	  for (i = 0; i < xstar.nd; i++)
	    xstar.t[i] = 1e4 * array[i];	// Convert to degrees
	}

    }

  for (i = 0; i < xstar.nd; i++)
    {
      xstar.depth[i + 1] = xstar.r[0] - xstar.r[i];
      xstar.v_offset[i + 1] = xstar.v[i];
    }

//spline(&xstar.depth,&xstar.v_offset,xstar.nd,1.e30,1.e30,&xstar.v_spline2);
  spline (&xstar.depth[0], &xstar.v_offset[0], xstar.nd, 1.e30, 1.e30,
	  &xstar.v_spline2[0]);
  spline (&xstar.depth[0], &xstar.rho_offset[0], xstar.nd, 1.e30, 1.e30,
	  &xstar.rho_spline2[0]);

  return (0);
}

int
thierry_read_set (fptr, array, nd)
     FILE *fptr;
     double array[];
     int nd;			// number of depth points
{
  int i, j, ncards;
  char line[132];
  double a[8];

  ncards = (nd + 7) / 8;
  for (i = 0; i < nd; i++)
    {
      if ((j = i % 8) == 0)
	{
	  fgets (line, 132, fptr);
	  printf ("%s\n", line);
	  sscanf (line, "%le %le %le %le %le %le %le %le", &a[0], &a[1],
		  &a[2], &a[3], &a[4], &a[5], &a[6], &a[7]);
	}
      array[i] = a[i % 8];
    }

  return (0);
}

/*
Note that for these models the velocities are specifified from the outside in in order to avoid having
to reorder the input matrices.
*/
double
xthierry_velocity (x, v)
     double x[];
     double v[];
{
  double r, speed;
  double length ();

  r = length (x);

  if (sane_check (r))
    {
      Error ("xthierry_velocity: sane_check failure\n");
    }

  if (r == 0)
    {
      v[0] = v[1] = v[2] = speed = 0.0;
      return (speed);
    }

  if (r < xstar.r[xstar.nd - 1])
    {
      speed = xstar.v[xstar.nd - 1] * r / xstar.r[xstar.nd - 1];
    }
  else if (r > xstar.r[0])
    {
      speed = xstar.v[0];
    }

  else
    {
//spline(&xstar.depth,&xstar.v_offset,xstar.nd,1.e30,1.e30,&xstar.v_spline2);

      sane_check (speed);
//splint(&xstar.depth,&xstar.v_offset,&xstar.v_spline2,xstar.nd,xstar.r[0]-r,&speed);
      splint (&xstar.depth[0], &xstar.v_offset[0], &xstar.v_spline2[0],
	      xstar.nd, xstar.r[0] - r, &speed);
      sane_check (speed);
    }


  sane_check (r);

  v[0] = speed * x[0] / r;
  v[1] = speed * x[1] / r;
  v[2] = speed * x[2] / r;


  if (sane_check (speed) || sane_check (v[0]) || sane_check (v[1])
      || sane_check (v[2]))
    {
      Error ("thierry_velocity: x %f %f %f speed %f r %f\n", x[0], x[1],
	     x[2], speed, r);
    }
  return (speed);
}

double
thierry_velocity (x, v)
     double x[];
     double v[];
{
  int n;
  double r, speed;
  double frac, z;
  double length ();

  n = 0;

  z = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
  if (sane_check (z))
    {
      printf ("%f\n", z);
      z = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
      printf ("%f\n", z);
    }
  r = length (x);

  if (sane_check (r))
    {
      printf ("%f\n", r);
      r = length (x);
      printf ("%f\n", r);
    };

  if (r == 0)
    {
      v[0] = v[1] = v[2] = speed = 0;
    }
  else
    {
      if (r > xstar.r[0])
	{
	  speed = xstar.v[0];
	}

      else
	{
//These models are organized from outside in 
	  n = 0;
	  while (r < xstar.r[n] && n < xstar.nd)
	    n++;
/*so xstar[n].r is just less than r in
which case xstar[n-1] is just greater than r
or n=xstar.nd); */

	  if (n == xstar.nd)
	    speed = xstar.v[(xstar.nd) - 1];
	  else if (n == 0)
	    speed = xstar.v[0];
	  else
	    {
	      frac = (r - xstar.r[n]) / (xstar.r[n - 1] - xstar.r[n]);
	      speed = (1. - frac) * xstar.v[n] + (frac) * xstar.v[n - 1];

	      sane_check (r);
	    }
	}


      sane_check (r);

      v[0] = speed * x[0] / r;
      v[1] = speed * x[1] / r;
      v[2] = speed * x[2] / r;
    }


  if (sane_check (speed) || sane_check (v[0]) || sane_check (v[1])
      || sane_check (v[2]))
    {
      Error ("thierry_velocity: x %f %f %f speed %f n %d r %f\n", x[0], x[1],
	     x[2], speed, n, r);
    }
  return (speed);
}

/*
thierry_vel_grad calculates the velocity radient tensor at any point in
the flow

The velocity gradient is defined as a 3 x 3 tensor such that

velgrad[i][j]= dv_i/dx_j

NB in c the rightmost index changes changes most rapidly so that
dv[i]= velgrad[i][j] * dx[j]
makes sense.

NB: Making ds too small can cause roundoff and/or precision errors.

        01dec   ksl     Added for python_40

*/
int
thierry_vel_grad (x, velgrad)
     double x[], velgrad[][3];
{
  double v0[3], v1[3];
  double dx[3], dv[3];
  double ds;
  int i, j;
  int vsub (), stuff_v ();
  thierry_velocity (x, v0);

  if (sane_check (v0[0]) || sane_check (v0[1]) || sane_check (v0[2]))
    {
      Error ("thiery_vel_grad: x %f %f %f v0 %f %f %f\n", x[0], x[1], x[2],
	     v0[0], v0[1], v0[2]);
    }

  ds = 1.e10;
  for (i = 0; i < 3; i++)
    {
      stuff_v (x, dx);
      dx[i] += ds;
      thierry_velocity (dx, v1);
      if (sane_check (v1[0]) || sane_check (v1[1]) || sane_check (v1[2]))
	{
	  Error ("sv_vel: dx %f %f %f v0 %f %f %f\n", dx[0], dx[1], dx[2],
		 v1[0], v1[1], v1[2]);
	}

      vsub (v1, v0, dv);
      for (j = 0; j < 3; j++)
	dv[j] /= ds;
      stuff_v (dv, &velgrad[i][0]);
    }

  return (0);
}



double
thierry_rho (x)
     double x[];
{

  int n;
  double r, rho, frac;
  double length ();
  r = length (x);
  n = xstar.nd - 1;
  while (r > xstar.r[n] && n >= 0)
    n--;
  if (n == xstar.nd - 1)
    rho = xstar.rho[n];
  else if (n < 0)
    rho = xstar.rho[0];
  else
    {
      frac = (r - xstar.r[n]) / (xstar.r[n + 1] - xstar.r[n]);
      rho = (1. - frac) * xstar.rho[n] + (frac) * xstar.rho[n + 1];
    }

  return (rho);
}

double
xthierry_rho (x)
     double x[];
{

  double r, rho;
  double length ();
  r = length (x);

  if (r < xstar.r[xstar.nd - 1])
    {
      rho = xstar.rho[xstar.nd - 1];
    }
  else if (r > xstar.r[0])
    {
      rho = 0.0;
    }

  else
    {
//spline(&xstar.depth,&xstar.rho_offset,xstar.nd,1.e30,1.e30,&xstar.rho_spline2);

//splint(&xstar.depth,&xstar.v_offset,&xstar.v_spline2,xstar.nd,xstar.r[0]-r,&speed);
      splint (&xstar.depth[0], &xstar.rho_offset[0], &xstar.rho_spline2[0],
	      xstar.nd, xstar.r[0] - r, &rho);
      sane_check (rho);
    }



  return (rho);
}
