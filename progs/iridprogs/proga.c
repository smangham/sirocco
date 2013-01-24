
/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  These routines are designed to match Daniel Proga's wind
	models onto sv

  Description:	These routines are designed to read the velocity and densities
	of Daniel Proga's wind models into the appropriate structures and
	allow one to calculate the density at any point in the gridded space.
	

  Arguments:		

  Returns:

  Notes:
	Daniel's models data are the results of calculations with Zeus.  The
	grid is in polar (r,theta) coordinates with theta = pi/2 corresponding
	to the disk, and 0 to the pole.  Recall that python uses cylintrical
	corrdinates with the z=0 being the disk plane.

  History:
	99dec	ksl	Began work
	00jan	ksl	Modified to read larger models and to make it slightly
			more general in terms of the input files.  Also added
			some checking on grid sizes, and set rho to zero outside
			the range of the model, and velocity at the edge of
			the model.
	04jun	ksl	Moved get_proga_wind_params from python.c to this file
	

 ************************************************************************/
#include	<stdio.h>
#include	<stdlib.h>
#include	<strings.h>
#include	<string.h>
#include	<math.h>
#include 	"log.h"
#include	"atomic.h"
#include 	"python.h"

#define LINE 132




/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	get_proga_wind_params gets input data for Daniel Proga's wind models
Arguments:		

Returns:
 
Description:	
Notes:
History:
 	99dec	ksl	Began work
**************************************************************/


int
get_proga_wind_params ()
{
  int get_proga ();





  Log
    ("Creating a wind model for a Cataclysmic Variable using a Proga Calc\n");

  get_proga ();
  geo.wind_rmin = 8.7e8;	/*Radius where wind begins */
  if (geo.wind_rmin < geo.rstar)
    {
      Error
	("get_stellar_wind_params: It is unreasonable to have the wind start inside the star!\n");
      Log ("Setting geo.wind_rmin to geo.rstar\n");
      geo.wind_rmin = geo.rstar;
    }

/* Assign the generic parameters for the wind the generic parameters of the wind */

  geo.wind_rmin = geo.rstar;
  geo.wind_rmax = geo.rmax;
  geo.wind_rho_min = 0;
  geo.wind_rho_max = geo.rmax;
  geo.wind_thetamin = 0.0;
  geo.wind_thetamax = 90. / RADIAN;

//      geo.wind_rho_min=4*8.7e8;
  //      geo.wind_rho_max=12.*8.7e8;
  //      geo.wind_thetamin=20.0/RADIAN;
  //      geo.wind_thetamax=65./RADIAN;

  geo.xlog_scale = 0.3 * geo.rstar;
  geo.zlog_scale = 0.3 * geo.rstar;

  return (0);
}

#define MAXPROGA 155


double proga_r[MAXPROGA];
double proga_theta[MAXPROGA];
int iproga_r, iproga_theta, iproga_mod;

typedef struct proga_mod
{
  double rho;
  double v[3];
}
proga_dummy, *ProgaPtr;

ProgaPtr proga_ptr;


int
get_proga ()
{

  FILE *fopen (), *fptr;
  char rfile[LINE], thetafile[LINE], datafile[LINE];
  char aline[LINE];
  int i, j, k;
  double r;
  double rho, theta, vr, vtheta, vphi;

/*Write something into the file name strings */
  strcpy (rfile, "grid1_r_big.dat");
  strcpy (thetafile, "grid_theta.dat");
  strcpy (datafile, "hdf062.221.dat");


  for (k = 0; k < MAXPROGA; k++)
    {
      proga_r[k] = 0;
      proga_theta[k] = 0;
    }

  proga_ptr = (ProgaPtr) calloc (sizeof (proga_dummy), MAXPROGA * MAXPROGA);

  if (proga_ptr == NULL)
    {
      printf
	("There is a problem in allocating memory for the proga structure\n");
      exit (0);
    }

  rdstr ("Proga_radii_file", rfile);
  if ((fptr = fopen (rfile, "r")) == NULL)
    {
      printf ("Could not open %s\n", rfile);
      exit (0);
    }

  iproga_r = 0;
  while (fgets (aline, LINE, fptr) != NULL)
    {
      sscanf (aline, "%d %lf", &i, &r);
      proga_r[i] = r;
      if (i > iproga_r)
	iproga_r = i;
    }

  iproga_r = i;
  printf ("Read %d values from grid_r.dat\n", i);
  fclose (fptr);

  rdstr ("Proga_theta_file", thetafile);
  if ((fptr = fopen (thetafile, "r")) == NULL)
    {
      printf ("Could not open %s\n", thetafile);
      exit (0);
    }

  i = 0;
  while (fgets (aline, LINE, fptr) != NULL)
    {
      sscanf (aline, "%d %lf", &i, &theta);
      proga_theta[i] = theta;
      i++;
    }

  iproga_theta = i;
  printf ("Read %d values from %s\n", i, thetafile);
  fclose (fptr);


  rdstr ("Proga_data_file", datafile);
  if ((fptr = fopen (datafile, "r")) == NULL)
    {
      printf ("Could not open %s\n", datafile);
      exit (0);
    }


  k = 0;
  while (fgets (aline, LINE, fptr) != NULL)
    {
      if (aline[0] != '#')
	{
	  sscanf (aline, "%d %d %lf %lf %lf %lf", &i, &j, &rho, &vr, &vtheta,
		  &vphi);
	  proga_ptr[i * MAXPROGA + j].rho = rho;
	  proga_ptr[i * MAXPROGA + j].v[0] = vr;
	  proga_ptr[i * MAXPROGA + j].v[1] = vtheta;
	  proga_ptr[i * MAXPROGA + j].v[2] = vphi;
	  k++;
	}
    }

  printf ("Read %d values from model %s\n", k, datafile);
  k = 82 * MAXPROGA + 19;
//printf("debug rho %e v %e %e %e\n",
  //      proga_ptr[k].rho,
  //      proga_ptr[k].v[0],
  //      proga_ptr[k].v[1],
  //      proga_ptr[k].v[2]);
  fclose (fptr);
  return (0);
}



/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	proga_velocity calculates the wind velocity at any position x
	in cartesian coordiantes
Arguments:		

Returns:
 
Description:	
Notes:
History:
 	99dec	ksl	Began work
	04dec	ksl	52a -- Mod to prevent divide by zero when 
			r=0.0, and NaN when xxx>1.0;
**************************************************************/


double
proga_velocity (x, v)
     double x[];
     double v[];
{
  double length ();
  int ii, jj;
  int im, jm;
  double f1, f2;
  double r, theta;
  double vr, vtheta, vphi;
  double speed;
  double xxx;

  if ((r = length (x)) == 0.0)
    {
      v[0] = v[1] = v[2] = 0.0;
      return (0.);
    };


  xxx = (sqrt (x[0] * x[0] + x[1] * x[1]) / r);
  if (xxx >= 1.)
    {
      theta = PI / 2.;
    }
  else
    theta = asin (xxx);

  printf ("x %.2g %.2g %.2g  -> r= %.2g theta = %.2g\n", x[0], x[1], x[2], r,
	  theta);
  im = jm = ii = jj = 0;
  while (proga_r[ii] < r && ii < iproga_r)
    ii++;
  while (proga_theta[jj] < theta && jj < iproga_theta)
    jj++;



  if (ii > 0 && ii < iproga_r)
    {				//r is in the normal range

      f1 = (r - proga_r[ii - 1]) / (proga_r[ii] - proga_r[ii - 1]);
      im = ii - 1;
    }
  else if (ii == iproga_r)
    {				// r is greater than anything in Proga's model

      f1 = 1;
      im = ii - 1;
    }
  else
    {
      f1 = 1;
      im = 0;
    }

  if (jj > 0 && jj < iproga_theta)
    {				// theta is inside the normal range

      f2 =
	(theta - proga_theta[jj - 1]) / (proga_theta[jj] -
					 proga_theta[jj - 1]);
      jm = jj - 1;
    }
  else if (jj == iproga_theta)
    {				//theta is more than the maximum theta

      f2 = 1;
      jm = jj - 1;
    }
  else
    {				//theta is less than the expected theta

      f2 = 1;
      jm = 0;
    }

  vr =
    (1. - f1) * ((1. - f2) * proga_ptr[im * MAXPROGA + jm].v[0] +
		 f2 * proga_ptr[im * MAXPROGA + jj].v[0]) + f1 * ((1. -
								   f2) *
								  proga_ptr[ii
									    *
									    MAXPROGA
									    +
									    jm].
								  v[0] +
								  f2 *
								  proga_ptr[ii
									    *
									    MAXPROGA
									    +
									    jj].
								  v[0]);

  vtheta =
    (1. - f1) * ((1. - f2) * proga_ptr[im * MAXPROGA + jm].v[1] +
		 f2 * proga_ptr[im * MAXPROGA + jj].v[1]) + f1 * ((1. -
								   f2) *
								  proga_ptr[ii
									    *
									    MAXPROGA
									    +
									    jm].
								  v[1] +
								  f2 *
								  proga_ptr[ii
									    *
									    MAXPROGA
									    +
									    jj].
								  v[1]);

  vphi =
    (1. - f1) * ((1. - f2) * proga_ptr[im * MAXPROGA + jm].v[2] +
		 f2 * proga_ptr[im * MAXPROGA + jj].v[2]) + f1 * ((1. -
								   f2) *
								  proga_ptr[ii
									    *
									    MAXPROGA
									    +
									    jm].
								  v[2] +
								  f2 *
								  proga_ptr[ii
									    *
									    MAXPROGA
									    +
									    jj].
								  v[2]);



  v[0] = vr * sin (theta) - vtheta * cos (theta);
  v[1] = vphi;
  v[2] = vr * cos (theta) + vtheta * sin (theta);

  speed = sqrt (v[0] * v[0] + v[1] * v[1] * v[2] * v[2]);

  if (sane_check (speed))
    {
      Error ("proga_velocity: v %e %e %e\n", v[0], v[1], v[2]);
    }
  return (speed);

}


double
proga_rho (x)
     double x[];
{
  double length ();
  int ii, jj;
  int im, jm;
  double r, theta, rrho;
  double f1, f2;
  r = length (x);
  theta = asin (sqrt (x[0] * x[0] + x[1] * x[1]) / r);
  printf ("x %.2g %.2g %.2g  -> r= %.2g theta = %.2g\n", x[0], x[1], x[2], r,
	  theta);

  if (r > proga_r[iproga_r])
    {
      printf (" r outside proga grid\n");
      rrho = 1.e-23;
      return (rrho);
    }

  im = jm = ii = jj = 0;
  while (proga_r[ii] < r && ii < iproga_r)
    ii++;
  while (proga_theta[jj] < theta && jj < iproga_theta)
    jj++;

  if (ii > 0)
    {
      f1 = (r - proga_r[ii - 1]) / (proga_r[ii] - proga_r[ii - 1]);
      im = ii - 1;
    }
  else
    f1 = 1;

  if (jj > 0)
    {
      f2 =
	(theta - proga_theta[jj - 1]) / (proga_theta[jj] -
					 proga_theta[jj - 1]);
      jm = jj - 1;
    }
  else
    f2 = 1;


//rrho=proga_ptr[ii*MAXPROGA+jj].rho;

  rrho =
    (1. - f1) * ((1. - f2) * proga_ptr[im * MAXPROGA + jm].rho +
		 f2 * proga_ptr[im * MAXPROGA + jj].rho) + f1 * ((1. -
								  f2) *
								 proga_ptr[ii
									   *
									   MAXPROGA
									   +
									   jm].
								 rho +
								 f2 *
								 proga_ptr[ii
									   *
									   MAXPROGA
									   +
									   jj].
								 rho);

  if (rrho < 1e-23)
    rrho = 1e-23;

  printf ("Grid point %d %d rho %e\n", ii, jj, rrho);

  return (rrho);
}
