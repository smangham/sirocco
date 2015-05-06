/* The routines in this file define and summarize the properties of the wind.  The routines here are specific to the simple keplerian
   description of a wind. It is only useful in the 2-d version of the code.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	get_sv_wind_params gets input data which is necessary for a Shlossman & Vitello 
	description of the wind
Arguments:		

Returns: 
Description:	
	The parameters, geo.sv...,  obtained here are only used in the routines in stellar_winds.c
	which calculate the velocity and density of the wind during the initialization process.
	Other portions of the structure, geo defined here are more general purpose.		
Notes:


History:
 	98dec	ksl	Coded and debugged as part of major change in IO structure required when
 				adding a spherical wind
        080518  ksl     60a - geo should contain only cgs units
	11aug	ksl	70b - kluge to get better xscale with compton torus
**************************************************************/
int
get_wind_keplerian_params (void)
{
	if(geo.coord_type != 1)
		Error ("wind_keplerian: Unsupported coordinate type %i- not cylindrical (1)\n", geo.coord_type);

 	geo.wind_thetamin = 0.;
 	geo.wind_thetamax = 0.;	//If both zero, set to maximum
 	geo.wind_rho_min = 2.8e9 / geo.rstar;
 	geo.wind_rho_max = 8.4e9 / geo.rstar;
 	rddoub ("wind_keplerian.diskmin(units_of_rstar)", &geo.wind_rho_min);
 	rddoub ("wind_keplerian.diskmax(units_of_rstar)", &geo.wind_rho_max);
 	geo.wind_rho_min *= geo.rstar;
	geo.wind_rho_max *= geo.rstar;


	geo.wind_keplerian_height = 1.0e3;
	geo.wind_keplerian_density = 1e0;
	geo.wind_keplerian_photon_bias = 0.1;
	rddoub ("wind_keplerian.density()", &geo.wind_keplerian_density);
	rddoub ("wind_keplerian.height()", &geo.wind_keplerian_height);
	rddoub ("wind_keplerian.photon_bias()", &geo.wind_keplerian_photon_bias);
	geo.wind_keplerian_photon_bias_const = geo.wind_keplerian_photon_bias /
		(exp (geo.wind_keplerian_photon_bias) - exp (-geo.wind_keplerian_photon_bias));

	geo.wind_rmin = geo.rstar;
	geo.wind_rmax = geo.rmax;
	geo.xlog_scale = geo.wind_rmin;
  	geo.zlog_scale = geo.wind_keplerian_height;
  	DFUDGE = geo.zlog_scale / 10.;

  	Log
    	("wind_keplerian: Creating wind for simple Keplerian disk of radius %g-%g\n",
     	geo.wind_rho_min, geo.wind_rho_max);
  return (0);
}

/***********************************************************
University of Southampton

Synopsis:
	double wind_keplerian_velocity(x,v) calulates the keplerian velocity
Arguments:		
	double x[]		the postion for which one desires the velocity
Returns:
	double v[]		the calculated velocity
	
	The amplitude of the velocity is returned 
	
History:
	11-2014	Cloned from sv_velocity by SWM
**************************************************************/
double
wind_keplerian_velocity (double x[], double v[])
{
  double r, speed;
  double xtest[3];

  r = sqrt (x[0] * x[0] + x[1] * x[1]);	//Convert position into radius

  v[0] = 0.;			//Zero R, Phi components
  v[2] = 0.;
  if (r > geo.wind_rho_min && r < geo.wind_rho_max && 
  		fabs (x[2]) < (geo.wind_keplerian_height + DFUDGE))	//If point is within the wind 
    v[1] = sqrt (G * geo.mstar / r);	//Simple keplerian Theta component
  else
    v[1] = 0;			//Or it's zero

  if (x[1] != 0.0)
    {
      project_from_cyl_xyz (x, v, xtest);
      stuff_v (xtest, v);
    }
  speed = (sqrt (v[0] * v[0] + v[1] * v[1]));
  return (speed);
}

/***********************************************************
University of Southampton

Synopsis:
	double wind_keplerian_rho(x) returns a flat density at any x
Arguments:		
	double x[]	the position where for the which one desires the denisty
Returns:
	The density at x is returned in gram/cm**3
	
History:
	11-2014	Cloned from proga_rho by SWM
**************************************************************/
double
wind_keplerian_rho (double x[])
{
  double r = sqrt (x[0] * x[0] + x[1] * x[1]);	//Convert position into radius
  if (r < geo.wind_rho_min || r > geo.wind_rho_max
      || fabs (x[2]) > geo.wind_keplerian_height)
    return (0.0);		//If the radius lies outside the wind, zero
  else
    return (geo.wind_keplerian_density);		//Else return flat value
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
 	cylind_volume(w) calculates the wind volume of a cylindrical cell
	allowing for the fact that some cells 

 Arguments:		
	WindPtr w;    the entire wind
 Returns:

 Description:
	This is a brute-froce integration of the volume 

	ksl--04aug--some of the machinations regarding what to to at the 
	edge of the wind seem bizarre, like a substiture for figuring out
	what one actually should be doing.  However, volume > 0 is used
	for making certain choices in the existing code, and so one does
	need to be careful.  
	
 History:
 	11-2014	Cloned from wind_cyl volume by SWM

 
**************************************************************/
#define RESOLUTION   1000	// Resolution for scanning partial disk section at

int
wind_keplerian_cyl_volumes (WindPtr w, int icomp)
{
  int i, j, n;
  double rmax, rmin;
  double zmin, zmax;

  for (i = 0; i < NDIM; i++)
    {
      for (j = 0; j < MDIM; j++)
	{
	  wind_ij_to_n (i, j, &n);	// Find wind index, calc wind cell volume
	  if (i == NDIM - 1)
	    {
	      w[n].inwind = W_NOT_INWIND;
	      w[n].vol = 0.0;
	    }
	  else
	    {
	      rmin = wind_x[i];
	      rmax = wind_x[i + 1];
	      zmin = wind_z[j];
	      zmax = wind_z[j + 1];
	      w[n].inwind = W_ALL_INWIND;
	      w[n].vol = 2 * PI * (rmax * rmax - rmin * rmin) * (zmax - zmin);
	    }
	}
    }
  return (0);
}


int				//http://www.nucleonica.net/wiki/images/8/89/MCNPvolI.pdf
wind_keplerian_randvec (PhotPtr pp, double r)
{
  double costheta, sintheta, x[3];
  double phi, sinphi, cosphi, d_rot;
  double p_biased;
  double k = geo.wind_keplerian_photon_bias, 
  		 c = geo.wind_keplerian_photon_bias_const;

  d_rot = 2. * PI * (rand () / MAXRAND);
  phi = 2. * PI * (rand () / MAXRAND);
  sinphi = sin (phi);
  cosphi = cos (phi);

  p_biased = rand () / MAXRAND;
  costheta = log ((k * p_biased / c) + exp (-k)) / k;
  sintheta = sqrt (1. - costheta * costheta);

  x[0] = r * costheta;
  x[1] = r * cosphi * sintheta;
  x[2] = r * sinphi * sintheta;

  pp->x[0] = cos (d_rot) * x[0] - sin (d_rot) * x[1];
  pp->x[1] = sin (d_rot) * x[0] + cos (d_rot) * x[1];
  pp->x[2] = x[2];
  pp->w *= 0.5 / (c * exp (k * costheta));

  //if(pp->np < 100) printf("POS: %g %g %g %g %g %g\n",pp->x[0]/r,pp->x[1]/r,pp->x[2]/r,pp->w,c,k);
  //Adjust weight by relative probality densities (.5 for random, c e^-kt for biased)
  return (0);
}
