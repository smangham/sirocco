
/*

                                       Space Telescope Science Institute

Synopsis:
	These are all of the routines necessary to define a KWD wind
Arguments:		

Returns:
 
Description:	
Notes:


History:
 	01mar	ksl	Added this possibility
	02jan	ksl	Rechecked all of the routines in an atttempt to
			determine why I do not produce KWD figure 7 well.
        11may   nsh     Added the option to have a multiple of sound speed as initial velocity
	15aug	ksl	Mods to accomodate multiple domains

 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/* This is a variable required to do the integration of local mdot in KWD */
double kn_lambda;               //kn_lambda required for the intgration

/***********************************************************
	Space Telescope Science Institute

 Synopsis:
	get_knigge_wind_params gets input data that is necessary for a Knigge's 
	description of the wind
Arguments:		

Returns:
 
Description:	
	The parameters, kn...,  obtained here are only used in the routines in stellar_winds.c
	which calculate the velocity and density of the wind during the initialization process.
	Other portions of the structure, geo defined here are more general purpose.		
Notes:


History:
 	01mar	ksl	Added this possibility
	04jun	ksl	Moved readding wind_mdot to this routine
        04Sep   SS      Changed finite disk case to cut off wind
                        at the top (rather than bottom) outer edge of disk.
	080518	ksl	60a - Fixed to account for fact that geo.wind_mdot 
			has been initialized in cgs units
**************************************************************/

int
get_knigge_wind_params (ndom)
     int ndom;
{
  double dmin;
  double disktheta, test;


  Log ("Creating KWD wind in domain %d\n", ndom);

  zdom[ndom].wind_mdot = 0.1 * geo.disk_mdot / (MSOL / YR);
  rddoub ("wind.mdot(msol/yr)", &zdom[ndom].wind_mdot);
  zdom[ndom].wind_mdot *= MSOL / YR;


  zdom[ndom].kn_r_scale = 7e10; /*Accleration length scale for wind */
  zdom[ndom].kn_alpha = 1.5;    /* Accleration scale exponent for wind */
  zdom[ndom].kn_v_infinity = 3; /* Final speed of wind in units of escape velocity */
  zdom[ndom].kn_lambda = 0.0;   /* Mass loss rate exponent */
  zdom[ndom].kn_dratio = dmin = 0.5 * sqrt (geo.diskrad / geo.rstar);   /* Center of collimation in units of the WD 
                                                                           radius. The value set here is for the minimum collimation, see KWD95.  The coefficient 0.5 is 
                                                                           approximate */
  zdom[ndom].kn_v_zero = 1.0;   /* NSH 19/04/11 Parameter which 
                                   can use to muliply the initial velocity of wind so that it
                                   is greater or less than the sound speed. */
  zdom[ndom].wind_rho_min = 1;  /* Innner and outer edges of the wind in stellar radii. These
                                   parameters were added by JM to allow one to duplicate the YSO
                                   paper */
  zdom[ndom].wind_rho_max = geo.diskrad / geo.rstar;


/* There is confusion in various papers concerning whether to use d or d/dmin.  In KWD95, d/dmin was
used but in later papers, e.g KD97 d in WD radii was used.  I believe d is more natural and so will use it, 
but one should remember that this differs from KWD.  To emphasize this we will calculate and log d/dmin.
Terminolgy is awful here. -- ksl 

As now represented kn_dratio is the distance to the focus point in stellar radii!

*/

  rddoub ("kn.d(in_units_of_rstar)", &zdom[ndom].kn_dratio);

  Log_silent ("dmin = %f so the ratio d/dmin here is %f  (%.2e %.2e) \n", dmin, zdom[ndom].kn_dratio / dmin, geo.diskrad, geo.rstar);


  rddoub ("kn.mdot_r_exponent", &zdom[ndom].kn_lambda); /* Mass loss rate exponent */
  rddoub ("kn.v_infinity(in_units_of_vescape)", &zdom[ndom].kn_v_infinity);     /* Final speed of wind in units of escape velocity */
  if (zdom[ndom].kn_v_infinity < 0)
  {
    Log ("Since kn_v_infinity is less than zero, will use SV prescription for velocity law.\n Velocity at base remains the sound speed\n");
  }

  rddoub ("kn.acceleration_length(cm)", &zdom[ndom].kn_r_scale);        /*Accleration length scale for wind */
  rddoub ("kn.acceleration_exponent", &zdom[ndom].kn_alpha);    /* Accleration scale exponent */
  rddoub ("kn.v_zero(multiple_of_sound_speed_at_base)", &zdom[ndom].kn_v_zero);

/* Assign the generic parameters for the wind the generic parameters of the wind. Note that one should not really
use these generic parameters in the rest of the routines here, as especially for the yso model they may have
to be modified -- ksl 04jun */


  /* JM 1502 -- added capability for user to adjust launch radii of KWD wind
     required to reproduce Stuart's X-ray models (Sim+ 2008,2010) 
     units are in stellar radii / WD radii / grav radii as in SV model 

     1509 -- ksl -- This had not be implemented quite right.  Now fixed.  Defaults
     are for standard KWD model */
  rddoub ("kn.rmin(in_units_of_rstar)", &zdom[ndom].wind_rho_min);
  rddoub ("kn.rmax(in_units_of_rstar)", &zdom[ndom].wind_rho_max);

  zdom[ndom].wind_rho_min *= geo.rstar;
  zdom[ndom].wind_rho_max *= geo.rstar;
  zdom[ndom].wind_thetamin = atan (1. / zdom[ndom].kn_dratio);

/* Somewhat paradoxically diskrad is in cm, while dn_ratio which is really d in KWD95 is 
in units of WD radii */
  zdom[ndom].wind_thetamax = atan (geo.diskrad / (zdom[ndom].kn_dratio * geo.rstar));

  /* Next lines added by SS Sep 04. Changed the wind shape so that the boundary touches the outer 
     corner of the disk rather than the intersection of the disk edge with the xy-plane. */

  if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
  {
    zdom[ndom].wind_thetamax = atan (geo.diskrad / (((zdom[ndom].kn_dratio * geo.rstar) + zdisk (geo.diskrad))));
  }


  zdom[ndom].rmin = zdom[ndom].wind_rho_min;
  /* JM 1710 -- removed for issue #305 
     //zdom[ndom].rmax = 10. * zdom[ndom].kn_r_scale;        // Set rmax to something reasoable. 
   */

  /* The change in the boundary of the wind (as corner of disk -- see above) 
     means that wind_rho_max nees to be redefined so that it is used correctly
     to compute the boundary of the wind elsewhere. */

  // XXX Next lines supercede definitions above and look wrong 
  if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
  {
    zdom[ndom].wind_rho_max = geo.diskrad - (zdisk (geo.diskrad) * tan (zdom[ndom].wind_thetamax));
  }

/* 70d - ksl - I changed the scaling to something that produced more cells
 * in the wind, at the cost of slightly less spatial resolution at the inner
 * edge of the wind
 */

  /* if modes.adjust_grid is 1 then we have already adjusted the grid manually */
  if (modes.adjust_grid == 0)
  {
    zdom[ndom].xlog_scale = geo.rstar;
    zdom[ndom].zlog_scale = geo.rstar;
  }


/*Now calculate the normalization factor for the wind*/

  test = geo.rstar;

  /* For non-flat disk some streamlines are missing (SS). */

  if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
  {
    disktheta = atan (zdisk (geo.diskrad) / geo.diskrad);
    test =
      zdom[ndom].kn_dratio * geo.rstar * sin (zdom[ndom].wind_thetamin) *
      cos (disktheta) / sin ((PI / 2.) - zdom[ndom].wind_thetamin - disktheta);
  }

  kn_lambda = zdom[ndom].kn_lambda;
  zdom[ndom].mdot_norm = qromb (kn_wind_mdot_integral, test, geo.diskrad, 1e-6);

  return (0);
}


/* The routines in this file define and summarize the properties of the wind.  The routines here are
   specific to the Knigge description of a wind. 
 */

/***********************************************************
	Space Telescope Science Institute

 Synopsis:
	double kn_velocity(x,v) calculates the v in cartesion coordinates
	of a Knigge wind from a position x in cartesian coordinates.  
Arguments:		
	ndom			The domain for the KN wind
	double x[]		the position where one desires the velocity
Returns:
	double v[]		the calculated velocity
	
	The amplitude of the velocity is returned 
	
Description:	

	This implements equations 17 & 18 of KWD95.  To get the constant
	velocity law set kn_alpha to 0.  

	58 -- If v_inf<0, then an SV law is used for the velocity along
	the streamline.	
		
Notes:

At the inner edge of the accretion disk the temperature of the disk
goes to zero, and so in principle the thermal velocity goes to zero.
The kfudge below keeps this from happening. What did Christian actually
do here?? 

04mar	ksl	The way in which the wind velocity is made to fill
		space, particular the space inside the wind cone,
		is not obvious.  The way in which this is done at
		present causes the div v calulation to be negative
		inside the windcone.  I have suppressed this after
		the first error in wind_div_v, but this could be
		a mistake.
History:
	01mar	ksl	Began coding
	02jan	ksl	Removed the implicit assumption that we always
			wanted the velocity "above" the disk plane.
	04aug	ksl	Modified to account for a disk with vertical
			extent, and to return the .
	05jul	ksl	Corrected error in the direction of ptest
	06oct	ksl	58 -- provided a very simple way to use a 
			SV velocity law (that is a law in which
			the scale height is always the place where
			the velocity is half the peak velocity) in
			what otherwise looks like s KWD model.  The
			SV law is used whenever v_infinity is less 
			than 0. Note that vzero is not changed from
			the KWD prescription.  
	08aug	ksl	Began conversion to allow domains
		
**************************************************************/

#define kfudge 1.01


double
kn_velocity (ndom, x, v)
     int ndom;
     double x[], v[];
{
  double r, rzero, theta;
  double ldist, zzz, v_escape, vl;
  double dd;
  double vzero;
  struct photon ptest;
  double xtest[3];
  double s;

  DomainPtr one_dom;


  one_dom = &zdom[ndom];

/*dd is the distance of the focus point along the z axis in cm */

  dd = geo.rstar * one_dom->kn_dratio;
  r = sqrt (x[0] * x[0] + x[1] * x[1]); // rho of the point we have been given

/* rzero is the position where the stream line  hits the xy plane. This could less than R_wd. */

  rzero = r / (1. + fabs (x[2] / dd));
  theta = atan (rzero / dd);

/* ldist is the poloidal distance along the streamline */

  ldist = sqrt ((r - rzero) * (r - rzero) + x[2] * x[2]);

/* 04aug -- ksl -- 52 Take the thickness of the disk into account if that is necessary.  */

  if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
  {
    xtest[0] = r;               // Define xtest in the +xz plane
    xtest[1] = 0;
    xtest[2] = fabs (x[2]);
    ptest.x[0] = rzero;         // Define ptest at the intersection of the streamline and x axis
    ptest.x[1] = 0.0;
    ptest.x[2] = EPSILON;
    ptest.lmn[0] = sin (theta); //56d -- ksl -- lmn is along the streamline toward xtest
    ptest.lmn[1] = 0.0;
    ptest.lmn[2] = cos (theta);
    s = ds_to_disk (&ptest, 1);
    move_phot (&ptest, s);      // Now test photon is at disk surface
    vsub (ptest.x, xtest, xtest);
    ldist = length (xtest);
    rzero = length (ptest.x);
  }


/* The if statement, below  limits the poloidal component of the velocity to be no more than 
the poloidal velocity at the inner edge of the wind. It is there for continuity reasons */

  if (rzero < geo.rstar)
    v_escape = sqrt (2. * G * geo.mstar / geo.rstar);
  else
    v_escape = sqrt (2. * G * geo.mstar / rzero);

/* Note that vzero for kwd the sound speed dos not depend on any of the kwd parameters */

  vzero = kn_vzero (rzero);

/*NSH 19/4/11 In the standard KWD model, the speed at the base of the wind is the
   * sound speed.  The next line modifies this by the factor given by kn_v_zero, which
   * is an input to the code. */

  vzero *= one_dom->kn_v_zero;

/* 578 -- 06oct -- ksl -- The next lines are modified to allow one to create a SV style
velocity law if kn_v_infinity is less than 0 */

  if (one_dom->kn_v_infinity > 0.0)
  {
    zzz = ldist / (ldist + one_dom->kn_r_scale);
    zzz = pow (zzz, one_dom->kn_alpha); // In Knigge's notation this is actually beta
    vl = vzero + (one_dom->kn_v_infinity * v_escape - vzero) * zzz;
  }
  else
  {
    zzz = pow (ldist / one_dom->kn_r_scale, one_dom->kn_alpha);
    vl = vzero + ((-one_dom->kn_v_infinity) * v_escape - vzero) * zzz / (1. + zzz);
  }


  v[0] = vl * sin (theta);

/* Unlike the polidal component we had been treating the WD as a point source for 
the azimuthal velocity.  One could do this more consistently, with the following
code.: 
  if (r > 0)
    v[1] = sqrt (G * geo.mstar * rzero) / r;  // Eqn 8 KWD
  else
    v[1] = 0;

As far as I could tell it made no difference  02apr-ksl  */

  if (rzero > geo.rstar)
    v[1] = sqrt (G * geo.mstar * rzero) / r;    // Eqn 8 KWD
  else if (r > 0)
    v[1] = sqrt (G * geo.mstar * geo.rstar) / r;
  else
    v[1] = 0;
  v[2] = vl * cos (theta);
  if (x[2] < 0)
    v[2] *= (-1);

/* 

04aug -- ksl --52 At this point we have calculated the velocity in the xz plane, which
is identical to the statement that we have calculated it in cylindrical coordinates.  
We now project back to xyz coordinates if necessary.  

06oct -- ksl -- In practice, the current version of python currently should never need to carry out  
this projection.  However, it is needed, if we were to move to full 3d calculation and for certain 
test programs.

*/
  if (x[1] != 0.0)
  {
    project_from_cyl_xyz (x, v, xtest);
    stuff_v (xtest, v);
  }

  return (sqrt (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));
}

/***********************************************************
	Space Telescope Science Institute

 Synopsis:
	double kn_rho(x) calculates the density of an kn_wind at a position x
Arguments:		
	ndom		The domain where kn params are stored
	double x[]	the position which one desires the density
Returns:
	The density at x is returned in gram/cm**3
	
Description:	
		
Notes:

History:
 	01mar	ksl	Began work.
	15aug	ksl	Began effort to allow domains
 
**************************************************************/

double
kn_rho (ndom, x)
     int ndom;
     double x[];
{
  double r, rzero;
  double dd;
  double vqr, vzero;
  double v[3], rho;
  struct photon ptest;
  double s, theta;

  DomainPtr one_dom;


  one_dom = &zdom[ndom];


  dd = geo.rstar * one_dom->kn_dratio;
  r = sqrt (x[0] * x[0] + x[1] * x[1]); //rho coordinate of the point we have been given
  rzero = r / (1. + fabs (x[2] / dd));  //rho at the base for this streamline
  /* If the disk is thick we need to modify the position of rzero */

  if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
  {
    theta = atan (rzero / dd);
    ptest.x[0] = rzero;
    ptest.x[1] = 0.0;
    ptest.x[2] = EPSILON;
    ptest.lmn[0] = cos (theta);
    ptest.lmn[1] = 0.0;
    ptest.lmn[2] = sin (theta);
    s = ds_to_disk (&ptest, 1);
    move_phot (&ptest, s);      // Now test photon is at disk surface
    rzero = sqrt (ptest.x[0] * ptest.x[0] + ptest.x[1] * ptest.x[1]);
  }
  kn_velocity (ndom, x, v);
  vqr = sqrt (v[0] * v[0] + v[2] * v[2]);       //poloidal velocity
  vzero = kn_vzero (rzero);     //polidal velocity at base
  if ((r * vqr) > 0)            // If statement to assure denominator is not zero
  {
    rho = kn_rho_zero (ndom, rzero) * (rzero * rzero * vzero) / (r * r * vqr);
  }
  else
    rho = 0.0;
  return (rho);
}



/***********************************************************
	Space Telescope Science Institute

 Synopsis:
	Calculate the sound speed (launch speed) of the wind at r
Arguments:		
 
Returns:
 
Description:	
	For KWD, the launch speed is the sound velocity
	
		
Notes:

	kn_vzero does not depend on the domain

History:
 	01mar      ksl	Coding began.
 
**************************************************************/


double
kn_vzero (r)
     double r;
{
  double tref, tdisk ();
  double t, teff ();
  double ratio, v;
  tref = tdisk (geo.mstar, geo.disk_mdot, geo.rstar);
  ratio = r / geo.rstar;
  if (ratio < kfudge)
    ratio = kfudge;
  t = teff (tref, ratio);
  v = 1.e6 * sqrt (t / 1e4);    // Frank, King & Raine 1985
  return (v);
}


/***********************************************************
	Space Telescope Science Institute

 Synopsis:
	double kn_wind_mdot_integral(r) is the integrand of SV model 
	for mdot as a function of radius
Arguments:		
	double r;
 
Returns:
 
Description:	
	
		
Notes:
       There is an extra factor of 2 because the wind emerges 
       from both sides of the disk

History:
 	01mar      ksl	Coding began.
	15aug	ksl	Added extrenal variable kn_wind_mdot_ndon
			in orred to integrate the surface 
 
**************************************************************/


double
kn_wind_mdot_integral (r)
     double r;
{
  double t;
  double x, ratio;
  double tref;


  tref = tdisk (geo.mstar, geo.disk_mdot, geo.rstar);
  ratio = r / geo.rstar;
  if (ratio < kfudge)
    ratio = kfudge;
  t = teff (tref, ratio);
  //x = 4. * PI * pow (t, 4. * zdom[kn_wind_mdot_ndom].kn_lambda) * r;
  x = 4. * PI * pow (t, 4. * kn_lambda) * r;
  return (x);
}


/***********************************************************
	Space Telescope Science Institute

 Synopsis:
	Calculate rho at the base of the wind
Arguments:		
 
Returns:
 
Description:	
	
		
Notes:

History:
 	01mar      ksl	Coding began.
 
**************************************************************/

double
kn_rho_zero (ndom, r)
     double r;
     int ndom;
{
  double tref;
  double t;
  double x, ratio;
  double vzero, dd, cosd;

  DomainPtr one_dom;


  one_dom = &zdom[ndom];


  tref = tdisk (geo.mstar, geo.disk_mdot, geo.rstar);
  ratio = r / geo.rstar;

  if (ratio < kfudge)
    ratio = kfudge;

  t = teff (tref, ratio);
  x = one_dom->wind_mdot * pow (t, 4. * one_dom->kn_lambda) / one_dom->mdot_norm;
  vzero = kn_vzero (r);
  dd = geo.rstar * one_dom->kn_dratio;
  cosd = dd / sqrt (r * r + dd * dd);
  x /= (vzero * cosd);
  return (x);
}
