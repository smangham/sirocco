
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
#include "balance_templates.h"

#define LINELENGTH 132
#define NP    10000

/* xsaha determines the abundances using the electron temperature and the Saha equation
   and then calculates the heating and cooling 

  02jun	ksl	Moved initializations of heading etc. after ionization calculation in
 		an attempt to get xsaha to behave more like python since some of ionization
	        calculations use the current value of the the total heating.	
 */

int
xsaha (www, p, nh, t_r, t_e, weight, ioniz_mode, freq_sampling)
     WindPtr www;
     PhotPtr p;
     double nh;
     double t_r;
     double t_e;
     double weight;
     int ioniz_mode;
     int freq_sampling;
{
  double total_fb ();
  double total_free ();
  double total_line_emission ();
  double total_emission ();
  double s;
  double alpha_recomb ();
  double planck (), kappa_ff ();
  double luminosity;
  double fb_cooling ();
//  double freqmin, freqmax;
  int n;

/* Now determine the abundances, concentrations, and partition functions based
on the conditions as currently defined.

ion_abundances is a python routine.  Currently (Sept 01) legal values are

0   on-the-spot approximation using existing t_e
1   LTE (e.g. concentrations)
2   Fixed concentrations
3   On the spot, with one_shot at updating t_e before calculating densities

It is far from clear that all of these will work at present with balance.  

01sept ksl

*/

  if (ioniz_mode < 0)
    {
      Log
	("xsaha: Calculating heating and cooling without changing abundances\n");
    }
  else
    {
      if (ioniz_mode == 0)
	{
	  Log
	    ("xsaha: Calculating abundances from On_the_spot or nebular approx\n");
	}
      else if (ioniz_mode == 1)
	{
	  Log ("xsaha: Calculating abundances from LTE\n");
	}
      else if (ioniz_mode == 2)
	{
	  Log ("xsaha: Fixed concentrations\n");
	}
      else if (ioniz_mode == 3)
	{
	  Log ("xsaha: Abundances using one shot (one_shot at t update\n");
	}
      else if (ioniz_mode == 4)
	{
	  Log ("xsaha: Test of detailed balance\n");
	}

      else
	{
	  Error ("xsaha: ioniz_mode unknown. Good luck!\n");
	}

      ion_abundances (www, ioniz_mode);
    }


/* initialize some parts of the wind ptr (see wind_rad_init) */

  www->j = www->ave_freq = www->lum = www->heat_tot = www->ntot = 0;
  www->heat_ff = www->heat_photo = www->heat_lines = 0.0;
  www->heat_z = 0;
  www->lum_z = 0.0;
  www->ntot = www->nioniz = 0;
  for (n = 0; n < nions; n++)
    {
      www->ioniz[n] = 0;
      www->recomb[n] = 0;
      www->heat_ion[n] = 0;
      www->lum_ion[n] = 0;
    }

  if (www->t_e != t_e)
    Error ("xsaha www->t_e %f != t_e %f. Why?\n", www->t_e, t_e);
  if (www->t_r != t_r)
    Error ("xsaha www->t_r %f != t_r %f. Why?\n", www->t_r, t_r);
  if (www->w != weight)
    Error ("xsaha www->w %f != weight %f. Why?\n", www->w, weight);

  s = www->vol;

/* Now calculate the heating  based on these abundances, etc.

There appear to be 3 types of heating, one in which the line is assumed
to have a simple flat topped profile, one in which the heating is calculated
using Python (i.e. there is no absorption until the photon scatters in a MC
sense), and one with a simple sobolev_line assumption (in which case one calculates
tau to define how much got absorbed.)

There appearted to be an error in the routine below which I fixed up, just so that
the choices actually matched what was said.  I fixed this.  01sept ksl

*/

/* Create photons which comprise a weighted BB distribution */

  init_bands (t_r, 0.0, 1e50, 0, &xband);
//  set_freq_lim (t_r, &freqmin, &freqmax);
//???? New  make_bb (p, t_r, weight);
  xbb (p, t_r, weight, xband.f1[0], xband.f2[0], 1);

  if (geo.rt_mode == 2)
    {
      for (n = 0; n < NPHOT; n++)
	{
	  sobolev_line_heating (www, &p[n], s);
	  radiation (www, &p[n], s);
	}
    }
  else if (geo.rt_mode == 1)
    {
      trans_phot (www, p, 0);
      www[0].heat_tot += www[0].heat_lines;
    }
  else
    {
      for (n = 0; n < NPHOT; n++)
	{
	  line_heating (www, &p[n], s);
	  radiation (www, &p[n], s);
	}
    }

  luminosity = total_emission (www, 0.0, INFINITY);
  num_recomb (&www[0], www->t_e);

  summary (www);

  return (0);
}


/* 
The routine cycle is intended to mimic what happens in Python in
a very simple way, updating the ion abundances, etc. in a single
cell, calculating the heating of the cell based on the updated
abundances and then calling one shot to try to adjust the temperature.

If t_r or t_e are changed then the routine will initialize everything
to LTE, but otherwise it will go straight to the heating, and updating
of abundances based on the new value of the heating.

NB: t_r, t_e, weight are only used for generating the photon field. 
The ionization calculation is based solely on what is contained in
www.

	01dec	ksl Updated to reflect new ways of calling routines
		    like concentrations.	
*/

double cy_tr = 0.0;
double cy_nh = 0.0;

int
cycle (www, p, nh, t_r, t_e, weight, mode, freq_sampling)
     WindPtr www;
     PhotPtr p;
     double nh;
     double t_r;
     double t_e;
     double weight;
     int mode;
     int freq_sampling;
{
  double total_fb ();
  double total_free ();
  double total_line_emission ();
  double total_emission ();
  double ne, s;
  double alpha_recomb ();
  double planck (), kappa_ff ();
  double luminosity;
//  double freqmin, freqmax;
  int n;


/* initialize some parts of the wind ptr (see wind_rad_init) */

  www->j = www->ave_freq = www->lum = www->heat_tot = www->ntot = 0;
  www->heat_ff = www->heat_photo = www->heat_lines = 0.0;
  www->heat_z = 0;
  www->lum_z = 0.0;
  www->ntot = www->nioniz = 0;
  for (n = 0; n < nions; n++)
    {
      www->ioniz[n] = 0;
      www->recomb[n] = 0;
      www->heat_ion[n] = 0;
      www->lum_ion[n] = 0;
    }

  if (www->t_e != t_e)
    Error ("xsaha www->t_e %f != t_e %f. Why?\n", www->t_e, t_e);
  if (www->t_r != t_r)
    Error ("xsaha www->t_r %f != t_r %f. Why?\n", www->t_r, t_r);
  if (www->w != weight)
    Error ("xsaha www->w %f != weight %f. Why?\n", www->w, weight);

  s = www->vol;


/* Now set up the initial concentrations.  This is to parallel what is
done in python where we initially begin with an LTE like assumption.
On successive calls to cycle this section should be jumped by the if
statement
 */
  if (t_r != cy_tr || nh != cy_nh)
    {
      Log
	("Cycle: Fixing abundances to Saha values since new t_r %g or nh %g\n",
	 t_r, nh);
      www->t_e = 0.9 * t_e;	// Lucy guess
      nebular_concentrations (www, 2);
      Log ("On the spot estimage  t_r of %g gives ne %g\n", t_r, ne);
      cy_tr = t_r;
      cy_nh = nh;
      www->dt_e = 0.0;		// Allow restart of oneshot limits
    }


/* Now calculate the heating */

/* First make some photons which comprise a BB distribution */

  init_bands (t_r, 0.0, 1e50, 0, &xband);
//  set_freq_lim (t_r, &freqmin, &freqmax);     //Needed even if you don't want to regenerate photons
//Error?? Uncommment line below to get photon distrubutions regenerated.
  xbb (p, t_r, weight, xband.f1[0], xband.f2[0], freq_sampling);


/* Shift values to old */
  www->heat_tot_old = www->heat_tot;
  www->dt_e_old = www->dt_e;

/* Next calculate the heating for this distribution */
  for (n = 0; n < NPHOT; n++)
    {
      line_heating (www, &p[n], s);
      radiation (www, &p[n], s);
    }

/* Now calculate the luminosity for these conditions */

//  luminosity = total_emission (www, freqmin, freqmax);
  luminosity = total_emission (www, xband.f1[0], xband.f2[0]);
  num_recomb (&www[0], www->t_e);

  summary (www);

/* Shift values to old */
  www->dt_e = www->t_e - www->t_e_old;	//Must store this before others
  www->t_e_old = www->t_e;
  www->t_r_old = www->t_r;
  www->lum_rad_old = www->lum_rad;

  one_shot (www, mode);


  Log ("Cycle: one_shot old t_r t_e %8.2g %8.2g, new t_r t_e %8.2g %8.2g\n",
       t_r, t_e, www->t_r, www->t_e);

//Error ?? -- Next statement is probably superfluous, since calc_te, called by oneshot calculated total emission
//  luminosity = total_emission (www, freqmin, freqmax);
  luminosity = total_emission (www, xband.f1[0], xband.f2[0]);
  num_recomb (&www[0], www->t_e);

  summary (www);

/* Convergence check */
  convergence (www);

  return (0);

}



/* 

The routine dumb_step increments te upward or downward depending on whether the total luminosity at
this temperature is less than or greater than the current value of the heating. It does not
recalculate the heating at any point.   The resulting te is incorporated into the windptr
	01dec	ksl	Updated to reflect new calls to routines like total_emission.  The routine
			seemed to work following these changes.
*/

int
dumb_step (www, te)
     WindPtr www;
     double *te;
{
  double lum_current, lum_hotter, lum_cooler;
  double tstep, hardness;
  double total_emission ();

  if (*te > 20000.)
    tstep = 0.1 * (*te);
  else
    tstep = 0.05 * (*te);

  www->t_e = *te;
  lum_current = total_emission (www, 0.0, INFINITY);
  num_recomb (&www[0], www->t_e);

  Log
    ("Current values of luminosities compared to previously calculated heating\n");
  Log
    ("t_e %8.2e lum_tot  %8.2e lum_lines  %8.2e lum_ff  %8.2e lum_fb     %8.2e %8.2e %8.2e %8.2e %8.2e\n",
     www->t_e, www->lum_rad, www->lum_lines, www->lum_ff, www->lum_fb,
     www->lum_ion[0], www->lum_ion[2], www->lum_ion[3], www->lum_z);
  Log
    ("t_r %8.2e heat_tot %8.2e heat_lines %8.2e heat_ff %8.2e heat_photo %8.2e %8.2e %8.2e %8.2e %8.2e\n",
     www->t_r, www->heat_tot, www->heat_lines, www->heat_ff,
     www->heat_photo, www->heat_ion[0], www->heat_ion[2], www->heat_ion[3],
     www->heat_z);

  www->t_e = *te + tstep;
  lum_hotter = total_emission (www, 0.0, INFINITY);
  num_recomb (&www[0], www->t_e);

  Log ("Results of raising t_e\n");
  Log
    ("t_e %8.2e lum_tot  %8.2e lum_lines  %8.2e lum_ff  %8.2e lum_fb     %8.2e %8.2e %8.2e %8.2e %8.2e\n",
     www->t_e, www->lum_rad, www->lum_lines, www->lum_ff, www->lum_fb,
     www->lum_ion[0], www->lum_ion[2], www->lum_ion[3], www->lum_z);
  Log
    ("t_r %8.2e heat_tot %8.2e heat_lines %8.2e heat_ff %8.2e heat_photo %8.2e %8.2e %8.2e %8.2e %8.2e\n",
     www->t_r, www->heat_tot, www->heat_lines, www->heat_ff,
     www->heat_photo, www->heat_ion[0], www->heat_ion[2], www->heat_ion[3],
     www->heat_z);

  www->t_e = *te - tstep;
  lum_cooler = total_emission (www, 0.0, INFINITY);
  num_recomb (&www[0], www->t_e);

  Log ("Results of lowering t_e\n");
  Log
    ("t_e %8.2e lum_tot  %8.2e lum_lines  %8.2e lum_ff  %8.2e lum_fb     %8.2e %8.2e %8.2e %8.2e %8.2e\n",
     www->t_e, www->lum_rad, www->lum_lines, www->lum_ff, www->lum_fb,
     www->lum_ion[0], www->lum_ion[2], www->lum_ion[3], www->lum_z);
  Log
    ("t_r %8.2e heat_tot %8.2e heat_lines %8.2e heat_ff %8.2e heat_photo %8.2e %8.2e %8.2e %8.2e %8.2e\n",
     www->t_r, www->heat_tot, www->heat_lines, www->heat_ff,
     www->heat_photo, www->heat_ion[0], www->heat_ion[2], www->heat_ion[3],
     www->heat_z);

  Log ("te %8.2g lum  %8.2e\n", *te, lum_current);
  Log ("te %8.2g lum  %8.2e\n", *te + tstep, lum_hotter);
  Log ("te %8.2g lum  %8.2e\n", *te - tstep, lum_cooler);

  Log ("tr %8.2g heat %8.2e\n", www->t_r, www->heat_tot);

  hardness = (lum_current - www->heat_tot) / (lum_current + www->heat_tot);

  if (fabs (hardness) < 0.1)
    {
      Log ("Convergence is achieved!\n");
    }
  else if (hardness > 0.0)
    {				// luminosity is greater the heating 

      if (lum_cooler < lum_hotter)
	{
	  *te -= tstep;
	  Log ("Dropping t_e to %8.2g\n", *te);
	}
      else
	{
	  *te += tstep;
	  Log ("Raising t_e to % 8.2g\n", *te);
	}
    }
  else
    {				//luminosity is less than heating

      if (lum_cooler > lum_hotter)
	{
	  *te -= tstep;
	  Log ("Dropping t_e to %8.2g\n", *te);
	}
      else
	{
	  *te += tstep;
	  Log ("Raising t_e to % 8.2g\n", *te);
	}
    }

  www->t_e = *te;
  return (0);
}

int
find_te (www)
     WindPtr www;
{
  double calc_te ();
  double luminosity, total_emission ();

  luminosity = total_emission (www, 0.0, INFINITY);
  num_recomb (&www[0], www->t_e);

  summary (www);

  www->t_e = calc_te (www, TMIN, 1.2 * www->t_r);

  luminosity = total_emission (www, 0.0, INFINITY);
  num_recomb (&www[0], www->t_e);

  summary (www);

  return (0);
}
