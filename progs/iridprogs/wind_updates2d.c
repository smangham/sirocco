/* The routines in this file are generic.  There is no dependence on a particlar wind model or
 * any coordinate system dependences.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "atomic.h"
#include "python.h"
#define LINELEN 132

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis: wind_update(w) updates the parameters in the wind that are 
	affected by radiation.

 Arguments:		
	WindPtr w;

Returns:
 
Description:	
	
Notes:
	There is not much point in calling this until you have propagated a few photons
	The routine should not be called until at least some photons have been propagated 
	through the wind because the radiation and electron temperatures are calculated from 
	the frequency moment
History:
 	97jan      ksl	Coding on python began.
 	97sep13	ksl	Added the possibility of fixing the concentrations of the ions to specific
 			values, but this really should not be called because in that case there
 			is little point in updating the wind since the only thing which wind_update
 			really does is to change the ionization fractions.  This change was really
 			for consistency.
 	97oct2	ksl	Modified relation between t_e and t_r to use CK's prescription derived from
 			Cloudy comparisons.
 	98may3	ksl	Fixed problems with the portion of the routine which placed new densities in
 			the regions just outside the wind.  These densities are used for interpolation
 			purposes.
 	98may16	ksl	Added some checks to indicate how great the changes in t_e and t_r in each
 			wind_update
 	98may17	ksl	Changed some of the logic to (hopefully) eliminate some calculations outside
 			the wind which should not be necessary.   Added some additional tracking of 
 			changes to see how things converge, and where the biggest changes in t_e 
 			and t_r are.
 	98may24	ksl	Changed code so that w[].vol > 0 is the primary criterion for attempting to
 			recalculate the ionization abundances
 	98may31 ksl	Corrected error in calculating dilution factor because cells are actually two 
 			cells one above and one below the plane of the disk.
 	98dec	ksl	Relocated section of code which chose between various ways of calculating
 			ion abundances 	to a subroutine ion_abundances.  Removed section of code
 			in which contained estimates of what t_e ought to be based on t_r, w, etc.
	99jan	ksl	Modified how t_e and t_r change in the event there were no photons in
			a cell.  Old solution was to set t_e and t_r to zero; new solution was
			to leave t_r at old value and drop t_e by the maximum amount 0.7*t_e_old
			being used in one_shot.  
        04May   SS      Removed extra factor of two in calculation of the radiation density since 
			the cell
                        volume now includes ares both above and below the plane.
	05apr	ksl	55d: Began work to make more compatible with 1-d coordinate systems.
	05apr	ksl	56: Moved functionality that extends the ion densityies
			to cylindrical.c, rtheta.c, and spherical.c in 
			continuing effort to isolate difference between
			coordinate system types.
	06may	ksl	57+: Began splitting of wind and plasma structure.  This routine loops
			over the entire grid and ultimately eveything here is likely to be
			contained in the plasma structure.  Therefore here we want to use
			plasmamain


**************************************************************/


int num_updates = 0;

int
wind_update (w)
WindPtr (w);
{
  int n, i, j;
  double trad, nh;
  double wtest, xsum, asum, psum, fsum, lsum;
  double volume;
  char string[LINELEN];
  double t_r_old, t_e_old, dt_r, dt_e;
  double t_r_ave_old, t_r_ave, t_e_ave_old, t_e_ave;
  int iave, nmax_r, nmax_e;
  int nplasma;
  int nwind;

  dt_r = dt_e = 0.0;
  iave = 0;
  nmax_r = nmax_e = -1;
  t_r_ave_old = t_r_ave = t_e_ave_old = t_e_ave = 0.0;


  for (n = 0; n < NPLASMA; n++)
    {
      nwind = plasmamain[n].nwind;
      volume = w[nwind].vol;
      /* Start with a call to the routine which normalises all the macro atom 
         monte carlo radiation field estimators. It's best to do this first since
         some of the estimators include temperature terms (stimulated correction
         terms) which were included during the monte carlo simulation so we want 
         to be sure that the SAME temperatures are used here. (SS - Mar 2004). */
      if (geo.rt_mode == 2 && geo.macro_simple == 0)	//test for macro atoms
	{
	  mc_estimator_normalise (nwind);
	}

      /* Store some information so one can determine how much the temps are changing */
      t_r_old = plasmamain[n].t_r;
      t_e_old = plasmamain[n].t_e;
      t_r_ave_old += plasmamain[n].t_r;
      t_e_ave_old += plasmamain[n].t_e;
      iave++;

      if (plasmamain[n].ntot < 100)
	{
	  Log
	    ("!!wind_update: Cell %d  with volume %8.2e has only %d photons\n",
	     n, volume, plasmamain[n].ntot);
	}

      if (plasmamain[n].ntot > 0)
	{
	  wtest = plasmamain[n].ave_freq;
	  plasmamain[n].ave_freq /= plasmamain[n].j;	/* Normalization to frequency moment */

	  if (sane_check (plasmamain[n].ave_freq))
	    {
	      Error ("wind_update: %d ave_freq %e j %e ntot %d\n",
		     n, wtest, plasmamain[n].j, plasmamain[n].ntot);
	    }

	  plasmamain[n].j /= (4. * PI * volume);	//Factor of 2 has been removed from this line (SS, May04)

	  trad = plasmamain[n].t_r =
	    H * plasmamain[n].ave_freq / (BOLTZMANN * 3.832);
	  plasmamain[n].w =
	    PI * plasmamain[n].j / (STEFAN_BOLTZMANN * trad * trad * trad *
				    trad);

	  if (sane_check (trad) || sane_check (plasmamain[n].w))
	    {
	      Error ("wind_update: %d trad %8.2e w %8.2g\n", n, trad,
		     plasmamain[n].w);
	      Error ("wind_update: ave_freq %8.2e j %8.2e\n",
		     plasmamain[n].ave_freq, plasmamain[n].j);
	      exit (0);
	    }
	}
      else
	{			// It is not clear what to do with no photons in a cell

	  plasmamain[n].j = 0;
	  trad = plasmamain[n].t_r;
	  plasmamain[n].t_e *= 0.7;
	  if (plasmamain[n].t_e < TMIN)
	    plasmamain[n].t_e = TMIN;
	  plasmamain[n].w = 0;
	}

      nh = plasmamain[n].rho * rho2nh;

      /* If geo.adiabatic is true, then alculate the adiabatic cooling using the current, i.e 
       * previous value of t_e.  Note that this may not be  best way to determien the cooling. 
       * Changes made here should also be reflected in wind2d.c.  At present, adiabatic colling
       * is not included in updates to the temperature, even if the adiabatic cooling is calculated
       * here. 04nov -- ksl */
      /* 05apr -- ksl -- The index being used was incorrect.  This has been fixed now */

      if (geo.adiabatic)
	plasmamain[n].lum_adiabatic =
	  adiabatic_cooling (&w[n], plasmamain[n].t_e);
      else
	plasmamain[n].lum_adiabatic = 0.0;


      /* Calculate the densities in various ways depending on the ioniz_mode */


      ion_abundances (&plasmamain[n], geo.ioniz_mode);


      //Perform checks to see how much temperatures have changed in this iteration
      if ((fabs (t_r_old - plasmamain[n].t_r)) > fabs (dt_r))
	{
	  dt_r = plasmamain[n].t_r - t_r_old;
	  nmax_r = n;
	}
      if ((fabs (t_e_old - plasmamain[n].t_e)) > fabs (dt_e))
	{
	  dt_e = plasmamain[n].t_e - t_e_old;
	  nmax_e = n;
	}
      t_r_ave += plasmamain[n].t_r;
      t_e_ave += plasmamain[n].t_e;
    }


  /* Now we need to updated the densities immediately outside the wind so that the density interpolation in resonate will work.
     In this case all we have done is to copy the densities from the cell which is just in the wind (as one goes outward) to the
     cell that is just inside (or outside) the wind. 

     SS asked whether we should also be extending the wind for other parameters, especially ne.  At present we do not interpolate
     on ne so this is not necessary.  If we did do that it would be required.

     In cylindrical coordinates, the fast dimension is z; grid positions increase up in z, and then out in r.
     In spperical polar coordinates, the fast dimension is theta; the grid increases in theta (measured)
     from the z axis), and then in r.
     In spherical coordinates, the grid increases as one might expect in r..
     *
   */

  if (geo.coord_type == CYLIND)
    cylind_extend_density (w);
  else if (geo.coord_type == RTHETA)
    rtheta_extend_density (w);
  else if (geo.coord_type == SPHERICAL)
    spherical_extend_density (w);
  else if (geo.coord_type == CYLVAR)
    cylvar_extend_density (w);
  else
    {
      Error ("Wind_update2d: Unknown coordinate type %d\n", geo.coord_type);
      exit (0);
    }

  /* Finished updatating region outside of wind */

  num_updates++;
  strcpy (string, "");
  sprintf (string, "# Wind update: Number %d", num_updates);


  /* Check the balance between the absorbed and the emitted flux */

  xsum = psum = lsum = fsum = 0;
  for (i = 0; i < NDIM2; i++)
    nplasma = w[n].nplasma;
  {
    if (sane_check (plasmamain[nplasma].heat_tot))
      Error ("wind_update: w\[%d}.heat_tot is %e\n", i,
	     plasmamain[nplasma].heat_tot);
    if (sane_check (plasmamain[nplasma].heat_photo))
      Error ("wind_update: w\[%d}.heat_photo is %e\n", i,
	     plasmamain[nplasma].heat_photo);
    if (sane_check (plasmamain[nplasma].heat_ff))
      Error ("wind_update: w\[%d}.heat_ff is %e\n", i,
	     plasmamain[nplasma].heat_ff);
    if (sane_check (plasmamain[nplasma].heat_lines))
      Error ("wind_update: w\[%d}.heat_lines is %e\n", i,
	     plasmamain[nplasma].heat_lines);
    xsum += plasmamain[nplasma].heat_tot;
    psum += plasmamain[nplasma].heat_photo;
    fsum += plasmamain[nplasma].heat_ff;
    lsum += plasmamain[nplasma].heat_lines;
  }

  asum = wind_luminosity (0.0, INFINITY);
  Log
    ("!!wind_update:  Absorbed flux   %8.2e  (photo %8.2e ff %8.2e lines %8.2e)\n",
     xsum, psum, fsum, lsum);
  Log
    ("!!wind_update: Wind luminosity  %8.2e (recomb %8.2e ff %8.2e lines %8.2e) after update\n",
     asum, geo.lum_fb, geo.lum_ff, geo.lum_lines);

  /* Print out some diagnositics of the changes in the wind update */
  t_r_ave_old /= iave;
  t_e_ave_old /= iave;
  t_r_ave /= iave;
  t_e_ave /= iave;
  wind_n_to_ij (nmax_r, &i, &j);
  Log ("!!wind_update: Max change in t_r %6.0f at cell %4d (%d,%d)\n", dt_r,
       nmax_r, i, j);
  Log ("!!wind_update: Ave change in t_r %6.0f from %6.0f to %6.0f\n",
       (t_r_ave - t_r_ave_old), t_r_ave_old, t_r_ave);
  wind_n_to_ij (nmax_e, &i, &j);
  Log ("!!wind_update: Max change in t_e %6.0f at cell %4d (%d,%d)\n", dt_e,
       nmax_e, i, j);
  Log ("!!wind_update: Ave change in t_e %6.0f from %6.0f to %6.0f\n",
       (t_e_ave - t_e_ave_old), t_e_ave_old, t_e_ave);

  check_convergence ();
/* Sumarize the raditive temperatures (ksl 04 mar)*/
  xtemp_rad (w);


  return (0);
}


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis: wind_rad_init(w) zeros those portions of the wind which contain the radiation properties 
	of the wind, i.e those portions which should be set to zeroed when the structure of the 
	wind has been changed or when you simply want to start off a calculation in a known state

 Arguments:		
	WindPtr w;

Returns:
 
Description:	
	
Notes:
	The routine is called by the main program at the beginning of each ionization calculation
	cycle.  It should zero all heating and radiation induced cooling in the wind array.  Since
	cooling is recalculated in wind_update, one needs to be sure that all of the appropriate
	cooling terms are also rezeroed there as well.
	

History:
 	97jan	ksl	Coding on python began.
 	98apr	ksl	Fixed major problem with this routine, previously I had only been 
 			reinitializing the first NDIM portions of the wind.
 	98may	ksl	Added initialization for .ioniz and .recomb
 	98sept	ksl	Added initialization for heating of h, he1, he2,z. and also for the remaining
 			luminosities.
	03jul	ksl	Updated to reflect modified wind structure with heat and cooling
			for all ions
	04mar	ss	Added additional quantities to zero for macro-atoms
        04July  SS      Modified to include setting the spontaneous recombination rates.
        04 Nov  SS      Modified to include setting the recomb_simple rates (needed for cooling in
                        simple continuua).
	06may	ksl	Modified to use wmain and plasmamain, in process of redefining all of
			the structures.  Most if not all of this will end up in plasmamain
	06jul	ksl	57+ -- Added to changes to allow for a separate macro structure
	06aug	ksl	57h -- Additional changes to allow for the fact that marcomain
			is not created at all if no macro atoms.

**************************************************************/


int
wind_rad_init ()
{
  int n, i;
  int njump;


  for (n = 0; n < NPLASMA; n++)
    {
      plasmamain[n].j = plasmamain[n].ave_freq = plasmamain[n].ntot = 0;
      plasmamain[n].heat_tot = plasmamain[n].heat_ff =
	plasmamain[n].heat_photo = plasmamain[n].heat_lines = 0.0;
      plasmamain[n].heat_z = 0.0;

      plasmamain[n].lum = plasmamain[n].lum_rad = plasmamain[n].lum_lines =
	plasmamain[n].lum_ff = 0.0;
      plasmamain[n].lum_fb = plasmamain[n].lum_z = 0.0;
      plasmamain[n].nrad = plasmamain[n].nioniz = 0;

      for (i = 0; i < nions; i++)
	{
	  plasmamain[n].ioniz[i] = plasmamain[n].recomb[i] =
	    plasmamain[n].heat_ion[i] = plasmamain[n].lum_ion[i] = 0.0;
	}

      /* Next blocks added by SS Mar 2004 to zero the Macro Atom estimators. */

      /* 57h -- 0608 -- These sections actually involve enough calculations that
         they are noticeable in term sof the overall speed.  One would if possible
         like to avoid this section, since it requires the creation of macromain, 
         even though macromain is not used -- ksl */


//OLD      for (i = 0; i < NLEVELS_MACRO; i++)
      for (i = 0; i < nlevels_macro; i++)	//57h
	{
	  for (njump = 0; njump < NBBJUMPS; njump++)
	    {
	      macromain[n].jbar[i][njump] = 0.0;	// mean intensity
	    }
	  for (njump = 0; njump < NBFJUMPS; njump++)
	    {
	      macromain[n].gamma[i][njump] = 0.0;
	      macromain[n].gamma_e[i][njump] = 0.0;
	      macromain[n].alpha_st[i][njump] = 0.0;	//stimulated recombination
	      macromain[n].alpha_st_e[i][njump] = 0.0;

	      /* Next block to set spontaneous recombination rates for next iteration. (SS July 04) */

	      if (njump < config[i].n_bfd_jump && plasmamain[n].t_e > 0.0)
		{

		  //04Jul--ksl-modified these calls to reflect changed alpha_sp
		  macromain[n].recomb_sp[i][njump] =
		    alpha_sp (&phot_top[config[i].bfd_jump[njump]],
			      &plasmamain[n], 0);
		  macromain[n].recomb_sp_e[i][njump] =
		    alpha_sp (&phot_top[config[i].bfd_jump[njump]],
			      &plasmamain[n], 2);
		}
	      else
		{
		  macromain[n].recomb_sp[i][njump] = 0.0;
		  macromain[n].recomb_sp_e[i][njump] = 0.0;
		}

	    }
	}

      for (i = 0; i < ntop_phot; i++)
	{
/* 57h -- recomb_simple is only required for we are using a macro atom approach, and only non-zero when
this particular phot_tob xsection is treated as a simple x-section. Stuart, is this correct?? I've added
checks so that macro_info is only 0 (false) or true (1), and so the logic of the next section can be 
simplified.  0608-ksl */
//OLD57h          if ((phot_top[i].macro_info != 1) && (!geo.macro_simple)) // 57h
//OLD57h            {
//OLD57h              plasmamain[n].recomb_simple[i] =
//OLD57h                alpha_sp (&phot_top[i], &plasmamain[n], 2);
//OLD57h            }
//OLD57h          else
//OLD57h            {
//OLD57h              plasmamain[n].recomb_simple[i] = 0.0;
//OLD57h            }
	  if (geo.macro_simple || phot_top[i].macro_info)
	    {

	      plasmamain[n].recomb_simple[i] = 0.0;
	    }
	  else
	    {			//we want a macro approach, but not for this ion so need recomb_simple
	      plasmamain[n].recomb_simple[i] =
		alpha_sp (&phot_top[i], &plasmamain[n], 2);
	    }
	}


      //zero the emissivities that are needed for the spectral synthesis step.
      plasmamain[n].kpkt_emiss = 0.0;
      plasmamain[n].kpkt_abs = 0.0;
//OLD      for (i = 0; i < NLEVELS_MACRO; i++)
      for (i = 0; i < nlevels_macro; i++)	//57h
	{
	  macromain[n].matom_abs[i] = 0.0;

	  macromain[n].matom_emiss[i] = 0.0;

	}

      /* End of added material. */
    }


  return (0);
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis: wind_rad_summary(w,filename,mode) writes out a file which sumarizes the
	the properties of the radiation field in each cell in the wind. 
Arguments:		
	WindPtr w;
	char filename[]		Name of the output file
	char mode[]		Either "w" to write a new file or "a" to append to
						and existing file
Returns:
 
Description:	
	
Notes:
	There is not much point in calling this until you have propagated a few photons
History:
 	97jan      ksl	Coding on python began.
	06may	ksl	57+: Began modificatioins to split structures.  Much of this
			will end up in the plasma structure.

**************************************************************/

int istat_wind_rad_summary = 0;

int
wind_rad_summary (w, filename, mode)
     WindPtr w;
     char filename[], mode[];
{
  FILE *fopen (), *fptr;
  int n;
  int nplasma;


  /* Open or reopen the file where you want to put the summary of the radiation properties of the wind */
  if (istat_wind_rad_summary == 0)
    {				/*First time through recreate the file regardless */
      if ((fptr = fopen (filename, "w")) == NULL)
	{
	  Error ("wind_rad_summary: Unable to open %s\n", filename);
	  exit (0);
	}

    }
  else
    {
      if ((fptr = fopen (filename, mode)) == NULL)
	{
	  Error ("wind_rad_summary: Unable to reopen %s with mode %s\n",
		 filename, mode);
	  exit (0);
	}
    }

  istat_wind_rad_summary++;

  fprintf (fptr,
	   "n       x         z         j       ntot    ave_freq    T_rad      W       Lum     Lum_abs\n");
  for (n = 0; n < NDIM2 - 1; n++)
    {
      nplasma = w[n].nplasma;
      fprintf (fptr,
	       "%-3d %8.3e %8.3e %8.3e %8d %8.3e %8.3e %8.3e %8.3e %8.3e\n",
	       n, w[n].x[0], w[n].x[2], plasmamain[n].j, plasmamain[n].ntot,
	       plasmamain[nplasma].ave_freq, plasmamain[nplasma].t_r,
	       plasmamain[nplasma].w, plasmamain[nplasma].lum,
	       plasmamain[nplasma].heat_tot);

    }

  fclose (fptr);

  return (0);
}
