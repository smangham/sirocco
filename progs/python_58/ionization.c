

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
 	ionization routines for the wind one cell at a time
	
 Arguments:		
	WindPtr w;

Returns:
 
Description:	
	The intent is that the routine ion_abundances is the steering routine for
	all calculations of the abundances
	
Notes:

History:
	98jul	ksl	Incorporated Christians mods with my comments 
	98nov	ksl	Created the steering routine ion_abundances
	00jan	ksl	Homogenized calls so that ionization_mod_on_the_spot
			and ionization_mod_on_the_spot_exact do not call the
			entire wind structure.  NOTE THAT THESE LAST ROUTINES
			HAVE NOT BEEN SERIOUSLY RECHECKED, SINCE THEY ARE NOT
			MUCH USED.
	00jan	ksl	Moved ionization_mod_on_the_spot and ionization_mod_on_the_spot_exact
			to "Legacy code" Legacy1.c
	01oct	ksl	Add calls to levels for calculation of level populations.  The
			calls are currently hardwired.
**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"

#include "python.h"

int
ion_abundances (w, mode)
     PlasmaPtr w;
     int mode;
{
  double nh, ne;
  int ireturn;


  if (mode == 0)
    {
/*on-the-spot approximation using existing t_e.   This routine does not attempt 
 * to match heating and cooling in the wind element! */

      if ((ireturn = nebular_concentrations (w, 2)))
	{
	  Error
	    ("ionization_abundances: nebular_concentrations failed to converge\n");
	  Error
	    ("ionization_abundances: j %8.2e t_e %8.2e t_r %8.2e w %8.2e\n",
	     w->j, w->t_e, w->w);
	}
    }
  else if (mode == 1)
    {
      // LTE using t_r

      ireturn = nebular_concentrations (w, 1);
    }
  else if (mode == 2)
    {				//  Hardwired concentrations

//?? ERROR -- For consistence fix_concentrations calls should becomve fix_concentrations(w), not as currently
//?? Error --  Also, when this is done, levels can be incorporated there 
      nh = w->rho * rho2nh;
      ireturn = fix_concentrations (nh, w->density, &ne);
      w->ne = ne;
      do_partitions (w, 0);
      levels (w, 0);		// levels from reduced BB
    }
  else if (mode == 3 || mode == 5)
    {
/* On the spot, with one_shot at updating t_e before calculating densities
*/

/* Shift values to old */
      w->dt_e_old = w->dt_e;
      w->dt_e = w->t_e - w->t_e_old;	//Must store this before others
      w->t_e_old = w->t_e;
      w->t_r_old = w->t_r;
      w->lum_rad_old = w->lum_rad;

      ireturn = one_shot (w, mode);

/* Convergence check */
      convergence (w);
    }
  else if (mode == 4)
    {				// Test method with detailed balance partially included
      ireturn = nebular_concentrations (w, 4);
    }
  else
    {
      Error
	("ion_abundances: Could not calculate abundances for mode %d\n",
	 mode);
      exit (0);
    }


  return (ireturn);

}


/***********************************************************
              Space Telescope Science Institute

 Synopsis: convergence checks to see whehter a single cell
	is or is not converging
	
 Arguments:		
	WindPtr w;

Returns:
 
Description:	
	
Notes:

History:
	06may	ksl	57+: Began modifications to reflect new
			structure definition.  Note that it
			is likelely that this entire routine
			will ultimatetely be replaced because
			everything here should only be in the wind
**************************************************************/
int
convergence (xplasma)
     PlasmaPtr xplasma;
{
  int trcheck, techeck, hccheck, whole_check, converging;
  double epsilon;

  trcheck = techeck = hccheck = converging = 0;
  epsilon = 0.05;

  if ((xplasma->converge_t_r =
       fabs (xplasma->t_r_old - xplasma->t_r) / (xplasma->t_r_old +
						 xplasma->t_r)) > epsilon)
    trcheck = 1;
  if ((xplasma->converge_t_e =
       fabs (xplasma->t_e_old - xplasma->t_e) / (xplasma->t_e_old +
						 xplasma->t_e)) > epsilon)
    techeck = 1;
  if ((xplasma->converge_hc =
       fabs (xplasma->heat_tot - xplasma->lum_rad) / (xplasma->heat_tot +
						      xplasma->lum_rad)) >
      epsilon)
    hccheck = 1;

  xplasma->converge_whole = whole_check = trcheck + techeck + hccheck;

  if (xplasma->dt_e_old * xplasma->dt_e < 0
      && fabs (xplasma->dt_e) > fabs (xplasma->dt_e_old))
    converging = 1;
  xplasma->converging = converging;

  if (converging == 1)
    {				// Not converging
      xplasma->gain *= 0.7;
      if (xplasma->gain < 0.1)
	xplasma->gain = 0.1;
    }
  else
    {
      xplasma->gain *= 1.1;
      if (xplasma->gain > 0.8)
	xplasma->gain = 0.8;
    }

  return (whole_check);
}





/***********************************************************
              Space Telescope Science Institute

 Synopsis:
   check_convergence -- Do a global check on how well the wind is converging
	
 Arguments:		
	WindPtr w;

Returns:
 
Description:	
	
Notes:
Notes--Eventually should absorb some of the calculations in wind_updates

History:
	05apr	ksl	55d -- Trivially eliminated MDIM from this
			in effort to make coordinate system independent.
			Whether something is in the wind or not is
			now totally determeined by the volume.
        06may	ksl	57+ -- Now using plasma structure
**************************************************************/
int
check_convergence ()
{
  int n;
  int nconverge, nconverging, ntot;
  double xconverge, xconverging;

  nconverge = nconverging = ntot = 0;

  for (n = 0; n < NPLASMA; n++)
    {
      ntot++;
      if (plasmamain[n].converge_whole == 0)
	nconverge++;
      if (plasmamain[n].converging == 0)
	nconverging++;

    }

  xconverge = ((double) nconverge) / ntot;
  xconverging = ((double) nconverging) / ntot;
  Log
    ("!!Check_converging: %d (%.3f) converged and %d (%.3f) converging of %d cells\n",
     nconverge, xconverge, nconverging, xconverging, ntot);
  return (0);
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	ionization_on_the_spot(w) calculates densities of ions in a single element of the wind
	according to equation 11 of Lucy & Mazzali.
	
 Arguments:		
	WindPtr w;

Returns:
 
Description:
	ionization_on_the_spot uses values of t_e, t_r, etc in the windptr w to determine what
	the ion concentrations ought to be.	
	
Notes:
	This routine does not attempt to match heating and cooling in the wind element!

	In Lucy and Mazzali, they suggest setting t_e to 0.9 t_r as an approximation.  This is
	no longer enforeced.  Also, at one point CK used Cloudy in an attempt to calculate
	a relationship between w,ne,t_r, etc. His parameterization was:
	
	w->t_e=3.409*pow(trad,0.9292)*pow(nh,-0.01992)*pow(w->w,0.1317);
	but that has been removed.

History:
	97	ksl	Coded as part of python effort

**************************************************************/
// ??? ionization_on_the_spot does nothing hardly Delete ??? ksl 01dec
//int
//ionization_on_the_spot (w)
//     WindPtr w;

//{
//  double nh;
//  int nebular_concentrations ();

//  nh = w->rho * rho2nh;


//  if (w->t_r > 10.)
//    {                         /* Then modify to an on the spot approx */
//      if (nebular_concentrations (w, 2))
//      {
//        Error
//          ("ionization_on_the_spot: nebular_concentrations failed to converge\n");
//        Error
//          ("ionization_on_the_spot: j %8.2e t_e %8.2e t_r %8.2e w %8.2e\n",
//           w->j, w->t_e, w->w);
//      }
 //     if (w->ne < 0 || INFINITY < w->ne)
//      {
//        Error ("ionization_on_the_spot: ne = %8.2e out of range\n", w->ne);
//      }
//    }
//  else
//    {
//      Error ("ionization_on_the_spot: t_r exceptionally small %g\n", w->t_r);
//      mytrap ();
//      exit (0);
//    }

//  return (0);
//}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	one_shot calculates new densities of ions in a single element of the wind
	according to equation 11 of Lucy & Mazzali, after having found the
	temperature which matches heating and cooling for the previous
	densities
	
 Arguments:		
	PlasmaPtr xplasma;

Returns:
 
Description:
	
Notes:
	This routine attempts to match heating and cooling in the wind element!
	To do this it calls calc_te.  Based on the returned value of te, the
	routine then calculates the concentrations in the on-the-spot approximation.

	IT SEEMS LIKELY some code could be eliminated by simply having this routine
	call the on-the-spot routine directly.

History:
	98	ksl	Coded as part of python effort
	02jul	ksl	Added mode variable so could try detailed balance
	06may	ksl	57+ -- Switched to use plasma structue

**************************************************************/

int
one_shot (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;

{
  double te_old, te_new, dte;
  double gain;


//OLD?  wxyz = w;   // Elimated this as it did not seem used anywhere 


  gain = xplasma->gain;

  te_old = xplasma->t_e;
  te_new = calc_te (xplasma, 0.7 * te_old, 1.3 * te_old);
  xplasma->t_e = (1 - gain) * te_old + gain * te_new;

  dte = xplasma->dt_e;

//  Log ("One_shot: %10.2f %10.2f %10.2f\n", te_old, te_new, w->t_e);


/* Modes in the driving routines are not identical to those in nebular concentrations.
The next lines are an attempt to mediate this problem.  It might be better internally
at least to define a flag for using one shot, and have the modes take on the
meaning in nebular concentrations.
*/

  if (mode == 3)
    mode = 2;
  else if (mode != 4)
    {

      Error ("one_shot: Sorry, Charlie, don't know how to process mode %d\n",
	     mode);
      exit (0);
    }

  if (xplasma->t_r > 10.)
    {				/* Then modify to an on the spot approx */
      if (nebular_concentrations (xplasma, mode))
	{
	  Error
	    ("ionization_on_the_spot: nebular_concentrations failed to converge\n");
	  Error
	    ("ionization_on_the_spot: j %8.2e t_e %8.2e t_r %8.2e w %8.2e\n",
	     xplasma->j, xplasma->t_e, xplasma->w);
	}
      if (xplasma->ne < 0 || INFINITY < xplasma->ne)
	{
	  Error ("ionization_on_the_spot: ne = %8.2e out of range\n",
		 xplasma->ne);
	}
    }
  else
    {
      Error ("ionization_on_the_spot: t_r exceptionally small %g\n",
	     xplasma->t_r);
      mytrap ();
      exit (0);
    }


  return (0);
}


/* calc_te determines and returns the electron temperature in the wind such that the energy emitted
   by the wind is equal to energy emitted.

   Description:
   calc_te does not modify any abundances.  It simply takes the current value of the heating in the
   cell and attempts to find the value of the electron temperature which will result in cooling which
   matches the heating.

   This approach is not entirely self consistent because if te is actually different then the
   abundances will be different and the heating will change as well.

   This routine is a kluge because it does not really deal with what happens if the cooling curve 
   has maxima and minima.

   xxxplasma is just a way to tranmit information to zero_emit

   History:

   98dec        ksl     Updated calls so that tmin and tmax were communicated externally,
			rather than hardwired
   01dec	ksl	Reversed 98dec decision as a result of massive changes for python38
   01dec	ksl	Assured that ww->t_e is updated here
   01dec	ksl	Added capability to modify the desired goal of calc_te from the full
			heating to something intermediate between the current value and the
			ultimate goal
   01dec	ksl	Rewrote to assure that boundaries will be bracketed properly and if
			not calc_te will handle
   04June       SS      Modified so that changes in the heating rate due to changes in the
                        temperature are included for macro atoms.
	06may	ksl	Modified for plasma structue
 */

PlasmaPtr xxxplasma;

double
calc_te (xplasma, tmin, tmax)
     PlasmaPtr xplasma;
     double tmin, tmax;
{
  double heat_tot;
  double z1, z2;
  int macro_pops ();
  xxxplasma = xplasma;
  heat_tot = xplasma->heat_tot;

  xplasma->t_e = tmin;
  z1 = zero_emit (tmin);
  xplasma->t_e = tmax;
  z2 = zero_emit (tmax);

  if ((z1 * z2 < 0.0))
    {				// Then the interval is bracketed 
      xplasma->t_e = zbrent (zero_emit, tmin, tmax, 50.);
    }
  else if (fabs (z1) < fabs (z2))
    {
      xplasma->t_e = tmin;
    }
  else
    xplasma->t_e = tmax;

  /* With the new temperature in place for the cell, get the correct value of heat_tot.
     SS June  04 */

  xplasma->heat_tot -= xplasma->heat_lines_macro;
  xplasma->heat_lines -= xplasma->heat_lines_macro;
  xplasma->heat_lines_macro = macro_bb_heating (xplasma, xplasma->t_e);
  xplasma->heat_tot += xplasma->heat_lines_macro;
  xplasma->heat_lines += xplasma->heat_lines_macro;

  xplasma->heat_tot -= xplasma->heat_photo_macro;
  xplasma->heat_photo -= xplasma->heat_photo_macro;
  xplasma->heat_photo_macro = macro_bf_heating (xplasma, xplasma->t_e);
  xplasma->heat_tot += xplasma->heat_photo_macro;
  xplasma->heat_photo += xplasma->heat_photo_macro;



  return (xplasma->t_e);

}


/* This is just a function which has a zero when total energy loss is equal to total energy gain */

double
zero_emit (t)
     double t;
{
  double difference;
  double total_emission ();
  int macro_pops ();
  double macro_bb_heating (), macro_bf_heating ();


  /* NOTE - IMPORTANT

     SS May 04: I'm removing adiabatic cooling for now. I'll need to be put
     back in once the heating/cooling tests are all sorted out. 

     SS June 04: Adiabatic cooling is still not switched on. But it is now 
     switched off in emission.c rather than here.

     SS July 04: Adiabatic cooling is now switched back on. These comments can
     all be deleted once this is tested (soon).
   */

  /* This block is not needed now - SS July 04
     if (xxxplasma->lum_adiabatic > 0)
     {
     Error("zero_emit: adiabatic cooling is switched on.\n");
     }
   */

  /*Original method */
  xxxplasma->t_e = t;


  /* Correct heat_tot for the change in temperature. SS June 04. */
  //macro_pops (xxxplasma, xxxplasma->ne);
  xxxplasma->heat_tot -= xxxplasma->heat_lines_macro;
  xxxplasma->heat_lines -= xxxplasma->heat_lines_macro;
  xxxplasma->heat_lines_macro = macro_bb_heating (xxxplasma, t);
  xxxplasma->heat_tot += xxxplasma->heat_lines_macro;
  xxxplasma->heat_lines += xxxplasma->heat_lines_macro;

  xxxplasma->heat_tot -= xxxplasma->heat_photo_macro;
  xxxplasma->heat_photo -= xxxplasma->heat_photo_macro;
  xxxplasma->heat_photo_macro = macro_bf_heating (xxxplasma, t);
  xxxplasma->heat_tot += xxxplasma->heat_photo_macro;
  xxxplasma->heat_photo += xxxplasma->heat_photo_macro;

  //  difference = (xxxplasma->heat_tot - total_emission (xxxplasma, 0., INFINITY));
  difference =
    xxxplasma->heat_tot - xxxplasma->lum_adiabatic -
    total_emission (&wmain[xxxplasma->nwind], 0., INFINITY);

  return (difference);
}
