


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  Description:	

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:


 

  History:
2004	ksl	Coded as better.c
2017	nsh - all cooling mechanisms gathered here, removed from the end of
				ionization.c. 

 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


double
cooling (xxxplasma, t)
     PlasmaPtr xxxplasma;
     double t;
{

  xxxplasma->t_e = t;





  if (geo.adiabatic)
  {
    if (wmain[xxxplasma->nwind].div_v >= 0.0)
    {
      /* This is the case where we have adiabatic cooling - we want to retain the old behaviour, 
         so we use the 'test' temperature to compute it. If div_v is less than zero, we don't do
         anything here, and so the existing value of adiabatic cooling is used - this was computed 
         in wind_updates2d before the call to ion_abundances. */
      xxxplasma->cool_adiabatic = adiabatic_cooling (&wmain[xxxplasma->nwind], t);
    }
  }

  else
  {
    xxxplasma->cool_adiabatic = 0.0;
  }


  /*81c - nsh - we now treat DR cooling as a recombinational process - still unsure as to how to treat emission, so at the moment
     it remains here */

  xxxplasma->cool_dr = total_fb (&wmain[xxxplasma->nwind], t, 0, VERY_BIG, FB_REDUCED, INNER_SHELL);

  /* 78b - nsh adding this line in next to calculate direct ionization cooling without generating photons */

  xxxplasma->cool_di = total_di (&wmain[xxxplasma->nwind], t);

  /* 70g compton cooling calculated here to avoid generating photons */

  xxxplasma->cool_comp = total_comp (&wmain[xxxplasma->nwind], t);

  /* we now call xtotal emission which computes the cooling rates for processes which can, in principle, make photons. */

  xxxplasma->cool_tot =
    xxxplasma->cool_adiabatic + xxxplasma->cool_dr + xxxplasma->cool_di +
    xxxplasma->cool_comp + xtotal_emission (&wmain[xxxplasma->nwind], 0., VERY_BIG);

  return (xxxplasma->cool_tot);
}




/***********************************************************
             Space Telescope Science Institute

Synopsis:  xtotal_emission (one, f1, f2) Calculate the total cooling of a single cell 
			due to free free - line and recombination processes.

Arguments:		


Returns:
 
Description:	
	

Notes:
	xtotal_emission gives the total enery loss due to photons.  It does
	not include other coooling sources, e. g. adiabatic expansion.

	It returns the total cooling, but also stores the luminosity due
	to various types of emssion, e.g ff and  lines into the
	Plasms cells. The fb cooling calculated here is *not* equal to
	the fb lumniosity and so this value is stored in cool_rr.

	Comment: The call to this routine was changed when PlasmaPtrs
	were introduced, but it appears that the various routines 
	that were called were not changed.

History:
	2017 - copied from total_emission and modified for use in calculatingh
			cooling alone.
 
 
**************************************************************/


double
xtotal_emission (one, f1, f2)
     WindPtr one;               /* WindPtr to a specific cell in the wind */
     double f1, f2;             /* The minimum and maximum frequency over which the emission is
                                   integrated */
{
  double t_e;
  int nplasma;
  double cooling;
  PlasmaPtr xplasma;

  cooling = 0.0;
  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];

  t_e = xplasma->t_e;           // Change so calls to total emission are simpler

  if (f2 < f1)
  {
    xplasma->cool_tot = xplasma->lum_lines = xplasma->lum_ff = xplasma->cool_rr = 0;    //NSH 1108 Zero the new cool_comp variable NSH 1101 - removed
  }
  else
  {
    if (geo.rt_mode == RT_MODE_MACRO)   //Switch for macro atoms (SS)
    {
      xplasma->cool_rr = total_fb_matoms (xplasma, t_e, f1, f2) + total_fb (one, t_e, f1, f2, FB_REDUCED, OUTER_SHELL); //outer shellrecombinations
      //The first term here is the fb cooling due to macro ions and the second gives
      //the fb cooling due to simple ions.
      //total_fb has been modified to exclude recombinations treated using macro atoms.
      //Note: This the fb_matom call makes no use of f1 or f2. They are passed for
      //now in case they should be used in the future. But they could
      //also be removed.
      // (SS)
      cooling = xplasma->cool_rr;
      xplasma->lum_lines = total_bb_cooling (xplasma, t_e);
      cooling += xplasma->lum_lines;
      /* total_bb_cooling gives the total cooling rate due to bb transisions whether they
         are macro atoms or simple ions. */
      xplasma->lum_ff = total_free (one, t_e, f1, f2);
      cooling += xplasma->lum_ff;


    }
    else                        //default (non-macro atoms) (SS)
    {
      /*The line cooling is equal to the line emission */
      cooling = xplasma->lum_lines = total_line_emission (one, f1, f2);
      /* The free free cooling is equal to the free free emission */
      cooling += xplasma->lum_ff = total_free (one, t_e, f1, f2);
      /*The free bound cooling is equal to the recomb rate x the electron energy - the boinding energy - this is computed 
         with the FB_REDUCED switch */
      cooling += xplasma->cool_rr = total_fb (one, t_e, f1, f2, FB_REDUCED, OUTER_SHELL);       //outer shell recombinations


    }
  }


  return (cooling);


}



/***********************************************************
                                       Space Telescope Science Institute

Synopsis:  
	adiabatic_cooling (one, t) determines the amount of 
	adiabatic cooling in a cell, in units of luminosity.

Arguments:		
	WindPtr one;	pointer to wind cell
	double t;		electron temperature

Returns:
 
Description:	
   Adiabatic cooling is clearly the amount of PdV work done by
   a fluid element, per unit time dt. Thus it is equal to P dV/dt.  
   The only real question is whether dV/dt is given by the volume * div v.
   div v here is the divergence of the velocity field.
	
Notes:
  JM 1401 -- I've rederived the expression for dv/dT. It follows
        directly from the continuity equation and is indeed equal 
        to volume * div_v. 

        Note also that this function should only be called
        if geo.adiabatic == 1, in which case it populates
        xplasma->cool_adiabatic. This is used in heating and cooling
        balance. We also use it as a potential destruction choice for 
        kpkts in which case the kpkt is thrown away by setting its istat 
        to P_ADIABATIC.


History:
	04nov	ksl	Stuart had switched adiabatic cooling off
  			by commenting out the entire program.  I have
  			reinstated the calculation, but created a new
  			variable geo.adiabatic to determine whether this
  			routine is called in python.  At present, even
  			if the routine is called it has no effect.  The
  			or more properly, one issue is how to meld adiabatic 
  			cooling with the macro atom approach.
	06may	ksl	57+ -- Adapted to include plasma structure
	11aug	ksl	70 - Adiabatic cooling returns the cooling but
			does not store it.
 
 
**************************************************************/


double
adiabatic_cooling (one, t)
     WindPtr one;
     double t;
{
  double cooling;
  int nplasma, nion;
  double nparticles;
  PlasmaPtr xplasma;

  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];

  //JM 1401 -- here was an old factor of 3/2 which KSL and JM believe to be incorrect. 
  //JM 1601 -- this should include the pressure from all particles, rather than just ne
  //cooling = xplasma->ne * BOLTZMANN * t * xplasma->vol * one->div_v;
  nparticles = xplasma->ne;

  for (nion = 0; nion < nions; nion++)
  {
    /* loop over all ions as they all contribute to the pressure */
    nparticles += xplasma->density[nion];
  }

  cooling = nparticles * BOLTZMANN * t * xplasma->vol * one->div_v;

  return (cooling);
}


/***********************************************************
                                       Space Telescope Science Institute

Synopsis:  wind_cooling (f1, f2) calculate the cooling rate of the entire 
wind between freqencies f1 and f2 

Arguments:		


Returns:
 
Description:	
	

Notes:

History:
	06may	ksl	57+ -- Changed to accommodate plasma structure.  
			Note that the call for wind_luminostiy has changed
			to remove the pointer call.  Ultimately rewrite
			using plasma structure alone.
	11aug	nsh	70 Modifications made to incorporate compton cooling
    11sep   nsh     70 Modifications in incorporate DR cooling (very approximate at the moment)
	12sep	nsh	73 Added a counter for adiabatic luminosity (!)
 	13jul	nsh	76 Split up adiabatic luminosity into heating and cooling.
    17aug   nsh xchanges to wind cooling to reflect the way we now
 
**************************************************************/

double
wind_cooling (f1, f2)
     double f1, f2;             /* freqmin and freqmax */
{
  double cool, lum_lines, cool_rr, lum_ff, cool_comp, cool_dr, cool_di, cool_adiab, heat_adiab; //1108 NSH Added a new variable for compton cooling 1408 NSH and for DI cooling
  //1109 NSH Added a new variable for dielectronic cooling
  //1307 NSH Added a new variable to split out negtive adiabatic cooling (i.e. heating).
  int n;
  double x;
  int nplasma;


  cool = lum_lines = cool_rr = lum_ff = cool_comp = cool_dr = cool_di = cool_adiab = heat_adiab = 0;    //1108 NSH Zero the new counter 1109 including DR counter 1408 and the DI counter
  for (n = 0; n < NDIM2; n++)
  {

    if (wmain[n].vol > 0.0)
    {
      nplasma = wmain[n].nplasma;
      cool += x = cooling (&plasmamain[nplasma], plasmamain[nplasma].t_e);      //1708 - changed this call - now computes cooling rather than luminosity - also we popultate a local
      //array called cool, rather than the xplasma array - this was overrwting the xplasma array with incorrect data. 
      lum_lines += plasmamain[nplasma].lum_lines;
      cool_rr += plasmamain[nplasma].cool_rr;
      lum_ff += plasmamain[nplasma].lum_ff;
      cool_comp += plasmamain[nplasma].cool_comp;       //1108 NSH Increment the new counter by the compton luminosity for that cell.
      cool_dr += plasmamain[nplasma].cool_dr;   //1109 NSH Increment the new counter by the DR luminosity for the cell.
      cool_di += plasmamain[nplasma].cool_di;   //1408 NSH Increment the new counter by the DI luminosity for the cell.

      if (geo.adiabatic)        //130722 NSH - slight change to allow for adiabatic heating effect - now logged in a new global variable for reporting.
      {

        if (plasmamain[nplasma].cool_adiabatic >= 0.0)
        {
          cool_adiab += plasmamain[nplasma].cool_adiabatic;
        }
        else
        {
          heat_adiab += plasmamain[nplasma].cool_adiabatic;
        }
      }

      else
      {
        cool_adiab = 0.0;
      }


      if (x < 0)
        Error ("wind_cooling: xtotal emission %8.4e is < 0!\n", x);

      if (recipes_error != 0)
        Error ("wind_cooling: Received recipes error on cell %d\n", n);
    }
  }

/* Deciding which of these to fill is tricky because geo contains
   values for geo.lum_wind and geo.f_wind separately.  Could be cleaner.
   ksl 98mar6 */

  //geo.lum_wind=lum;
  geo.lum_lines = lum_lines;
  geo.cool_rr = cool_rr;
  geo.lum_ff = lum_ff;
  geo.cool_comp = cool_comp;    //1108 NSH The total compton luminosity of the wind is stored in the geo structure
  geo.cool_dr = cool_dr;        //1109 NSH the total DR luminosity of the wind is stored in the geo structure
  geo.cool_di = cool_di;        //1408 NSH the total DI luminosity of the wind is stored in the geo structure
  geo.cool_adiabatic = cool_adiab;
  geo.heat_adiabatic = heat_adiab;

//  cool = cool+ cool_comp+cool_dr+cool_di+cool_adiab; //1708 NSH we no longer need to add these things on - its done in cooli ng

  return (cool);
}
