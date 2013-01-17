

/* This file contains various ways in which to calculate the
   ionization of a plasma.  The guiding difference between
   whether to put routines here or in ionization.c appears to
   be that here one takes the conditions, radiative weight as
   given and calculate the densities of each ionic species, while
   ionization.c is involved in trying to match heating and cooling.

	00Jan02	ksl	Removed some of the routines which allowed
   			special calculation of H and He 
	01dec	ksl	Major modifications to all routines.  (a) update
			calls to concentrations and nebular_concentrations
			so WindPtrs are passed.
	02jun	ksl	Modified calls yet again.  This time the goal
			is to have all the modes etc interpreted in
			nebular concentrations, so one can actually tell what 
			is happening in nebular concentrations, and to
			avoid doing any real work in nebular concentrations.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "atomic.h"
#include "python.h"

/***********************************************************
                                       Space Telescope Science Institute

  Synopsis:   

int
nebular_concentrations (ww, mode)  modifies the densities of ions, levels, and
	partition functions of ions within a cell of the wind based upon the mode, 
	and other data contained within the WindPtr itself, such as t_r, t_e, w, 
	based on the "mode".  

        This version of the code is a first, quick and dirty, attempt to include the
        fraction of recombinations going into the ground state as a correction factor
        following Mazzali and Lucy (CK).  In other words, ionization equilibrium here
        is calculated from equation 11 of Mazzali & Lucy.
  
  Arguments:		
     WindPtr ww;
     int mode;			// 0=saha using tr, 1=saha using te, 2= Lucy & Mazzali


  Returns:
 	nebular concentrations alters portions of the wind ptr.  Exactly how things 
	are changed depends on the mode.
 
 	nebular_concentrations returns 0 if it converged to an answer, -1 otherwise.  On an
 	abnormal return the density array and ne are left unchanged.
 	
 Description:	

	As of 01dec (python39), nebular concentrations serves as the steering routine
	for all ionization calculations.  This change was made to push down all of the
	various choices about calculation of ionization to a single routine, so that
	no longer needs to include case statements or logic for how to handle the various
	choices elsewhre.  The Lucy & Mazzali calculation is still carried out within
	nebular concentrations although that is historical at present.

	It remains true that for the Lucy and Mazalli calculation, concentrations is
	called initially, a choice which today seems bizarre...since the LM calculation
	can be carried equally well straight away.
	
	My pre-01dec notes said:

	This is a first attempt to fudge the ionization equation in the wind to account 
	for the fact that one is not sitting in a BB radiation field.  The routine assumes one
	has already put the concentrations of of the ions under LTE conditions for a temperature
	t_r into density.  It uses the same dumb procedure as for "concentrations."  

 
  Notes:
	It's probably OK but it would be worthwhile worrying whether I have introduced 
	any problems by setting the floor to any species and ne of DENSITY_MIN	

	Mode is not consistent across all of the parts of python, so don't assume that
	mode here is geo.ionization_mode...for example.


  History:
	1997	ksl	Coded and debugged as part of Python effort. 
	97aug27	ksl	Eliminated the goto statement in the main iteration loop and
 			added numerous comments.
	97oct2	ck	Implemented Mazzali & Lucy eqn 11.
	98may	ksl	Moved ne calculation to separate subroutine
	98may	ksl	Moved actual calculation for Mazzali & Lucy to a separate routine for
 			each element
	01dec	ksl	Changed the way in which nebular concentrations was called, so that all
			the variables are taken from the windptr except for the mode. In addition
			nebular_concentrations was made the sole call from other parts of python, and
			so for example saha should never be called elsewhere.  This last step wasw
			to encapsulate the determination of abundances (given known conditions.)
	02jun	ksl	Moved Lucy-Mazzali calculation to a separate routine.  Now nebular
			concentrations is simply a driving routine
	06may	ksl	57+ -- Modified to use plasma structure.  There are no volumes involved
			so use plasma stucture in call to routines.


**************************************************************/


int
nebular_concentrations (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;			// 0=saha using tr, 1=saha using te, 2= Lucy & Mazzali
{
  double get_ne ();
  int lucy_mazzali1 ();
  int m;
  int concentrations (), lucy (), dlucy ();



  if (mode == 0)
    {				// LTE all the way -- uses tr
      do_partitions (xplasma, mode);
      levels (xplasma, mode);
      m = concentrations (xplasma, mode);
      return (m);

    }
  else if (mode == 1)
    {				//LTE but uses temperature te
      do_partitions (xplasma, mode);
      levels (xplasma, mode);
      m = concentrations (xplasma, mode);
      return (m);

    }
  else if (mode == 2)		// This is the standard LM method
    {

      do_partitions (xplasma, mode);
      levels (xplasma, mode);
      m = concentrations (xplasma, 0);	// Saha equation using t_r

      m = lucy (xplasma);	// Main routine for running LucyMazzali

      return (m);

    }
  else if (mode == 4)
    {
      /* This is a new method which tries to incorprate some detailed
       * calculations based on ionization rates */
      do_partitions (xplasma, 2);
      levels (xplasma, 2);

      if (xplasma->nioniz > 0)
	m = dlucy (xplasma);	// Main routine for running detailed balance
      else
	{
	  m = concentrations (xplasma, 0);	// Saha equation using t_r
	  m = lucy (xplasma);	//Only run detailed balance if one already has ionizing photons 
	  // version of LucyMazzali
	}

      return (m);

    }


// If reached this point the program does not understand what is desired.
  Error ("nebular_concentrations: Unknown mode %d\n", mode);
  exit (0);

}


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:   
int
concentrations (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;			//   0=saha using tr, 1=saha using te
  	determines the concentrations of individual ions assuming a form of the 
	Saha equantion.  The atomic data required must have been read in with the 
	routine "get_atomic_data".  Most of the actual information to calculate
	the concentrations is already contained in ww
  
 Arguments:		

     WindPtr ww;
     int mode;			//   0=saha using tr, 1=saha using te

 Returns:
	concentrations returns 0 if it converged; -1 if it did not converge.  If it does
	not converge, both density and *ne are likely to be garbage.

	Other results...e.g. the densities.. are stored in ww.
 
Description:	
 	
	The Saha eqn. is  (e.g. Motz, p 110, eqn 4.14

  	 N_r+1 N_e / N_r =  (2*pi*m*kT)**1.5/h**3 2 Z_r+1(T)/Z_r * exp (-X/kT)
   
	where N_r+1 is the number density of the r+1 ionization state
	and   Z_r+1 is the partition function of that state.  

	The program works by perfoming a crude iteration on ne.  It could undoubtedly 
	be improved, if one wanted better accuracy or higher speed.


 
Notes:
	
	So that one doesn't get zero-divides in other programs, e.g. nebular_concentrations, 
	a floor is set on the density of any ion or the electron density.



History:
 	1997      ksl	Coded and debugged as part of Python effort.  All of the partition functions
 				have been set to 1.
 	97aug27 ksl	Eliminated the goto statement in the iteration of the Saha equation and
 				added numberous comments.
 	97aug27 ksl      Updated code to include an estimate of the partition function based on
 				the ground state multiplities.  Add an error return if the function failed
 				to return.
 	98may	ksl	Moved ne calculation to separate subroutine
 	98may	ksl	Moved saha calculation for an individual element to a separate routine
	00mar	ksl	Fixed error in initial estimate of ne so that it is correct at high temperatures
			Basically the problem was due to the fact that at high enough temperatures
			a subtraction was being performed which could lead to an initial value of
			0 for ne from which the code never recovered.  Note--I believe there still
			may be a small problem with this section of the code because it does not
			take account of the partition function in the initial estimate.  On the
			other hand the effect of this error should be minor.
	01jul	ksl	Added a real partition function to routine
	01dec	ksl	Changed variables used to call saha, and made a subroutine of nebular
			concentrations (python_39)
	02jul	ksl	Included the electron multiplicity in the SAHA constant.
      
        04Apr   SS      Added a call to macro_pops so that macro atom ionization fractions are used
                        in place of others.
	06may	ksl	57+ -- Switched to plasma structue wind no volumes
**************************************************************/



#define SAHA 4.82907e15		/* 2* (2.*PI*MELEC*k)**1.5 / h**3  (Calculated in constants) */
#define MAXITERATIONS	20
#define FRACTIONAL_ERROR 0.03
#define THETAMAX	 1e4
#define MIN_TEMP         100.


int
concentrations (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;			//   0=saha using tr, 1=saha using te
{
  int nelem, nion, niterate;
  double xne, xxne, xnew, xsaha;
  double theta, x;
  double get_ne ();
  double t, nh;
  int saha ();


  // This needs to be moved up into nebular_concentrations given that we
  // do not want these routines to do mode shifting.
  //
  if (mode == 0)
    {
      t = xplasma->t_r;
    }
  else if (mode == 1)
    {
      t = xplasma->t_e;
    }
  else if (mode == 2)
    {
      t = sqrt (xplasma->t_e * xplasma->t_r);
    }
  else
    {
      Error ("Concentrations: Unknown mode %d\n", mode);
      mytrap ();
      exit (0);
    }

  nh = xplasma->rho * rho2nh;	//LTE

  /* make an initial estimate of ne based on H alone,  Our guess
     assumes ion[0] is H1.  Note that x below is the fractional
     ionization of H and should vary from 0 to 1
   */

  if (t < MIN_TEMP)
    t = MIN_TEMP;		/* fudge to prevent divide by zeros */

  xsaha = SAHA * pow (t, 1.5);

  theta = xsaha * exp (-ion[0].ip / (BOLTZMANN * t)) / nh;

  if (theta < THETAMAX)
    {
      x = (-theta + sqrt (theta * theta + 4 * theta)) / 2.;
      xne = xxne = x * nh;
    }
  else
    xne = xxne = nh;

  if (xne < 1.e-6)
    xne = 1.e-6;		/* fudge to assure we can actually calculate
				   xne the first time through the loop */

  /* At this point we have an intial estimate of ne. */


  niterate = 0;
  while (niterate < MAXITERATIONS)
    {

      /* Assuming a value of ne calculate the relative densities of each ion for each element */
      for (nelem = 0; nelem < nelements; nelem++)
	{
	  saha (xne, nh, t, nelem, xplasma->density);
	}

      /* New (SS Apr 04) call the routine that will replace the macro atoms ionization fractions and
         level populations with those computed with the Monte Carlo estimators. geo.macro_iniz_mode
         controls the use of macro_pops. */

      if (geo.macro_ioniz_mode == 1)
	{
	  macro_pops (xplasma, xne);
	}

      /*Set some floor so future divisions are sensible */
      for (nion = 0; nion < nions; nion++)
	{
	  if (xplasma->density[nion] < DENSITY_MIN)
	    xplasma->density[nion] = DENSITY_MIN;
	}

      /* Now determine the new value of ne from the ion abundances */
      xnew = get_ne (xplasma->density);	/* determine the electron density for this density distribution */

      if (xnew < DENSITY_MIN)
	xnew = DENSITY_MIN;	/* fudge to keep a floor on ne */

      if (fabs ((xne - xnew) / (xnew)) < FRACTIONAL_ERROR || xnew < 1.e-6)
	break;

      xne = (xnew + xne) / 2.;	/* Make a new estimate of xne */
      niterate++;
    }

  if (niterate == MAXITERATIONS)
    {
      Error ("concentrations: failed to converge t %.2g nh %.2g xnew %.2g\n",
	     t, nh, xnew);
      Error ("concentrations: xxne %e theta %e\n", xxne, theta);
      return (-1);
    }

  xplasma->ne = xnew;
  return (0);
}

/* 
   Calculate the saha densities for a single element
   98may        ksl     Abstracted from concentrations
   01jul	ksl	Added partition functions more properly
 */
int
saha (ne, nh, t, nelem, density)
     double ne, nh, t;
     int nelem;
     double density[];
{
  int first, last, nion;
  double xsaha;
  double sum, a, b;
  double partition ();
  double big;

  xsaha = SAHA * pow (t, 1.5);
  first = ele[nelem].firstion;	/*first and last identify the position in the array */
  last = first + ele[nelem].nions;	/*  So for H which has 2 ions, H1 and H2, first will generally
					   be 0 and last will be 2 so the for loop below will just be done once for nion = 1 */
  if (first < 0 || first >= nions || last < 0 || last > nions)
    {
      Error
	("saha: Confusion for element %d with first ion %d and last ion %d\n",
	 nelem, first, last);
      exit (0);
    }

  sum = density[first] = 1.;
  big = pow (10., 250. / (last - first));

  for (nion = first + 1; nion < last; nion++)
    {
      b = xsaha * partition (nion, ne, t, 1., 1)
	* exp (-ion[nion - 1].ip / (BOLTZMANN * t)) / (ne *
						       partition (nion - 1,
								  ne, t, 1.,
								  1));
      if (b > big)
	b = big;		//limit step so there is no chance of overflow

      a = density[nion - 1] * b;

      sum += density[nion] = a;
      if (density[nion] < 0.0)
	mytrap ();
      sane_check (sum);
    }

  a = nh * ele[nelem].abun / sum;
  for (nion = first; nion < last; nion++)
    {
      density[nion] *= a;
      sane_check (density[nion]);
    }


  return (0);
}

/* lucy finishes the execution a lucy-mazalli concentration scheme.

Notes: 


Concentrations should have been called before this routine is executed.

Procedurally, the routine is analogous to concentrations()

History:

	02jun	ksl	Made separate routine, removing it from nebular concentrations

        04Apr   SS      If statement added to stop this routine from changing macro atom
	                populations.
        04May   SS      "If" statement modified for compatibility with the "macro_simple" option
	                (i.e. all ions treated as simple).

*/

#define MIN_FUDGE  1.e-10
#define MAX_FUDGE  10.

int
lucy (xplasma)
     PlasmaPtr xplasma;
{
  int nelem, nion, niterate;
  double xne, xnew;
  double newden[NIONS];
  double t_r, nh;
  double t_e, www;

  t_r = xplasma->t_r;
  t_e = xplasma->t_e;
  www = xplasma->w;


  /* Initally assume electron density from the LTE densities */

  xne = xplasma->ne;
  if (xne < DENSITY_MIN)
    {
      Log
	("Warning:nebular_concentrations: Very low ionization: ne initially %8.2e\n",
	 xne);
      xne = DENSITY_MIN;
    }

  nh = xplasma->rho * rho2nh;	//LTE -- Not clear needed at this level

  /* Begin iteration loop to find ne */
  niterate = 0;
  while (niterate < MAXITERATIONS)
    {
      for (nelem = 0; nelem < nelements; nelem++)
	{
	  lucy_mazzali1 (nh, t_r, t_e, www, nelem, xplasma->ne,
			 xplasma->density, xne, newden);
	}
      for (nion = 0; nion < nions; nion++)
	{
	  /*Set some floor so future divisions are sensible */
	  if (newden[nion] < DENSITY_MIN)
	    newden[nion] = DENSITY_MIN;
	}
      xnew = get_ne (newden);
      if (xnew < DENSITY_MIN)
	xnew = DENSITY_MIN;

      /* Check to see whether the search for xne has converged and if so exit loop */
      if (fabs ((xne - xnew) / (xnew)) < FRACTIONAL_ERROR || xnew < 1.e-6)
	break;

      /* else star another iteration of the main loop */
      xne = (xnew + xne) / 2.;	/* Make a new estimate of xne */
      niterate++;
    }
  /* End of main iteration loop */

  if (niterate == MAXITERATIONS)
    {
      Error
	("nebular_concentrations: failed to converge:nh %8.2e www %8.2e t_e %8.2e  t_r %8.2e \n",
	 nh, www, t_e, t_r);
      return (-1);
    }

/* Finally transfer the calculated densities to the real density array */

  xplasma->ne = xnew;
  for (nion = 0; nion < nions; nion++)
    {
      /* If statement added here to suppress interference with macro populations (SS Apr 04) */
      if (ion[nion].macro_info == 0 || geo.macro_ioniz_mode == 0
	  || geo.macro_simple == 1)
	{
	  xplasma->density[nion] = newden[nion];
	}
    }
  return (0);
}

/***********************************************************
                                       Space Telescope Science Institute

  Synopsis:   
	lucy_mazzali1(nh,t_r,t_e,www,nelem,ne,density,xne,newden) calculates ion
	densities of a single element (defined by nelem) according to a modified on
	the spot aprroximation, in this case equation 11 of Mazzali and Lucy which
	attempts to consider ionzations from excieted and ground state levels.
 
  Arguments:		

	double nh,ne		number density of H atoms and the electron density
	double t_r,t_e		temperature of the radiation field and the electrons in degrees
	double density[]	densities of individual ions under pure Saha conditions
	double www		dilution factor for the radiation field
	int nelem;		the element number for which the calculation is carried out

  Returns:
 	double newden[]		the modified concentrations of the individual ions. Array 
 				elements of density correspond to those ions in the structure ions.
 
 	nebular_concentrations returns 0 if it converged to an answer, -1 otherwise.  On an
 	abnormal return the density array and ne are left unchanged.
 	
 Description:	

	This is a first attempt to fudge the ionization equation in the wind to account 
	for the fact that one is not sitting in a BB radiation field.  The routine assumes one
	has already put the concentrations of of the ions under LTE conditions for a temperature
	t_r into density.  It uses the same dumb procedure as for "concentrations."  

 
  Notes:
	It's probably OK but it would be worthwhile worrying whether I have introduced 
	any problems by setting the floor to any species and ne of DENSITY_MIN	


  History:
 
 98may	ksl	Coded as a a separate routine for each element

**************************************************************/


int
lucy_mazzali1 (nh, t_r, t_e, www, nelem, ne, density, xne, newden)
     double nh, t_r, t_e, www;
     int nelem;
     double ne, density[], xne, newden[];
{
  double fudge, dummy, interpfrac;
  double fudge2, q;
  double sum, a;
  int ilow, ihi;
  int first, last, nion;
  double numerator, denominator;

  if (t_r > MIN_TEMP)
    {
      fudge = www * sqrt (t_e / t_r);

      /* now get the right place in the ground_frac tables  CK */
      dummy = t_e / TMIN - 1.;
      ilow = dummy;		/* have now truncated to integer below */
      ihi = ilow + 1;		/*these are the indeces bracketing the true value */
      interpfrac = (dummy - ilow);	/*this is the interpolation fraction */
      if (ilow < 0)
	{
	  ilow = 0;
	  ihi = 0;
	  interpfrac = 1.;
	}
      if (ihi > 19)
	{
	  ilow = 19;
	  ihi = 19;
	  interpfrac = 1.;
	}

    }
  else
    {
      Error ("lucy_mazzali1: t_r too low www %8.2e t_e %8.2e  t_r %8.2e \n",
	     www, t_e, t_r);
      fudge = 0.0;
      interpfrac = 0.0;
      ihi = ilow = 0;
    }

  if (fudge < MIN_FUDGE || MAX_FUDGE < fudge)
    {
      Log
	("Warning: lucy_mazzali1: fudge>10 www %8.2e t_e %8.2e  t_r %8.2e \n",
	 www, t_e, t_r);
    }

  /* Initialization of fudges complete */

  first = ele[nelem].firstion;	/*first and last identify the postion in the array */
  last = first + ele[nelem].nions;	/*  So for H which has 2 ions, H1 and H2, first will generally
					   be 0 and last will be 2 so the for loop below will just be done once for nion = 1 */


  while (density[first] < 1.1 * DENSITY_MIN)
    {
      newden[first] = DENSITY_MIN;
      first++;
    }
  while (density[last - 1] < 1.1 * DENSITY_MIN)
    {
      newden[last - 1] = DENSITY_MIN;
      last--;
    }

  sum = newden[first] = 1.;
  for (nion = first + 1; nion < last; nion++)
    {
      numerator = newden[nion - 1] * fudge * (ne) * density[nion];
      denominator = density[nion - 1] * xne;
      q = numerator / denominator;

      /* now apply the Mazzali and  Lucy fudge, i.e. zeta + w (1-zeta)
         first find fraction of recombinations going directly to
         the ground state for this temperature and ion */
      fudge2 =
	ground_frac[nion - 1].frac[ilow] +
	interpfrac * (ground_frac[nion - 1].frac[ihi] -
		      ground_frac[nion - 1].frac[ilow]);
      /*Since nion-1 points at state i-1 (saha: n_i/n_i-1) we want ground_frac[nion-1].
         Note that we NEVER access the last ion for any element in that way
         which is just as well since you can't have recombinations INTO
         a fully ionized state -- The corresponding lines in recombination.out 
         are just dummies to make reading the data easier */
      fudge2 = fudge2 + www * (1.0 - fudge2);
      newden[nion] = fudge2 * q;
      sum += newden[nion];
// This is the equation being calculated
//              sum+=newden[nion]=newden[nion-1]*fudge*(*ne)*density[nion]/density[nion-1]/xne;
    }
  a = nh * ele[nelem].abun / sum;
  for (nion = first; nion < last; nion++)
    newden[nion] *= a;



  return (0);
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:   
  	fix_concentrations(nh,density,ne) sets the concentrations of individual ions
  	assuming hard coded ionization fractions
  
 Arguments:		

	double nh		number density of H atoms

 Returns:
 	double density[]	the concentrations of the individual ions. Array elements of 
 				density correspond to those ions in the structure ions.
 
	double *ne		the derived electron density
	
Description:	
 	

Notes:
	To change the ions that you want to fix, simply update con_force,
	AND MAKE SURE NFORCE IS SET TO THE NUMBER OF IONS

History:
 	97sep13	ksl	Adapted from concentrations
 	98may	ksl	Moved ne calculation to separate subroutine
	01mar	ksl	Moved specification of fixed concentrations to a file
	04dec	ksl	54e -- previously a specific file fixed.dat was
			read.  With this change, the concentration file
			name is read as part of the inputs

**************************************************************/

#define LINELENGTH	132
int fix_con_start = 0;
struct force_con
{
  int z, istate;
  float frac;
}
con_force[10];
int nforce;

int
fix_concentrations (nh, density, ne)
     double nh, density[];
     double *ne;
{
  int nelem, nion;
  int n;
  double get_ne ();
  FILE *fopen (), *cptr;
  char line[LINELENGTH];


  /* Define the element and ion which will be present in the wind */


  if (fix_con_start == 0)
    {
      if ((cptr = fopen (geo.fixed_con_file, "r")) == NULL)
	{
	  Error
	    ("fix_concentrations: Could not open %s to read concentrations\n",
	     geo.fixed_con_file);
	  exit (0);
	}

      nforce = 0;
      while (fgets (line, LINELENGTH, cptr) != NULL && nforce < 10)
	{
	  if ((n =
	       sscanf (line, "%d %d %f", &con_force[nforce].z,
		       &con_force[nforce].istate,
		       &con_force[nforce].frac)) != 3)
	    {
	      Error ("fix_concentrations: Read %d items from %s\n", n, line);
	    }
	  Log ("Fixing element %d istate %d to fraction %f\n",
	       con_force[nforce].z, con_force[nforce].istate,
	       con_force[nforce].frac);
	  nforce++;
	}

      fix_con_start = 1;
      fclose (cptr);
    }

  /* Set all the ion abundances to 0 and *ne to 0 */
  *ne = 0;
  for (nion = 0; nion < nions; nion++)
    density[nion] = 0;
  /* Search for matches and if one is found set the density of that ion to be the density of that
     atom */
  for (n = 0; n < nforce; n++)
    {
      nion = 0;
      while (nion < nions &&
	     !(ion[nion].z == con_force[n].z
	       && ion[nion].istate == con_force[n].istate))
	nion++;
      if (nion < nions)
	{			/* There was a match */
	  /* Find the element associated with the ion */
	  nelem = 0;
	  while (nelem < nelements && ele[nelem].z != con_force[n].z)
	    nelem++;
	  /* Increment the ion density and the electron density */
	  density[nion] = nh * ele[nelem].abun * con_force[n].frac;
	}
    }

  *ne = get_ne (density);
  return (0);
}


double
get_ne (density)
     double density[];
{
  int n;
  double ne;
  ne = 0;
  for (n = 0; n < nions; n++)
    {
      ne += density[n] * (ion[n].istate - 1);
    }
  return (ne);
}


/* do_partitions is a driving routine for partition.  Mode here is identical
to that used by nebular_concentrations, e.g

0 = LTE using t_r
1 = LTE using t_e
2 = Lucy and Mazzali

The modes here are not identical those in partition!!!

	06may	ksl	57+	Switched to using PlasmaPtr
 */

int init_partition = 0;

int
do_partitions (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;
{
  int n, pmode;
  double partition ();
  double t, weight;

  if (mode == 0)
    {				// LTE using t_r
      pmode = 1;
      t = xplasma->t_r;
      weight = 1;
    }
  else if (mode == 1)
    {				// LTE using t_e
      pmode = 1;
      t = xplasma->t_e;
      weight = 1;
    }
  else if (mode == 2)
    {				// Non LTE calculation with radiative weights
      pmode = 1;
      t = xplasma->t_r;
      weight = xplasma->w;
    }
  else
    {
      Error ("do_partitions: Unknown mode %d\n", mode);
      exit (0);
    }

  for (n = 0; n < nions; n++)
    {
      xplasma->partition[n] = partition (n, xplasma->ne, t, weight, pmode);
    }

  return (0);
}

/*
Calculate the partition function for this element.  

If there is no level info, we use the ground state multiplicity

This is the "improved version" of the code.

pmodes are as follows:
	0	Use ground state multiplicity
	1       Non-LTE (or if w=1 LTE) calculatation with internal levels and failing 
		that ion g's ground state multiplicities
	2	Use topbase (Hubeny's subroutine) to calculate the partition 
		function (this is LTE only).  
		Fall back to levels and ground state if necessary

Important! pmodes are not the same as the modes used by do_partitions

NB: Almost an identical calculation occurs in levels.  It is likely that partition
should subsume "levels". In the mean time, problems with one have to be fixed
in both places !! 01dec 

NB:  There is an implicit assumption that the very first level of an ion
is the ground state for that ion

At present, the partition function calculated rely on the "lte" levels.  This
probably should be rewritten with case statements or something to make it 
easier to use any option, and to have more graceful degradation with less
and less data.  

At present the probably should remove Hubeny/Topbase routine 
altogether should probably not be used.  As written it is incompatible
with most of the rest of python.  It would be be possible to use the topbase
routine to make a weighted blackbody assumption

partition (t=whatever) =  opfrac (t=0) + w  (opfrac (t=whatever) - opfrac( t=0))

I haven't done this, because it is likely to make it harder to match heating and
cooling but...it might be the right thing to do.
	

History:
	01jul	ksl	Coded and debugged
	01aug	ksl	Original and limited version of partition function
			replaced with Opacity Project Routine.
	01dec	ksl	Recoded to favor internally calculated g after recognition
			that the topbase partition functions were actually only
			useful in an LTE case.  This is because one cannot
			obtain detailed balance for fb photons in a reduced
			radiation field with the LTE g's
	01dec12	ksl	Modified to account for split of so-called non-LTE and 
			LTE levels
	02may	ksl	Added a kluge to get around normal setup of partition modes.
			If some of these are adopted, either need additional
			setups in do_partitions or better approach to ionization
			calls, which are needlessly complex at this point.


*/

double
partition (nion, ne, t, w, pmode)
     int nion;
     double ne, t, w;		// electron density, temperature, weight
     int pmode;

{
  int n, m;
  int m_ground;			//added by SS Jan 05
  double z, kt;
  double u, xne, xt, frac;
  int istate, iat;
  int opfrac_ ();		// returns 0 if the partition function can be calculated, -1 otherwise

  if (init_partition == 0)
    {
      iat = 0;
      opfrac_ (&iat, &istate, &xt, &xne, &u, &frac);
      init_partition = 1;
    }

  iat = ion[nion].z;
  istate = ion[nion].istate;
  xne = ne;
  xt = t;

/* Next step is goes around pmode call from do_partitions and uses geo.partition_mode instead.  Normally
geo.partition_mode should be set to -1, and in this case the partition function will be calculated in
the normal way, i.e. as the routine is called.  But one can for test reason want to force use of either
opfrac or the ground state multiplicities.  
*/
  if (geo.partition_mode > -1)
    {
      pmode = geo.partition_mode;
    }
  if (pmode == 2 && opfrac_ (&iat, &istate, &xt, &xne, &u, &frac) == 0)	// Calculate partition function using Ivan's routine
    {
      z = u;
    }

  else if (pmode == 1 && ion[nion].nlevels > 0)	// Calculate data on levels using a weighed BB assumption
    {
      kt = BOLTZMANN * t;
      m = ion[nion].firstlevel;
      m_ground = m;		//store ground state - in case energy neq 0 (SS)
      z = config[m].g;		//  Makes explicit assumption that first level is ground

      for (n = 1; n < ion[nion].nlevels; n++)
	{
	  m++;
	  z +=
	    w * config[m].g * exp ((-config[m].ex + config[m_ground].ex) /
				   kt);
	}
    }
  else if (pmode == 1 && ion[nion].nlte > 0)	// Calculate using "non-lte" levels
    {
      kt = BOLTZMANN * t;
      m = ion[nion].first_nlte_level;
      m_ground = m;		//store ground state - in case energy neq 0 (SS)
      z = config[m].g;		// This statement makes an explicit assumption that first level is ground

      for (n = 1; n < ion[nion].nlte; n++)
	{
	  m++;
	  z +=
	    w * config[m].g * exp ((-config[m].ex + config[m_ground].ex) /
				   kt);
	}
    }
  else
    {
      z = ion[nion].g;
    }

  return (z);
}

/* dlucy is a new version of the lucy scheme which allows for some detailed
 * ionization equilibrium calculations

Notes: 


Concentrations should have been called before this routine is executed.

Procedurally, the routine is analogous to concentrations()

History:

	02jun	ksl	Made separate routine, removing it from nebular concentrations
	06may	ksl	Swittch to plasma structure
*/


int
dlucy (xplasma)
     PlasmaPtr xplasma;
{
  int nelem, nion, niterate;
  double xne, xnew;
  double newden[NIONS];
  double t_r, nh;
  double t_e, www;
  int detailed_balance ();

  t_r = xplasma->t_r;
  t_e = xplasma->t_e;
  www = xplasma->w;


  /* Initally assume electron density from the LTE densities */

  xne = xplasma->ne;
  if (xne < DENSITY_MIN)
    {
      Log
	("Warning:nebular_concentrations: Very low ionization: ne initially %8.2e\n",
	 xne);
      xne = DENSITY_MIN;
    }

  nh = xplasma->rho * rho2nh;	//LTE -- Not clear needed at this level

  /* Carry out the detailed balance calculations first */

  for (nelem = 0; nelem < 2; nelem++)
    {
      detailed_balance (xplasma, nelem, newden);
      wind_update_after_detailed_balance (xplasma, nelem, newden);

    }


  /* Begin iteration loop to find ne from remaining elements */
  niterate = 0;
  while (niterate < MAXITERATIONS)
    {
      for (nelem = 2; nelem < nelements; nelem++)
	{
	  saha (xne, nh, t_r, nelem, xplasma->density);
	  lucy_mazzali1 (nh, t_r, t_e, www, nelem, xplasma->ne,
			 xplasma->density, xne, newden);
	}
      for (nion = 0; nion < nions; nion++)
	{
	  /*Set some floor so future divisions are sensible */
	  if (newden[nion] < DENSITY_MIN)
	    newden[nion] = DENSITY_MIN;
	}
      xnew = get_ne (newden);
      if (xnew < DENSITY_MIN)
	xnew = DENSITY_MIN;

      /* Check to see whether the search for xne has converged and if so exit loop */
      if (fabs ((xne - xnew) / (xnew)) < FRACTIONAL_ERROR || xnew < 1.e-6)
	break;

      /* else star another iteration of the main loop */
      xne = (xnew + xne) / 2.;	/* Make a new estimate of xne */
      niterate++;
    }
  /* End of main iteration loop */

  if (niterate == MAXITERATIONS)
    {
      Error
	("nebular_concentrations: failed to converge:nh %8.2e www %8.2e t_e %8.2e  t_r %8.2e \n",
	 nh, www, t_e, t_r);
      return (-1);
    }

/* Finally transfer the calculated densities to the real density array */

  xplasma->ne = xnew;
  for (nion = 0; nion < nions; nion++)
    {
      xplasma->density[nion] = newden[nion];
    }
  return (0);
}
