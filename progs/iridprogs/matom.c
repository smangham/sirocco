

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "my_linalg.h"

/*****************************************************************************/

/* As for similar routines in recomb.c, in order to use the integrator the 
   following external structures are used (SS)*/

/* This relates to the alpha_sp routines at the end of this file */

struct topbase_phot *cont_ext_ptr;	//continuum pointer passed externally
double temp_ext;		//temperature passed externally
int temp_choice;		//choice of type of calcualation for alpha_sp

// Note this should be moved to be closer to the routines


/************************************************************
                                    Imperial College London
Synopsis:
	Matom is the core of the implementation of Leon Lucy's Macro Atom
	approach to dealing with radiation-matter interactions. It is called whenever
	a photon packet activates a macro atom. As input it takes the process by
	which activation occurred (identified by the label "nres") which allow it 
	to deduce the level that has been excited. It then calculates all the 
	Macro Atom jumping/deactivaion probabilities following Lucy and so deduces
	the process by which deactivation occurs. At output "nres" identifies this 
	process and the packet information has been updates. 

Arguments:

       WindPtr w                   the ptr to the structure defining the wind
       PhotPtr p                   the packet at the point of activation
       int nres                    the process which activates the Macro Atom
  
Returns:
       int nres                    the process by which deactivation occurs
       PhotPtr p                   the packet following deactivation
       int escape                  to tell us whether the de-activation of the
                                   macro atom is via an r-packet (escape = 1) or
                                   a k-packet (escape = 0)

Description:


Notes: 

SS -- Currently this routine does not include collisional processes. These should be added.

ksl-- There is a lot of arithmetic required to keep track of indices.  I would
be inclined to figure out a way to avoid this.  One possibility would be to record the
go_to level in an array when one is calculating the emission probabilities. This would
make it easier to add new processes, I suspect. 

It would be possible to convert the for loop to a while statement.  This would avoid the
break statement in the middle.  

ksl -- 06aug -- At present, we do not have a way to dynamically allocated the structures
associated with matom, namely macromain, considerations one might do are documeented
in gridwind, where macromain is allocated. 


History:
          Feb 04  SS   Coding began.
          Mar 04  SS   Minor changes made (based on comments from ksl)
	04apr	ksl	A. Eliminated some unnecessary if statements.  These are indicated
			by the words extra.  Stuart should remove the lines assuming he
			agrees.
			B. Changed the logic of the major do loop so that the break is moved
			higher in the do loop.  Stuart, the intent is to make the
			code more readable, and to reduce the level of nesting. 
			C. Corrected an errror that on upward bb transitions seemed to leave
			the macro-atom in the same state as previously.
			D. Shamelessly modified some of the comments to make it more
			straightforward for me to understand.
			E. m = m + 1 --> m++ in several places
			F. Modified selections to allow simple lines to be read in but
			not to affect matom.  nlines_macro is the number on macro-lines
			read in and these are the first elements in the line structure
			G. Modified the way matom determines that a transition is a photo
			ionization transition bo be consitent with change made in
			resonate.
          Apr 04   SS   Modifications to include collisional jumps/deactivation included.
          May 04   SS   Bug corrected for identification of line (nres is the place of the
                        line in the ORDERED list, not the input list).
          Jun 04   SS   Modified to carry the "escape" variable to identify de-activation
                        via r-packets and k-packets - avoids the need to call kpkt from
                        within this routine
          Jun 04   SS   Adding collisional ionization as an upwards jump possibility. No collisional
                        recombination for the moment.
          Jun 04   SS   Now putting in collisional recombination too.
          July04   SS   Modifying so that this routine does not call alpha_sp but rather uses 
                        pre-computed (stored values) for the spontaneous recombination coefficient.
                        This is an effort to speed up the code.
	06may	ksl	57+ -- Adapted for use with plsama structure.  Changed call to
			eliminate passing entire w array
	06jun	ksl	57g -- Split macro variables into a separate structure. The structue
			which is createdin gridwind, is only created if there are macroatoms.
	06jul	ksl	57h -- In the process of speeding up the program in the simple 
			non-macro atom case, I tried to make sure that matom and the 
			derivative routiens were not called at all.  Note that at present
			the macroatom case is quite slow, due to what is happening in
			matom.  I am suspicious that it could be speeded up a lot.

************************************************************/

int
matom (p, nres, escape)
     PhotPtr p;
     int *nres;
     int *escape;
{
  struct lines *line_ptr;
  struct topbase_phot *cont_ptr;
  int uplvl, uplvl_old;
  double jprbs[2 * (NBBJUMPS + NBFJUMPS)];
  double eprbs[NBBJUMPS + NBFJUMPS];
  double pjnorm, penorm;
  double threshold, run_tot;
  double sp_rec_rate;
  int n, m;
  int njumps;
  int nbbd, nbbu, nbfd, nbfu;
  double t_e, ne;
  double bb_cont, choice, bf_cont;
  WindPtr one;
  double rad_rate, coll_rate;
  PlasmaPtr xplasma;
  MacroPtr mplasma;

  one = &wmain[p->grid];	//This is to identify the grid cell in which we are
  xplasma = &plasmamain[one->nplasma];
  check_plasma (xplasma, "matom");

  mplasma = &macromain[xplasma->nplasma];


  t_e = xplasma->t_e;		//electron temperature 
  ne = xplasma->ne;		//electron number density


  /* The first step is to identify the configuration that has been excited. */

/* ksl -- This approach will not work in a situation where there are extra non-macro lines */
/* ksl -- bug, I think.  I have run into some situations where the routine goes more than
max jumpts.  The one time I tracked this it occurred in a situation in which jumps were
not jumping to a new upper level.  I would have guested this was impossible */
/* SS Apr 04 -- both these issues should now be sorted out: this routine should never be called
   for non-macro lines. The bug is fixed. Jumps should always lead to a new level. */


  if (*nres < NLINES)		//this means that it was a line excitation CHECK
    {
      uplvl = lin_ptr[*nres]->nconfigu;
    }
  else if (*nres > NLINES)	//this means it was photoionisation
    {
      uplvl = phot_top[*nres - NLINES - 1].uplev;
    }
  else
    {
      Error ("matom: upper level not identified. nres = %d\n", *nres);
      exit (0);
    }

  /* Now follows the main loop to govern the macro atom jumps. Keeps jumping until
     an emission occurs (at which point it breaks). */


  for (njumps = 0; njumps < MAXJUMPS; njumps++)
    {
      /*  The excited configuration is now known. Now compute all the probabilities of deactivation
         /jumping from this configuration. Then choose one. */


      nbbd = config[uplvl].n_bbd_jump;	//store these for easy access -- number of bb downward jumps
      nbbu = config[uplvl].n_bbu_jump;	// number of bb upward jump from this configuration
      nbfd = config[uplvl].n_bfd_jump;	// number of bf downared jumps from this transition
      nbfu = config[uplvl].n_bfu_jump;	// number of bf upward jumps from this transiion


      // Start by setting everything to 0

      m = 0;			//m counts the total number of possible ways to leave the level
      pjnorm = 0.0;		//stores the total jump probability
      penorm = 0.0;		//stores the total emission probability

      for (n = 0; n < nbbd + nbfd; n++)
	{
	  eprbs[n] = 0;		//stores the individual emission probabilities SS
	  jprbs[n] = 0;		//stores the individual jump probabilities SS
	}
      for (n = nbbd + nbfd; n < nbbd + nbfd + nbbu + nbfu; n++)
	{
	  jprbs[n] = 0;		/*slots for upward jumps */
	}
      /* Finished zeroing. */


      /* bb */

      /* First downward jumps. (I.e. those that have emission probabilities. */

      /* For bound-bound decays the jump probability is A-coeff * escape-probability * energy */
      /* At present the escape probability is only approximated (p_escape). This should be improved. */
      /* The collisional contribution to both the jumping and deactivation probabilities are now added (SS, Apr04) */

      for (n = 0; n < nbbd; n++)
	{
	  line_ptr = &line[config[uplvl].bbd_jump[n]];
	  bb_cont =
	    (a21 (line_ptr) * p_escape (line_ptr, xplasma)) +
	    (q21 (line_ptr, t_e) * ne);

	  jprbs[m] = bb_cont * config[line_ptr->nconfigl].ex;	//energy of lower state

	  eprbs[m] = bb_cont * (config[uplvl].ex - config[line[config[uplvl].bbd_jump[n]].nconfigl].ex);	//energy difference

	  if (jprbs[m] < 0.)	//test (can be deleted eventually SS)
	    {
	      Error ("Negative probability (matom, 1). Abort.");
	      exit (0);
	    }
	  if (eprbs[m] < 0.)	//test (can be deleted eventually SS)
	    {
	      Error ("Negative probability (matom, 2). Abort.");
	      exit (0);
	    }

	  pjnorm += jprbs[m];
	  penorm += eprbs[m];
	  m++;
	}

      /* bf */
      for (n = 0; n < nbfd; n++)
	{
	  /* For bf recombination the jump probability is alpha_sp * nne * energy where
	     alpha_sp is the spontaneous recombination rate coeff. This is computed as
	     needed at the moment - perhaps this is too slow - should compute it and
	     store it? Stimulated recombination is included as negative photoionisation
	     (see below) (SS) */

	  /* Calls to alpha_sp are now replaced with use of recomb_sp which is the currently
	     stored value of alpha_sp. */

	  cont_ptr = &phot_top[config[uplvl].bfd_jump[n]];	//pointer to continuum
	  sp_rec_rate = mplasma->recomb_sp[uplvl][n];	//need this twice so store it
	  bf_cont = (sp_rec_rate + q_recomb (cont_ptr, t_e) * ne) * ne;

	  jprbs[m] = bf_cont * config[phot_top[config[uplvl].bfd_jump[n]].nlev].ex;	//energy of lower state
	  eprbs[m] = bf_cont * (config[uplvl].ex - config[phot_top[config[uplvl].bfd_jump[n]].nlev].ex);	//energy difference
	  if (jprbs[m] < 0.)	//test (can be deleted eventually SS)
	    {
	      Error ("Negative probability (matom, 3). Abort.");
	      exit (0);
	    }
	  if (eprbs[m] < 0.)	//test (can be deleted eventually SS)
	    {
	      Error ("Negative probability (matom, 4). Abort.");
	      exit (0);
	    }
	  pjnorm += jprbs[m];
	  penorm += eprbs[m];
	  m++;
	}



      /* Now upwards jumps. */

      /* bb */
      /* For bound-bound excitation the jump probability is B-coeff times Jbar with a correction 
         for stimulated emission. To avoid the need for recalculation all the time, the code will
         be designed to include the stimulated correction in Jbar - i.e. the stimulated correction
         factor will NOT be included here. (SS) */
      /* There is no emission probability for upwards transitions. */
      /* Collisional contribution to jumping probability added. (SS,Apr04) */

      for (n = 0; n < nbbu; n++)
	{
	  line_ptr = &line[config[uplvl].bbu_jump[n]];

	  jprbs[m] = ((b12 (line_ptr) * mplasma->jbar_old[uplvl][n]) + (q12 (line_ptr, t_e) * ne)) * config[uplvl].ex;	//energy of lower state
	  if (jprbs[m] < 0.)	//test (can be deleted eventually SS)
	    {
	      Error ("Negative probability (matom, 5). Abort.");
	      exit (0);
	    }
	  pjnorm += jprbs[m];
	  m++;
	}

      /* bf */
      for (n = 0; n < nbfu; n++)
	{
	  /* For bf ionization the jump probability is just gamma * energy
	     gamma is the photoionisation rate. Stimulated recombination also included. */
	  cont_ptr = &phot_top[config[uplvl].bfu_jump[n]];	//pointer to continuum

	  jprbs[m] = (mplasma->gamma_old[uplvl][n] - (mplasma->alpha_st_old[uplvl][n] * xplasma->ne * den_config (xplasma, cont_ptr->uplev) / den_config (xplasma, cont_ptr->nlev)) + (q_ioniz (cont_ptr, t_e) * ne)) * config[uplvl].ex;	//energy of lower state
	  if (jprbs[m] < 0.)	//test (can be deleted eventually SS)
	    {
	      Error ("Negative probability (matom, 6). Abort?\n");
	      jprbs[m] = 0.0;
	      //exit (0);
	    }
	  pjnorm += jprbs[m];
	  m++;
	}

      /* Probabilities of jumping (j) and emission (e) are now known. The normalisation factor pnorm has
         also been recorded. the integer m now gives the total number of possibilities too. 
         now select what happens next. Start by choosing the random threshold value at which the
         event will occur. */

      threshold = ((rand () + 0.5) / MAXRAND);


      if ((pjnorm / (pjnorm + penorm)) < threshold)
	break;			// An emission occurs and so we leave the for loop.

      uplvl_old = uplvl;

// Continue on if a jump has occured 

      /* Use a running total to decide where event occurs. */
      run_tot = 0;
// Note threshold could be 0 essentially, what keeps n from being -1 after this loop
      n = 0;
      threshold = ((rand () + 0.5) / MAXRAND);
      threshold = threshold * pjnorm;
      while (run_tot < threshold)
	{
	  run_tot += jprbs[n];
	  n++;
	}
      n = n - 1;

      /* n now identifies the jump that occurs - now set the new level. */
      if (n < nbbd)
	{			/* bb downwards jump */
	  uplvl = line[config[uplvl].bbd_jump[n]].nconfigl;
	}
      else if (n < (nbbd + nbfd))
	{			/* bf downwards jump */
	  uplvl = phot_top[config[uplvl].bfd_jump[n - nbbd]].nlev;
	}
      else if (n < (nbbd + nbfd + nbbu))
	{			/* bb upwards jump */
	  uplvl = line[config[uplvl].bbu_jump[n - nbbd - nbfd]].nconfigu;
	}
      else if (n < (nbbd + nbfd + nbbu + nbfu))
	{			/* bf upwards jump */
	  uplvl =
	    phot_top[config[uplvl].bfu_jump[n - nbbd - nbfd - nbbu]].uplev;
	}
      else
	{
	  Error ("Trying to jump but nowhere to go! Matom. Abort");
	  exit (0);
	}

/* ksl: Check added to verify that the level actually changed */
      if (uplvl_old == uplvl)
	{
	  Error ("matom: uplvl did not change with jump: %d %d\n", uplvl, n);

	}

    }

  /* When is gets here either the sum has reached maxjumps: didn't find an emission: this is
     an error and stops the run OR an emission mechanism has been chosen in which case all is well. SS */

  if (njumps == MAXJUMPS)
    {
      Error ("Matom: jumped %d times with no emission. Abort.\n", MAXJUMPS);
      exit (0);
    }

  /* If it gets here then an emission has occurred. SS */
  /* Use a running total to decide where event occurs. SS */

  run_tot = 0;
  n = 0;
  threshold = ((rand () + 0.5) / MAXRAND);
  threshold = threshold * penorm;	//normalise to total emission prob.
  while (run_tot < threshold)
    {
      run_tot += eprbs[n];
      n++;
    }
  n = n - 1;
  /* n now identifies the jump that occurs - now set nres for the return value. */
  if (n < nbbd)
    {				/* bb downwards jump */
      /* With collisions included we need to decide whether the deactivation is radiative (r-packet)
         or collisional (k-packet). Get a random number and then use the ratio of the collisional
         emission probability to the (already known) collisional+radiative probability to decide whether
         collisional or radiative deactivation occurs. */
      choice = ((rand () + 0.5) / MAXRAND);	// the random number

      line_ptr = &line[config[uplvl].bbd_jump[n]];	//pointer for the bb transition

      rad_rate = a21 (line_ptr) * p_escape (line_ptr, xplasma);
      coll_rate = q21 (line_ptr, t_e) * ne;

      if (choice > (coll_rate / (rad_rate + coll_rate)))
	{
	  *escape = 1;		//don't need to follow this routine with a k-packet
	  /* This is radiative deactivation. */
	  *nres = line[config[uplvl].bbd_jump[n]].where_in_list;
	  p->freq = line[config[uplvl].bbd_jump[n]].freq;
	  /* This is the comoving frequency - changed to rest frequency by doppler */
	}
      else
	{			/* This is collisional deactivation. In this case a k-packet is created and so it must be processed by
				   the appropriate routine. (SS, Apr04) */
	  /* This is now done by returning escape = 0 rather than calling kpkt -- SS June 04 */
	  *escape = 0;

	  //OLD   kpkt (w, p, nres);
	  /* That deals with the k-packet and follows the energy until it finds its way 
	     back to being an r-packet. kpkt will get the frequency of the new r-packet too
	     (SS) */
	}
    }
  else if (n < (nbbd + nbfd))
    {				/* bf downwards jump */
      /* With collisional recombination included we need to decide whether to deactivate
         radiatively or make a k-packet. */

      cont_ptr = &phot_top[config[uplvl].bfd_jump[n - nbbd]];

      rad_rate = mplasma->recomb_sp[uplvl][n - nbbd];	//again using recomb_sp rather than alpha_sp (SS July 04)
      coll_rate = q_recomb (cont_ptr, t_e) * ne;

      choice = ((rand () + 0.5) / MAXRAND);	// the random number

      if (choice > (coll_rate / (rad_rate + coll_rate)))
	{			//radiative deactivation
	  *escape = 1;
	  *nres = config[uplvl].bfd_jump[n - nbbd] + NLINES + 1;
	  /* continuua are indicated by nres > NLINES */
	  p->freq = phot_top[config[uplvl].bfd_jump[n - nbbd]].freq[0]
	    -
	    (log (1. - (rand () + 0.5) / MAXRAND) * xplasma->t_e / H_OVER_K);
	  /* Co-moving frequency - changed to rest frequency by doppler */
	  /*Currently this assumed hydrogenic shape cross-section - Improve */
	}
      else
	{			//collisional deactivation
	  *escape = 0;
	}

    }
  else
    {
      Error
	("Trying to emitt from Macro Atom but no available route (matom). Abort.");
      exit (0);
    }
  return (0);
}

/********************************************************************************/

/*
  b12  very similar to a21 - returns the b12 Einstein coefficient.
  History:
  2004feb       coded by S Sim
*/

// !! ksl OK, Probably should be moved to atomic.c for consistency, but that can be done
// later.  
#define B12_CONSTANT 5.01983e25

struct lines *b12_line_ptr;
double b12_a;

double
b12 (line_ptr)
     struct lines *line_ptr;
{
  double freq;

  if (b12_line_ptr != line_ptr)
    {
      freq = line_ptr->freq;
      b12_a = B12_CONSTANT * line_ptr->f / freq;
      b12_line_ptr = line_ptr;
    }

  return (b12_a);
}


/*****************************************************************************/

/*
alpha_sp - to govern the calculation of the spontaneous recombination rate.

The rate is given by (4 pi /c2) (gu/gl) (h2/2 pi m k T)^(3/2) 
times the integral of   a(nu) nu2 exp [(chi- h nu)/kT].

This sort of thing is already evaluated in functions like sigma_phot_topbase but
not quite as I want it here. So have written new code here rather than altering the 
other parts of the code. In the future it may be possible to eliminate the older
parts. (SS)

Modified May 04 (SS) to work for continua of simple elements too.

04jul30	ksl	Modified so that one does not need to have multiple versions
		of the code depending on how whether the integrand is 
		multiplied by 1, f/fthresh, or f/fthresh-1.  This was
		done eliminate alpha_sp_e as a separate set of routines
		and to assure that bf rates are positive 
			ichoice = 0   --> spontanous recombination
			ichoice = 1   --> energy weighted recombination 
			ichoice = 2   --> the difference between energy_weighted
					and spontaneous

	06may	ksl	57+ -- Modified to use plasma structure
*/
#define ALPHA_SP_CONSTANT 5.79618e-36

double
alpha_sp (cont_ptr, xplasma, ichoice)
     struct topbase_phot *cont_ptr;
     PlasmaPtr xplasma;
     int ichoice;
{
  double alpha_sp_value;
  double fthresh, flast;
  double qromb ();
  double alpha_sp_integrand ();

  temp_choice = ichoice;
  temp_ext = xplasma->t_e;	//external for use in alph_sp_integrand
  cont_ext_ptr = cont_ptr;	//"
  fthresh = cont_ptr->freq[0];	//first frequency in list
  flast = cont_ptr->freq[cont_ptr->np - 1];	//last frequency in list
  alpha_sp_value = qromb (alpha_sp_integrand, fthresh, flast, 1e-4);

  /* The lines above evaluate the integral in alpha_sp. Now we just want to multiply 
     through by the appropriate constant. */
  if (cont_ptr->macro_info == 1 && geo.macro_simple == 0)
    {
      alpha_sp_value = alpha_sp_value * config[cont_ptr->nlev].g
	/ config[cont_ptr->uplev].g * pow (xplasma->t_e, -1.5);
    }
  else				//case for simple element
    {
      alpha_sp_value = alpha_sp_value * config[cont_ptr->nlev].g / ion[cont_ptr->nion + 1].g * pow (xplasma->t_e, -1.5);	//g for next ion up used
    }

  alpha_sp_value = alpha_sp_value * ALPHA_SP_CONSTANT;

  return (alpha_sp_value);
}



/******************************************************************************/

/* alpha_sp_integrand. This returns the integrand for alpha_sp at a chosen
   frequency*/

double
alpha_sp_integrand (freq)
     double freq;		//frequency 
{
  double fthresh;
  double x;
  double integrand;
  double tt;

  fthresh = cont_ext_ptr->freq[0];
  tt = temp_ext;

  if (freq < fthresh)
    return (0.0);		// No recombination at frequencies lower than the threshold freq occur

  x = sigma_phot_topbase (cont_ext_ptr, freq);	//this is the cross-section
  integrand = x * freq * freq * exp (H_OVER_K * (fthresh - freq) / tt);

  if (temp_choice == 1)
    return (integrand * freq / fthresh);	//energy weighed case
  if (temp_choice == 2)
    return (integrand * (freq - fthresh) / fthresh);	// difference case
  return (integrand);		//spontanoues case
}

/*****************************************************************************/
/****************************************************************************/


/* kpkt
This deals with the elimination of k-packets. Whenever a k-packet is created it is
immediately destroyed insitu by this routine. At output "nres" identified the process
that destroys the k-packet and the packet information has been updated in the same
way as at the end of matom. */

/************************************************************
                                    Imperial College London
Synopsis:

Arguments:

       WindPtr w                   the ptr to the structure defining the wind
       PhotPtr p                   the packet at the point of activation
       int nres                    the process which activates the Macro Atom

Returns:
       int nres                    the process by which deactivation occurs
       PhotPtr p                   the packet following deactivation
       int escape                  indicated whether removal of k-packet was by creating
                                   an r-packet (escape = 1) or by exciting a macro atom
                                   (escape = 0)

Description:


Notes: 


History:
          Mar 04  SS   Coding began.
          Apr 04  SS   Various improvements including the addition of ff and collisions.
          May 04  SS   Minor changes made to collisional cooling rates (bug fixed)
                       and fb cooling for simple continua.
          May 04  SS   Modified to work for case with all "simple" ions.
          May 04  SS   Modified to use "scattering probability" formalism for 
                       simple ion cooling rates.
          Jun 04  SS   Modified to include the "escape" variable to identify excitation of
                       macro atoms. This removes the need to call matom from within this routine.
          Jun 04  SS   Modified to include collisional ionization as a cooling term.
          July04  SS   Modified to use recomb_sp(_e) rather than alpha_sp(_e) to speed up.
	06may	ksl	57+ -- Modified to use new plasma array.  Eliminated passing
			entire w array
          
************************************************************/

int
kpkt (p, nres, escape)
     PhotPtr p;
     int *nres;
     int *escape;
{

  int i;
  int ulvl;
  double cooling_bf[NTOP_PHOT];
  double cooling_bf_col[NTOP_PHOT];	//collisional cooling in bf transitions
  double cooling_bb[NLINES];
  struct topbase_phot *cont_ptr;
  struct lines *line_ptr;
  double cooling_normalisation;
  double destruction_choice;
  double electron_temperature;
  double cooling_bbtot, cooling_bftot, cooling_bf_coltot;
  double lower_density, upper_density;
  double cooling_ff;
  WindPtr one;
  PlasmaPtr xplasma;
  MacroPtr mplasma;

  double coll_rate, rad_rate;
  double q_ioniz ();

  /* For the moment only bf processes are included. Idea is to calculated the cooling
     terms for all the processes and then choose one at random to destroy the k-packet and
     turn it back into a photon bundle. 
     The loop here runs over all the bf processes, regardless of whether they are macro or simple
     ions. The reason that we can do this is that in the macro atom method stimulated recombination
     has already been considered before creating the k-packet (by the use of gamma-twiddle) in "scatter"
     and so we only have spontaneous recombination to worry about here for ALL cases. */

  /* Apr 04: now adding the possibility of removing k-packets by collisional excitation of a
     macro atom. To do this there's now a new array "cooling_bb" which has the collisional
     cooling rate for bb transitions. (SS) */

  /* Apr 04: now adding the possibility of removing k-packets by ff emission. (SS) */


  one = &wmain[p->grid];
  xplasma = &plasmamain[one->nplasma];
  check_plasma (xplasma, "kpkt");
  mplasma = &macromain[xplasma->nplasma];

  electron_temperature = xplasma->t_e;
  cooling_normalisation = 0.0;
  cooling_bftot = 0.0;
  cooling_bbtot = 0.0;
  cooling_ff = 0.0;
  cooling_bf_coltot = 0.0;

  for (i = 0; i < ntop_phot; i++)
    {
      cont_ptr = &phot_top[i];
      ulvl = cont_ptr->uplev;
      if (cont_ptr->macro_info == 1 && geo.macro_simple == 0)
	{
	  upper_density = den_config (xplasma, ulvl);
	  /* SS July 04 - for macro atoms the recombination coefficients are stored so use the
	     stored values rathert than recompue them. */
	  cooling_bf[i] = upper_density * H * cont_ptr->freq[0] * (mplasma->recomb_sp_e[ulvl][cont_ptr->down_index]);	// - modified the definition of _sp_e to be the difference
	  //           xplasma->recomb_sp[ulvl][cont_ptr->down_index]);
	}
      else
	{
	  upper_density = xplasma->density[cont_ptr->nion + 1];

	  /* SS Nov 04 - testing for speed - next line need to go back in. */

	  cooling_bf[i] = upper_density * H * cont_ptr->freq[0]
	    * (xplasma->recomb_simple[i]);
	  // cooling_bf[i] = 0;
	}


      /* Note that the electron density is not included here -- all cooling rates scale
         with the electron density so I've factored it out. */
      if (cooling_bf[i] < 0)
	{
	  Error ("kpkt: bf cooling rate negative. Density was %g\n",
		 upper_density);
	  Error ("alpha_sp(cont_ptr, xplasma,2) %g \n",
		 alpha_sp (cont_ptr, xplasma, 2));
	  Error ("i, ulvl, ntop_phot, nion %d %d %d %d\n", i, ulvl, ntop_phot,
		 cont_ptr->nion);
	  Error ("nlev, z, istate %d %d %d \n", cont_ptr->nlev, cont_ptr->z,
		 cont_ptr->istate);
	  Error ("freq[0] %g\n", cont_ptr->freq[0]);
	  //      exit(0);
	  cooling_bf[i] = 0.0;
	}
      else
	{
	  cooling_bftot += cooling_bf[i];
	}
      cooling_normalisation += cooling_bf[i];
      if (cont_ptr->macro_info == 1 && geo.macro_simple == 0)
	{
	  /* Include collisional ionization as a cooling term in macro atoms. Don't include
	     for simple ions for now.  SS */

	  lower_density = den_config (xplasma, cont_ptr->nlev);
	  cooling_bf_col[i] = lower_density * H * cont_ptr->freq[0]
	    * q_ioniz (cont_ptr, electron_temperature);

	  cooling_bf_coltot += cooling_bf_col[i];

	  cooling_normalisation += cooling_bf_col[i];

	}


    }

  for (i = 0; i < nlines; i++)
    {
      line_ptr = &line[i];
      if (line_ptr->macro_info == 1 && geo.macro_simple == 0)
	{			//It's a macro atom line and so the density of the upper level is stored
	  cooling_bb[i] =
	    den_config (xplasma, line_ptr->nconfigl) * q12 (line_ptr,
							    electron_temperature)
	    * line_ptr->freq * H;
	  /* Note that the electron density is not included here -- all cooling rates scale
	     with the electron density so I've factored it out. */
	}
      else
	{			//It's a simple line. Get the upper level density using two_level_atom

	  two_level_atom (line_ptr, xplasma, &lower_density, &upper_density);
	  coll_rate = q21 (line_ptr, electron_temperature)
	    * (1. - exp (-H_OVER_K * line_ptr->freq / electron_temperature));

	  cooling_bb[i] =
	    (lower_density * line_ptr->gu / line_ptr->gl -
	     upper_density) * coll_rate / (exp (H_OVER_K * line_ptr->freq /
						electron_temperature) -
					   1.) * line_ptr->freq * H;

	  rad_rate = a21 (line_ptr) * p_escape (line_ptr, xplasma);

	  /* Now multiply by the scattering probability - i.e. we are only going to consider bb cooling when
	     the photon actually escapes - we don't to waste time by exciting a two-level macro atom only so that
	     it makes another k-packet for us! (SS May 04) */


	  cooling_bb[i] *= rad_rate / (rad_rate + (coll_rate * xplasma->ne));
	}
      if (cooling_bb[i] < 0)
	{
	  Error ("kpkt: bb cooling rate negative. Abort");
	  exit (0);
	}
      else
	{
	  cooling_bbtot += cooling_bb[i];
	}
      cooling_normalisation += cooling_bb[i];
    }


/* 57+ -- This might be modified later since we "know" that xplasma cannot be for a grid with zero
volume.  Recall however that vol is part of the windPtr */
  if (one->vol > 0)
    {
      cooling_ff =
	total_free (one, xplasma->t_e, 0.0,
		    INFINITY) / one->vol / xplasma->ne;
    }
  else
    {
      /* SS June 04 - This should never happen, but sometimes it does. I think it is because of 
         photons leaking from one cell to another due to the push-through-distance. It is sufficiently
         rare (~1 photon in a complete run of the code) that I'm not worrying about it for now but it does
         indicate a real problem somewhere. */

      cooling_ff = 0.0;
      Error ("kpkt: A scattering event in cell %d with vol = 0???\n",
	     one->nwind);
//Diagnostic      return(-1);  //57g -- Cannot diagnose with an exit
      exit (0);
    }


  if (cooling_ff < 0)
    {
      Error ("kpkt: ff cooling rate negative. Abort.");
      exit (0);
    }
  else
    {
      cooling_normalisation += cooling_ff;
    }

  /* The cooling rates for the recombination and collisional processes are now known. 
     Choose which process destroys the k-packet with a random number. */

  destruction_choice = ((rand () + 0.5) / MAXRAND) * cooling_normalisation;


  if (destruction_choice < cooling_bftot)
    {				//destruction by bf
      for (i = 0; i < ntop_phot; i++)
	{
	  if (destruction_choice < cooling_bf[i])
	    {
	      /* Having got here we know that desturction of the k-packet was via the process labelled
	         by i. Let's just check that i is a sensible number. */

	      if (i > ntop_phot - 1)
		{
		  Error
		    ("kpkt (matom.c): trying to destroy k-packet in unknown process. Abort.\n");
		  exit (0);
		}

	      /* If it gets here, all seems fine. Now set nres for the destruction process. */

	      *nres = i + NLINES + 1;
	      *escape = 1;	//record that an r-packet is made - no need to excite a macro atom again


	      /* Now (as in matom) choose a frequency for the new packet. */

	      p->freq = phot_top[i].freq[0]
		-
		(log (1. - (rand () + 0.5) / MAXRAND) * xplasma->t_e /
		 H_OVER_K);
	      /* Co-moving frequency - changed to rest frequency by doppler */
	      /*Currently this assumed hydrogenic shape cross-section - Improve */

	      /* k-packet is now eliminated. All done. */

	      return (0);
	    }
	  else
	    {
	      destruction_choice = destruction_choice - cooling_bf[i];
	    }
	}
    }
  else if (destruction_choice < (cooling_bftot + cooling_bbtot))
    {				//this means that a collisional destruction has occured - this results in 
      //a macro atom being excited. Choose which macro atom and level to excite
      destruction_choice = destruction_choice - cooling_bftot;
      for (i = 0; i < nlines; i++)
	{
	  if (destruction_choice < cooling_bb[i])
	    {			//This is the bb collision which removes the k-packet
	      *nres = line[i].where_in_list;	//label for bb process 
	      if (line[i].macro_info == 1 && geo.macro_simple == 0)	//line is for a macro atom
		{

		  /* SS June 04 - in this case the escape = 0 flag is now returned to tell macro_gov that
		     a macro atom should be excited, rather than making a call to matom here. */

		  *escape = 0;

		  // OLD   matom (w, p, nres);
		}
	      else		//line is not for a macro atom - use simple method
		{
		  /* Since the cooling rate now accounts for the scattering fraction we know that if we
		     get here we want a line emission, not just an excited macro atom. (SS May 04) */
		  *escape = 1;	//No need for re-exciting a macro atom.
		  p->freq = line[i].freq;
		}
	      /* When it gets here the packet is back to an
	         r-packet and the emission mechanism is identified by nres
	         i.e. that's it finished. (SS, Apr 04). */

	      return (0);
	    }
	  else
	    {
	      destruction_choice = destruction_choice - cooling_bb[i];
	    }
	}
    }
  else if (destruction_choice < (cooling_bftot + cooling_bbtot + cooling_ff))
    {				//this is a ff destruction

      /* For the moment I've just set this to work in the frequency range 400000AA up to
         twice the He II edge. The range should probably be fixed to something that
         depends on the current interval freqmin -> freqmax, but at present
         freqmin and freqmax are not known here - they's need to be passed from "main"
         somehow. Perhaps they could be made global? (SS, May04). */

      *escape = 1;		//we are making an r-packet not exciting a macro atom

      *nres = -2;

      p->freq = one_ff (one, 7.5e12, 2.626e16);	//get frequency of resulting energy packet

      return (0);		//from one_ff
    }
  else
    {
      /* We want destruction by collisional ionization in a macro atom. */
      destruction_choice =
	destruction_choice - cooling_bftot - cooling_bbtot - cooling_ff;

      for (i = 0; i < ntop_phot; i++)
	{
	  if (destruction_choice < cooling_bf_col[i])
	    {
	      /* Having got here we know that destruction of the k-packet was via the process labelled
	         by i. Let's just check that i is a sensible number. */

	      if (i > ntop_phot - 1)
		{
		  Error
		    ("kpkt (matom.c): trying to destroy k-packet in unknown process. Abort.\n");
		  exit (0);
		}

	      /* If it gets here, all seems fine. Now set nres for the destruction process. */

	      *nres = i + NLINES + 1;
	      *escape = 0;	//collisional ionization makes an excited macro atom


	      /* k-packet is now eliminated. All done. */

	      return (0);
	    }
	  else
	    {
	      destruction_choice = destruction_choice - cooling_bf_col[i];
	    }
	}
    }



  Error ("matom.c: Failed to select a destruction process in kpkt. Abort.\n");
  Error
    ("matom.c: cooling_bftot %g, cooling_bbtot %g, cooling_ff %g, cooling_bf_coltot %g\n",
     cooling_bftot, cooling_bbtot, cooling_ff, cooling_bf_coltot);

  exit (0);

  return (0);
}



/************************************************************
                                    Imperial College London
Synopsis:
       fake_matom_bb is the macro atom routine that deals with line events involving
       simple ions (i.e. ions for which a full macro atom treatment is not employed.
       When this routine is called a simple line has absorbed a packet. This routine 
       creates a fake two-level macro atom and determines whether the packet energy
       is simply re-emitted in the line or is thermalised. If it is thermalised it
       turns into a k-packet and the appropriate routine is called. 


Arguments:

       WindPtr w                   the ptr to the structure defining the wind
       PhotPtr p                   the packet at the point of activation
       int nres                    the process which activates the Macro Atom

Returns:
       int nres                    the process by which deactivation occurs
       PhotPtr p                   the packet following deactivation
       int escape                  identifies whether the macro atom deactivated via an
                                   r-packet (escape = 1) or a k-packet (escape = 0). If
                                   a k-packet then the call to this routine should be
                                   followed by a call to kpkt.

Description:


Notes: 

History:
          Apr 04  SS   Coding began.
          Jun 04  SS   Modified to return escape = 1 for r-packet and 2 for k-packet
                       to avoid call to k-packet within this routine.
	06may	ksl	57+ -- Eliminationg passing entire w structure

************************************************************/

int
fake_matom_bb (p, nres, escape)
     PhotPtr p;
     int *nres;
     int *escape;
{
  double kprb, rprb;
  WindPtr one;
  PlasmaPtr xplasma;
  struct lines *line_ptr;
  double electron_temperature;
  double normalisation;
  double choice;

  one = &wmain[p->grid];	//record where we are
  xplasma = &plasmamain[one->nplasma];

  line_ptr = lin_ptr[*nres];	//record the line pointer
  electron_temperature = xplasma->t_e;

  /* Upon calling we know that the upper level of our fake two level macro
     atom is excited. Since it's only two-levels there are no jumping probabilities
     and the only decision that needs to be made it does it de-excite by an r-packet
     in the line or a collisionally produced k-packet. To decide this we compute the
     two emission probabilities. rprb is the probability of the r-packet and kprb
     for the k-packet. */

  /* Since both k- and r-packet emission have the same delta energy associated with them then
     we can ignore the energy factor that is otherwise needed in the de-activation probabilities.
     The upper level population is also factored out. (SS) */

  rprb = a21 (line_ptr) * p_escape (line_ptr, xplasma);
  // Einstein-A * escape probability

  kprb =
    q21 (line_ptr,
	 electron_temperature) * xplasma->ne * (1. -
						exp (-H_OVER_K *
						     line_ptr->freq /
						     electron_temperature));
  // Collision de-excitation rate coeff. * electron density 
  /* The extra factor of (1 - exp (-h nu / k T)) is explained in KSL's notes on Python. It appears because we are
     using the "scattering fraction" formalism for simple ions. (SS May 04). */

  normalisation = kprb + rprb;

  kprb = kprb / normalisation;

  rprb = rprb / normalisation;

  /* Now just use a random number to decide what happens. */

  choice = ((rand () + 0.5) / MAXRAND);

  /* If "choice" is less than rprb then we have chosen a radiative decay - for this fake macro atom there 
     is only one line so there's nothing to do - the energy is re-radiated in the line and that's it. We
     don't need to change nres. Otherwise we've chosen a collisional destruction and so we need to get
     kpkt to deal with the resulting k-packet and give us a new value of nres following thermalisation and
     re-emission of the energy. */

  if (choice < rprb)
    {
      *escape = 1;		//it's an r-packet de-activation
      p->freq = line_ptr->freq;	// As in matom, set this to the comoving frequency
    }
  else
    {
      /* Rather than calling kpkt here we are now going to return the instruction to macro_gov to
         make a k-packet. */

      *escape = 0;

    }

  /* That's it - we've dealt with the simple line. */

  return (0);

}


/************************************************************
                                    Imperial College London
Synopsis:
       fake_matom_bf is the macro atom routine that deals with photoionisation 
       events involving simple ions (i.e. ions for which a full macro atom treatment 
       is not employed).
       When this routine is called a photoionisation has absorbed a packet. 
       The idea of this routine is to deal with the subsequenct
       by creating a fake two-level atom. 
       However, in the absense of collisional recombination (or something
       similar) there's only one thing that can happen - radiative 
       recombination. Therefore there's no need to do anything here unless
       collisional recombination is to be introduced at some point in the
       future.
       All this routine does for now is choose a new frequency for the 
       emitted r-packet.


Arguments:

       WindPtr w                   the ptr to the structure defining the wind
       PhotPtr p                   the packet at the point of activation
       int nres                    the process which activates the Macro Atom

Returns:
       int nres                    the process by which deactivation occurs
       PhotPtr p                   the packet following deactivation
       int escape                  in principle this tells us whether de-activation is
                                   via an r-packet or a k-packet. For this routine at the 
                                   moment only r-packets are possible so it always returns
                                   escape = 1

Description:


Notes:    

History:
          Apr 04  SS   Coding began.
          Jun 04  SS   Modified to include "escape" being set to 1
	06may	ksl	57+ -- Modified for new structure.  Have not
			fixed the call to fake_atom

************************************************************/

int
fake_matom_bf (p, nres, escape)
     PhotPtr p;
     int *nres;
     int *escape;
{
  WindPtr one;
  PlasmaPtr xplasma;

  one = &wmain[p->grid];	//record where we are
  xplasma = &plasmamain[one->nplasma];

  /* All that can happen is radiative recombination (no collisional recombination
     for now. Just choose frequency for re-emission). */

  *escape = 1;			//always an r-packet here

  p->freq = phot_top[*nres - NLINES - 1].freq[0]
    - (log (1. - (rand () + 0.5) / MAXRAND) * xplasma->t_e / H_OVER_K);

  /*Currently this assumes hydrogenic shape cross-section - Improve */

  return (0);

}



  /************************************************************
                                    Imperial College London
Synopsis:

       macro_pops uses the Monte Carlo estimators to compute a set
       of level populations for levels of macro atoms.


Arguments: one -> pointer to the wind 
           xne -> current value for electron density in this shell

Returns:   Should compute the fractional level populations for 
          macro atoms and store them in "levden" array. The ion fractions
          are also computed and stored in w[n].density[nion]


Description:


Notes:   Although it seems to be working, I've not tested this routine
         carefully yet - SS 27 Apr 04.

         This routine uses a matrix inversion method to get the level populations.
         For now the matrix solver used is the LU decomposition method provided
         by the Gnu Scientific Library (GSL, which is free). This requires the
         include files that I've added to the top of this file and access to the
         GSL "library" file (I've added the library into the Makefile too). I found
         GSL to be very easy to install but if there are problems in the future we
         may need to switch to another matrix solver. (SS, Apr 04)
         

History:
          Apr 04  SS   Coding began.
          May 04  SS   Checking/testing began. Corrected some problems with the
                       array indexing in the last part of the calculation.
          June04  SS   Adding collisional ionization/recombination rate.
	  Aug 05  SS   Modified to include fast collisional transition between
	               ground states and levels that should be metastable
                       - some treatment of metastables is needed for 
                       doing CNO macro atoms.
	06may	ksl	57+ -- Modified to reflect plasma structue

************************************************************/

int
macro_pops (xplasma, xne)
     PlasmaPtr xplasma;
     double xne;
{

  int index_element, index_ion, index_lvl;
  int n_macro_lvl;
  double rate;
  double rate_matrix[100][100];
  int matrix_to_conf[100];
  int conf_to_matrix[100];
  struct lines *line_ptr;
  struct topbase_phot *cont_ptr;
  int nn, mm;
  int index_bbu, index_bbd, index_bfu, index_bfd;
  int lower, upper;
  double this_ion_density;
  double inversion_test;
  double q_ioniz (), q_recomb ();
  int s;			//NEWKSL
  double *a_data, *b_data;
  gsl_permutation *p;		//NEWKSL
  gsl_matrix_view m;		//NEWKSL
  gsl_vector_view b;		//NEWKSL
  gsl_vector *populations;	//NEWKSL
  int index_fast_col;

  MacroPtr mplasma;
  mplasma = &macromain[xplasma->nplasma];

  /* Start with an outer loop over elements: there are no rates that couple
     levels of different elements so we can always separate them out. */

  for (index_element = 0; index_element < nelements; index_element++)
    {

      /* Zero all elements of the matrix before doing anything else. */

      for (nn = 0; nn < 100; nn++)
	{
	  for (mm = 0; mm < 100; mm++)
	    {
	      rate_matrix[mm][nn] = 0.0;
	    }
	}

      /* See if this element uses a macro atom treatment or is a simple element.
         For now I'm assuming that either all ions of a given element are
         treated using the macro atom method, or else none are (mixing and 
         matching is probably a bad idea because of the way in which bf 
         processes couple different ionisation stages). */

      /* The check is against the first ion of the element. */

      if (ion[ele[index_element].firstion].macro_info == 1
	  && geo.macro_simple == 0)
	{

	  /* Having established that the ion requires a macro atom treatment we
	     are going to construct a matrix of rates between the levels and 
	     invert that matrix to get the level populations. The first thing we need 
	     to do is work out how many levels we are dealing with in total. This is
	     easily done by summing up the number of levels of each ion. */

	  n_macro_lvl = 0;

	  for (index_ion = ele[index_element].firstion;
	       index_ion <
	       (ele[index_element].firstion + ele[index_element].nions);
	       index_ion++)
	    {
	      for (index_lvl = ion[index_ion].first_nlte_level;
		   index_lvl <
		   ion[index_ion].first_nlte_level + ion[index_ion].nlte;
		   index_lvl++)
		{
		  /* I want to be able to easily go from knowing the index of a level in the 
		     configurations structure to its position in the rates matrix. So I'm making
		     two arrays here that allow the mapping between these indices to be done easily.
		   */
		  conf_to_matrix[index_lvl] = n_macro_lvl;
		  matrix_to_conf[n_macro_lvl] = index_lvl;
		  n_macro_lvl++;
		}
	    }

	  //      Log("n_macro_lvl %d\n", n_macro_lvl);

	  /* We now know how many levels there are and therefore how big the matrix we
	     need to invert will be. */

	  /* Now we want to populate the matrix with all the rates between the levels. */

	  for (index_ion = ele[index_element].firstion;
	       index_ion <
	       (ele[index_element].firstion + ele[index_element].nions);
	       index_ion++)
	    {
	      index_lvl = ion[index_ion].first_nlte_level;

	      /* The next loop judges whether or not a level is to be fixed in population relative to ground
	         star. The input radiative lifetime is used to judge this at the moment. If the lifetime was set
	         to be long (essentially infite) then a very fast collisional transition is put in to dominate
	         all other rates into and out of this level. 

	         Whether this is really the best thing to do I don't know, but it's an improvement over ignoring 
	         this issue altogether! SS Aug 05 */

	      for (index_fast_col = index_lvl;
		   index_fast_col <
		   ion[index_ion].first_nlte_level + ion[index_ion].nlte - 1;
		   index_fast_col++)
		{
		  {
		    if (config[index_fast_col + 1].rad_rate > 1.e15)
		      {
			fast_line.gl = config[index_lvl].g;
			fast_line.gu = config[index_fast_col + 1].g;
			fast_line.freq =
			  (config[index_fast_col + 1].ex -
			   config[index_lvl].ex) / H;
			fast_line.f = 1e4;
			rate = q12 (&fast_line, xplasma->t_e) * xne;
			lower = conf_to_matrix[index_lvl];
			upper = conf_to_matrix[index_fast_col + 1];
			rate_matrix[lower][lower] += -1. * rate;
			rate_matrix[upper][lower] += rate;
			rate = q21 (&fast_line, xplasma->t_e) * xne;
			rate_matrix[upper][upper] += -1. * rate;
			rate_matrix[lower][upper] += rate;
			rate = q12 (&fast_line, (0.5 * xplasma->t_e));
		      }
		  }
		}

	      for (index_lvl = ion[index_ion].first_nlte_level;
		   index_lvl <
		   ion[index_ion].first_nlte_level + ion[index_ion].nlte;
		   index_lvl++)
		{
		  /* Now add contribution for all the bb and bf processes between the levels. This is
		     done by looping over the numbers "bbu, bbd, bfu, bfd" which tell us how many
		     processes there are. */


		  for (index_bbu = 0;
		       index_bbu < config[index_lvl].n_bbu_jump; index_bbu++)
		    {
		      /* These are bb upwards jumps. The rate in such a jump is given by
		         Jbar which has been computed as a Monte Carlo estimator. I'm also
		         including a collisional term (which depends on ne). */

		      line_ptr = &line[config[index_lvl].bbu_jump[index_bbu]];
		      rate =
			b12 (line_ptr) *
			mplasma->jbar_old[index_lvl][index_bbu];
		      rate += q12 (line_ptr, xplasma->t_e) * xne;

		      /* This is the rate out of the level in question. We need to add it
		         to the matrix in two places: firstly as a -ve contribution to the
		         diagonal and secondly as a +ve contribution for the off-diagonal
		         corresponding to the level populated by this process. */

		      /* Get the matix indices for the upper and lower states of the jump. */

		      lower = conf_to_matrix[index_lvl];
		      upper = conf_to_matrix[line_ptr->nconfigu];

		      rate_matrix[lower][lower] += -1. * rate;
		      rate_matrix[upper][lower] += rate;
		    }

		  for (index_bbd = 0;
		       index_bbd < config[index_lvl].n_bbd_jump; index_bbd++)
		    {
		      /* These are bb downwards jumps. The rate in such a jump is given by
		         the A-value. I'm also
		         including a collisional term (which depends on ne). */

		      line_ptr = &line[config[index_lvl].bbd_jump[index_bbd]];
		      rate = (a21 (line_ptr) * p_escape (line_ptr, xplasma));
		      rate += q21 (line_ptr, xplasma->t_e) * xne;

		      /* This is the rate out of the level in question. We need to add it
		         to the matrix in two places: firstly as a -ve contribution to the
		         diagonal and secondly as a +ve contribution for the off-diagonal
		         corresponding to the level populated by this process. */

		      /* Get the matix indices for the upper and lower states of the jump. */

		      upper = conf_to_matrix[index_lvl];
		      lower = conf_to_matrix[line_ptr->nconfigl];

		      rate_matrix[upper][upper] += -1. * rate;
		      rate_matrix[lower][upper] += rate;
		    }

		  for (index_bfu = 0;
		       index_bfu < config[index_lvl].n_bfu_jump; index_bfu++)
		    {
		      /* These are bf upwards jumps. The rate in such a jump is given by
		         gamma which has been computed as a Monte Carlo estimator. */

		      cont_ptr =
			&phot_top[config[index_lvl].bfu_jump[index_bfu]];
		      rate = mplasma->gamma_old[index_lvl][index_bfu];
		      rate += q_ioniz (cont_ptr, xplasma->t_e) * xne;

		      /* This is the rate out of the level in question. We need to add it
		         to the matrix in two places: firstly as a -ve contribution to the
		         diagonal and secondly as a +ve contribution for the off-diagonal
		         corresponding to the level populated by this process. */

		      /* Get the matix indices for the upper and lower states of the jump. */

		      lower = conf_to_matrix[index_lvl];
		      upper = conf_to_matrix[cont_ptr->uplev];

		      rate_matrix[lower][lower] += -1. * rate;
		      rate_matrix[upper][lower] += rate;

		      /* Now deal with the stimulated emission. */
		      /* Lower and upper are the same, but now it contributes in the
		         other direction. */

		      rate =
			mplasma->alpha_st_old[index_lvl][index_bfu] * xne;

		      rate_matrix[upper][upper] += -1. * rate;
		      rate_matrix[lower][upper] += rate;


		    }

		  for (index_bfd = 0;
		       index_bfd < config[index_lvl].n_bfd_jump; index_bfd++)
		    {
		      /* These are bf downwards jumps. The rate in such a jump is given by
		         the alpha value. */
		      cont_ptr =
			&phot_top[config[index_lvl].bfd_jump[index_bfd]];
		      /* Get new values of the recombination rates and store them. */
		      mplasma->recomb_sp[index_lvl][index_bfd] =
			alpha_sp (cont_ptr, xplasma, 0);
		      mplasma->recomb_sp_e[index_lvl][index_bfd] =
			alpha_sp (cont_ptr, xplasma, 2);
		      rate = mplasma->recomb_sp[index_lvl][index_bfd] * xne;
		      rate += q_recomb (cont_ptr, xplasma->t_e) * xne * xne;


		      /* This is the rate out of the level in question. We need to add it
		         to the matrix in two places: firstly as a -ve contribution to the
		         diagonal and secondly as a +ve contribution for the off-diagonal
		         corresponding to the level populated by this process. */

		      /* Get the matix indices for the upper and lower states of the jump. */

		      upper = conf_to_matrix[index_lvl];
		      lower = conf_to_matrix[cont_ptr->nlev];

		      rate_matrix[upper][upper] += -1. * rate;
		      rate_matrix[lower][upper] += rate;
		    }
		}
	    }

	  /* The rate matrix is now filled up. Since the problem is not closed as it stands, the next
	     thing is to replace one of the rows of the matrix (say the first row) with the contraint
	     that the sum of all the populations is 1.0 (this will let us get population fractions). */

	  for (index_lvl = 0; index_lvl < n_macro_lvl; index_lvl++)
	    {
	      rate_matrix[0][index_lvl] = 1.0;
	    }

	  /* Now we can just invert the matrix to get the fractional level populations. */


	  /********************************************************************************/
	  /* The block that follows (down to next line of ***s) is to do the
	     matix inversion. It uses something called LU decomposition - the
	     code for doing this is taken from the GSL manual with very few modifications. */


	  //  Error("n_macro_lvl %d \n", n_macro_lvl);

	  /* Replaced inline array allocaation with calloc, which will work with older version of c compilers */
//OLD     double a_data[n_macro_lvl * n_macro_lvl];

	  a_data =
	    (double *) calloc (sizeof (rate), n_macro_lvl * n_macro_lvl);

	  for (nn = 0; nn < n_macro_lvl; nn++)
	    {
	      for (mm = 0; mm < n_macro_lvl; mm++)
		{
		  a_data[nn * n_macro_lvl + mm] = rate_matrix[nn][mm];
		  //Error("a_data, rate_matrix %g %g %d \n", a_data[nn*n_macro_lvl + mm],rate_matrix[nn][mm], nn*n_macro_lvl + mm); 
		}
	      //                  Error ("%g %g %g %g %g \n", a_data[nn*n_macro_lvl],a_data[nn*n_macro_lvl+1],a_data[nn*n_macro_lvl+2],a_data[nn*n_macro_lvl+3],a_data[nn*n_macro_lvl+4]); 
	    }


	  /* Replaced inline array allocaation with calloc, which will work with older version of c compilers */

	  b_data = (double *) calloc (sizeof (rate), n_macro_lvl);
	  b_data[0] = 1.0;

	  m = gsl_matrix_view_array (a_data, n_macro_lvl, n_macro_lvl);	//KSLNEW

	  b = gsl_vector_view_array (b_data, n_macro_lvl);	//KSLNEW

	  populations = gsl_vector_alloc (n_macro_lvl);	//KSLNEW



	  p = gsl_permutation_alloc (n_macro_lvl);	//NEWKSL

	  gsl_linalg_LU_decomp (&m.matrix, p, &s);

	  gsl_linalg_LU_solve (&m.matrix, p, &b.vector, populations);

	  gsl_permutation_free (p);

	  free (a_data);	//NEWKSL
	  free (b_data);	//NEWKSL

	  /**********************************************************/

	  //for (nn=0; nn < n_macro_lvl; nn++)
	  //{
	  //  Error ("population %d %g \n", nn, gsl_vector_get(populations, nn));
	  // }

	  /* MC noise can cause population inversions (particularly amongst highly excited states)
	     which are never a good thing and most likely unphysical.
	     Therefor let's follow Leon's procedure (Lucy 2003) and remove inversions. */



	  for (index_ion = ele[index_element].firstion;
	       index_ion <
	       (ele[index_element].firstion + ele[index_element].nions);
	       index_ion++)
	    {
	      for (index_lvl = ion[index_ion].first_nlte_level;
		   index_lvl <
		   ion[index_ion].first_nlte_level + ion[index_ion].nlte;
		   index_lvl++)
		{		/* Start loop with lowest level of the ion. For each level in turn check to see if there's a population 
				   inversion i.e. is  upper_pop > lower_pop * g_upper / g_lower. If it is then replace upper_pop with
				   lower_pop * g_upper / g_lower. We loop over all levels higher than the currently chosen lower level. */
		  for (nn = index_lvl + 1;
		       nn <
		       (ion[index_ion].first_nlte_level +
			ion[index_ion].nlte); nn++)
		    {
		      inversion_test = gsl_vector_get (populations, conf_to_matrix[index_lvl]) * config[nn].g / config[index_lvl].g * 0.999999;	//include a correction factor 
		      if (gsl_vector_get (populations, conf_to_matrix[nn]) >
			  inversion_test)
			{
			  gsl_vector_set (populations, conf_to_matrix[nn],
					  inversion_test);
			}
		    }
		}
	    }

	  //      for (nn=0; nn < n_macro_lvl; nn++)
	  // {
	  //   Error ("population %d %g \n", nn, gsl_vector_get(populations, nn));
	  //}




	  /* The populations are now known. The populations need to be stored
	     firstly as ion populations and secondly as fractional
	     level populations within an ion. Get the ion
	     populations and write them to one->density[nion]. The level populations
	     are to be put in "levden". */

	  nn = 0;
	  mm = 0;
	  for (index_ion = ele[index_element].firstion;
	       index_ion <
	       (ele[index_element].firstion + ele[index_element].nions);
	       index_ion++)
	    {
	      this_ion_density = 0.0;
	      for (index_lvl = ion[index_ion].first_nlte_level;
		   index_lvl <
		   ion[index_ion].first_nlte_level + ion[index_ion].nlte;
		   index_lvl++)
		{
		  this_ion_density += gsl_vector_get (populations, nn);
		  nn++;
		}

	      //              Log("this ion density %g\n", this_ion_density);

	      /* Write the level populations to the 
	         levden array. These are fractional level populations within an ion. */

	      for (index_lvl = ion[index_ion].first_nlte_level;
		   index_lvl <
		   ion[index_ion].first_nlte_level + ion[index_ion].nlte;
		   index_lvl++)
		{

		  // The next couple of lines let you see how much the populations are changin on screen if there are in. 
		  //Error("Delta n / n %g\n",((one->levden[config[index_lvl].nden]) -
		  //(gsl_vector_get(populations, conf_to_matrix[mm])/this_ion_density))/
		  //    (one->levden[config[index_lvl].nden]));
		  //Error("Old Level %d %d %g %g\n", conf_to_matrix[index_lvl], index_lvl, 
		  //one->levden[config[index_lvl].nden], config[index_lvl].g);


		  xplasma->levden[config[index_lvl].nden] +=
		    gsl_vector_get (populations,
				    conf_to_matrix[index_lvl]) /
		    this_ion_density;

		  //Log("conf_to_matrix[mm] %d %d\n", conf_to_matrix[mm], mm, index_lvl);

		  xplasma->levden[config[index_lvl].nden] /= 2;


		  mm++;


		  // The next line was in so that I could watch the populations on screen.
		  //Log("New Level %d %d %g\n", conf_to_matrix[index_lvl], index_lvl, 
		  //    one->levden[config[index_lvl].nden]);
		}

	      xplasma->density[index_ion] +=
		this_ion_density * ele[index_element].abun * xplasma->rho *
		rho2nh;

	      xplasma->density[index_ion] /= 2;


	      // The next line was in so that I could watch the populations on screen.
	      //Log("Ion %d %g %g \n", index_ion, this_ion_density, one->density[index_ion]);
	      //Log(" \n");

	    }

	}
    }

  //  exit(0);
  return (0);
  /* All done. (SS, Apr 04) */
}

  /************************************************************
                                    Imperial College London
Synopsis:

       macro_gov is a new routine that will sit at a higher level in the code than either matom or kpkt
       and will govern the passage of packets between these routines. At the moment, since matom and kpkt 
       call each other it call all get rather confusing and loads of nested subroutine calls ensue. 
       macro_gov removes this by calling matom and kpkt which on return tell macro_gov to either return
       an r-packet to resonate or make another call to either kpkt or matom as appropriate.



Arguments:   
       WindPtr w                   the ptr to the structure defining the wind
       PhotPtr p                   the packet at the point of activation
       int nres                    the process which activates the Macro Atom
       int matom_or_kpkt           to tell us if we should initially excite a
                                   Macro Atom or create a k-packet. For Macro Atom
                                   it is set to 1 for a k-packet it is set to 2.
       int which_out               set to 1 if return is via macro atom and 2 if via kpkt

           

Returns:   Will return an r-packet after (possibly) several calls to matom and kpkt

       int nres                    the process by which re-emission occurs
       PhotPtr p                   the packet following re-emission



Description:


Notes:   I've written this to be as general as possible so that if we want to improve the treatment of simple
         ions it should not need to be re-written.

         During the spectrum calculation the emission of r-packets within the spectral region of interest
         is done using emissivities obtained during the ionization cycles. Therefore whenever an r-packet
         is converted into a k-packet or an excited macro atom that ends the need to follow the packet any
         further. To deal with this, this routine sets the weights of such packets to zero. Upon return to 
         trans_phot these packets will then be thrown away.

         
History:
          SS   June 04 - coding began
          SS   June 04 - modified to kill photons that interact with macro atoms or
                         become k-packets during the specrum calculation part of the code.
          SS   June 04 - returns 1 if last interaction was with full macro atom. 0 otherwise.
          SS   Sep  04 - added return of "which_out"
************************************************************/

int
macro_gov (p, nres, matom_or_kpkt, which_out)
     PhotPtr p;
     int *nres;
     int matom_or_kpkt;
     int *which_out;

{
  int escape;			//this tells us when the r-packet is escaping 

  escape = 0;			//start with it not being ready to escape as an r-packet

  while (escape == 0)
    {
      if (matom_or_kpkt == 1)	//excite a macro atom - dependig on simple/macro choice call different routines
	{

	  if (*nres > (-1) && *nres < NLINES && geo.macro_simple == 0
	      && lin_ptr[*nres]->macro_info == 1)
	    /* This is a bb line of a macro atoms and we want the full macro atom treatment. */
	    {

	      if (geo.matom_radiation == 1)
		{
		  /* During the spectrum cycles we want to throw these photons away. */
		  p->w = 0.0;
		  escape = 1;	//This doesn't matter but it breaks us out of this loop
		}
	      else
		{

		  matom (p, nres, &escape);

		  if (escape == 1)
		    {
		      /* It is going to escape as a r-packet that was created by de-activation of a macro atom. 
		         Therefore, if the frequency is suitable it should be recorded as a macro atom emission event for use in the 
		         computation of the k-packet emissivity needed for the final spectrum calculation. */
		      //   if (p->freq > em_rnge.fmin && p->freq < em_rnge.fmax)
		      //        {
		      /*We do want to include it. Emissivities are stored by level so we need to identify the upper level
		         that produced this emission. */
		      //          if (*nres > (-1) && *nres < NLINES)
		      //            {
		      /* It was bb. */
		      //              upper = lin_ptr[*nres]->nconfigu;
		      //            }
		      //  else if (*nres > NLINES)
		      //    {
		      /* It was bf. */
		      //      upper = phot_top[*nres - NLINES - 1].uplev;
		      //    }
		      //  else
		      //    {
		      //      Error
		      //        ("macro_gov: cannot identify deactivation nres. Abort.\n");
		      //      exit (0);
		      //    }
		      /* Now record emissivity. */
		      //w[p->grid].matom_emiss[upper] += p->w;
		      //w[p->grid].matom_nemiss[upper] += 1;
		      //}
		      *which_out = 1;
		      return (1);
		    }
		}
	    }
	  else if (*nres > (-1) && *nres < NLINES
		   && (geo.macro_simple == 1
		       || lin_ptr[*nres]->macro_info == 0))
	    /* This is a bb line for which we don't want a full macro atom treatment. */
	    {
	      fake_matom_bb (p, nres, &escape);
	    }
	  else if (*nres > NLINES
		   && phot_top[*nres - NLINES - 1].macro_info == 1
		   && geo.macro_simple == 0)
	    /* This is a bf continuum of a macro atom and we want the full treatment. */
	    {

	      if (geo.matom_radiation == 1)
		{
		  /* During the spectrum cycles we want to throw these photons away. */
		  p->w = 0.0;
		  escape = 1;	//This doesn't matter but it breaks us out of this loop
		}
	      else
		{
		  matom (p, nres, &escape);

		  if (escape == 1)
		    {
		      /* It is going to escape as a r-packet that was created by de-activation of a macro atom. 
		         Therefore, if the frequency is suitable it should be recorded as a macro atom emission event for use in the 
		         computation of the k-packet emissivity needed for the final spectrum calculation. */
		      //if (p->freq > em_rnge.fmin && p->freq < em_rnge.fmax)
		      //{
		      /*We do want to include it. Emissivities are stored by level so we need to identify the upper level
		         that produced this emission. */
		      //  if (*nres > (-1) && *nres < NLINES)
		      //    {
		      /* It was bb. */
		      //      upper = lin_ptr[*nres]->nconfigu;
		      //    }
		      //  else if (*nres > NLINES)
		      //    {
		      /* It was bf. */
		      //      upper = phot_top[*nres - NLINES - 1].uplev;
		      //    }
		      //  else
		      //    {
		      //      Error
		      //        ("macro_gov: cannot identify deactivation nres. Abort.\n");
		      //      exit (0);
		      //    }
		      /* Now record emissivity. */
		      //w[p->grid].matom_emiss[upper] += p->w;
		      //w[p->grid].matom_nemiss[upper] += 1;
		      //}
		      *which_out = 1;
		      return (1);
		    }
		}
	    }
	  else if (*nres > NLINES
		   && (phot_top[*nres - NLINES - 1].macro_info == 0
		       || geo.macro_simple == 1))
	    /* This is a bf continuum but we don't want the full macro atom treatment. */
	    {
	      fake_matom_bf (p, nres, &escape);
	    }

	  matom_or_kpkt = 2;	//if it did not escape then it must have had a
	  //de-activation by collision processes -> need a k-packet 
	  //next
	}
      else if (matom_or_kpkt == 2)	//deal with a k-packet
	{
	  if (geo.matom_radiation == 1)
	    {
	      /* During the spectrum cycles we want to throw these photons away. */
	      p->w = 0.0;
	      escape = 1;	//This doesn't matter but it breaks us out of the loop.
	    }
	  else
	    {
	      kpkt (p, nres, &escape);

	      //if (escape == 1)
	      //{
	      /* It is going to escape as a r-packet that was last a k-packet. Therefore, if the frequency
	         is suitable it should be recorded as a k-packet emission event for use in the computation of
	         the k-packet emissivity needed for the final spectrum calculation. */
	      //if (p->freq > em_rnge.fmin && p->freq < em_rnge.fmax)
	      // {
	      //  w[p->grid].kpkt_emiss += p->w;
	      // }
	      //}
	    }

	  matom_or_kpkt = 1;	//if it did not escape then the k-packet must have been 
	  //destroyed by collisionally exciting a macro atom - 
	  //excite a macro atom next
	}
      else
	{
	  Error ("macro_gov: Unknown choice for next action. Abort.\n");
	  exit (0);
	}
    }

  //When it gets here an escpae has taken place. That's it done.

  *which_out = 2;
  return (0);

  //All done. 
}





/************************************************************
                                    Imperial College London
Synopsis:

       get_kpkt_f returns the specific luminosity in the band needed for the computation of the
       spectrum. It gets the total energy radiated by the process k-packet -> r-packet in the
       required wavelength range.


Arguments:   
       WindPtr w                   the ptr to the structure defining the wind

                                  
           

Returns:   The energy radiated by the process k-packet -> r-packet in the wind in the 
           wavelength range required for the specrum calculation.



Description:


Notes:  
        
         
History:
          June 04 SS - coding began
	06may	ksl	57+ -- Recoded to use plasma array

************************************************************/

double
get_kpkt_f ()
{
  int n;
  double lum;

  lum = 0.0;

  for (n = 0; n < NPLASMA; n++)
    {
      lum += plasmamain[n].kpkt_emiss;
    }


  return (lum);
}

/* All done. */

/************************************************************
                                    Imperial College London
Synopsis:

       get_matom_f returns the specific luminosity in the band needed for the computation of the
       spectrum. It gets the total energy radiated by the deactivation of macro atoms in the
       required wavelength range.


Arguments:   
       WindPtr w                   the ptr to the structure defining the wind

                                              

Returns:   The energy radiated by the deactivation of macro atoms in the wind in the 
           wavelength range required for the specrum calculation.



Description:


Notes:  
        
         
History:
          June 04 SS - coding began
          Sep  04 SS - significant modification to improve the treatment of macro
                       atoms in spectral synthesis steps. 
	06may	ksl	57+ -- Recoded to use plasma structure

************************************************************/

double
get_matom_f ()
{
  int n, m;
  int mm, ss;
  double lum;
  int level_emit[NLEVELS_MACRO], kpkt_emit;
  int n_tries;
  struct photon ppp;
  double contribution;
  int nres, which_out;

  which_out = 0;
  lum = 0.0;
  n_tries = 1000;
  geo.matom_radiation = 0;

  for (n = 0; n < NPLASMA; n++)
    {
      for (m = 0; m < nlevels_macro + 1; m++)
	{
	  if (m == nlevels_macro
	      || (m < nlevels_macro && macromain[n].matom_abs[m] > 0))
	    {


	      /* We know "matom_abs" is the amount of energy absobed by
	         each macro atom level in each cell. We now want to determine what fraction of
	         that energy re-appears in the frequency range we want and from which macro atom
	         level it is re-emitted (or if it appears via a k-packet). */

	      for (mm = 0; mm < NLEVELS_MACRO; mm++)
		{
		  level_emit[mm] = 0;
		}
	      kpkt_emit = 0;
	      for (ss = 0; ss < n_tries; ss++)
		{
		  if (m < nlevels_macro)
		    {
		      /* Dealing with excitation of a macro atom level. */
		      /* First find a suitable transition in which to excite the macro atom (to get it
		         in the correct starting level. */

		      /* ?? Stuart, why is this loop up to nlines, instead of nlines_macro, since we are only dealing with 
		       * macro_levels here.  It looks like the reason is because you expect n in some cases to 
		       * correspond to a photoionzation level but surely there are quicker ways to get there
		       * ksl 04dec.
		       */
		      nres = 0;
		      while ((nres < nlines))
			{
			  if (lin_ptr[nres]->nconfigu != m)
			    {
			      nres += 1;
			    }
			  else
			    {
			      break;
			    }
			}


		      if (nres > (nlines - 1))
			{
			  nres = NLINES + 1;
			  while (nres < (NLINES + 1 + ntop_phot))
			    {
			      if (phot_top[nres - NLINES - 1].uplev != m)
				{
				  nres += 1;
				}
			      else
				{
				  break;
				}
			    }
			}

		      if (nres > NLINES + ntop_phot)
			{
			  Error ("Problem in get_matom_f (1). Abort. \n");
			  exit (0);
			}

		      ppp.nres = nres;
		      //ppp.grid = n; Removed by SS May 06
		      ppp.grid = plasmamain[n].nwind;
		      ppp.w = 0;

		      macro_gov (&ppp, &nres, 1, &which_out);


		      /* Now a macro atom has been excited and followed until an r-packet is made. Now, if that 
		         r-packet is in the correct frequency range we record it. If not, we throw it away. */
		    }
		  else if (m == nlevels_macro)
		    {
		      /* kpkt case. */

		      nres = -2;	//will do 

		      ppp.nres = nres;
		      //ppp.grid = n; Removed by SS May 06
		      ppp.grid = plasmamain[n].nwind;
		      ppp.w = 0;

		      macro_gov (&ppp, &nres, 2, &which_out);

		      /* We have an r-packet back again. */
		    }


		  if (ppp.freq > em_rnge.fmin && ppp.freq < em_rnge.fmax)
		    {
		      if (which_out == 1)
			{
			  if (nres < 0)
			    {
			      Error ("Negative out from matom?? Abort.\n");
			      exit (0);
			    }

			  /* It was a macro atom de-activation. */
			  if (nres < nlines)
			    {	/* It was a line. */
			      level_emit[lin_ptr[nres]->nconfigu] += 1;
			    }
			  else
			    {
			      level_emit[phot_top[nres - NLINES - 1].
					 uplev] += 1;
			    }
			}
		      else if (which_out == 2)
			{
			  /* It was a k-packet de-activation. */
			  kpkt_emit += 1;
			}
		      else
			{
			  Error
			    ("Packet didn't emerge from matom or kpkt??? Abort. \n");
			  exit (0);
			}
		    }
		}


	      /* Now we've done all the runs for this level de-activating so we can add the contributions
	         to the level emissivities, the k-packet emissivity and the total luminosity. */

	      for (mm = 0; mm < nlevels_macro; mm++)
		{
		  contribution = 0;
		  if (m < nlevels_macro)
		    {
		      macromain[n].matom_emiss[mm] += contribution =
			level_emit[mm] * macromain[n].matom_abs[m] / n_tries;
		    }
		  else if (m == nlevels_macro)
		    {
		      macromain[n].matom_emiss[mm] += contribution =
			level_emit[mm] * plasmamain[n].kpkt_abs / n_tries;
		    }
		  lum += contribution;
		}

	      if (m < nlevels_macro)
		{
		  plasmamain[n].kpkt_emiss +=
		    kpkt_emit * macromain[n].matom_abs[m] / n_tries;
		}
	      else if (m == nlevels_macro)
		{
		  plasmamain[n].kpkt_emiss +=
		    kpkt_emit * plasmamain[n].kpkt_abs / n_tries;
		}
	    }
	}
//OLD   }
    }


  geo.matom_radiation = 1;

  return (lum);
}

/* All done. */



/************************************************************
                                    Imperial College London
Synopsis:

     photo_gen_kpkt produces photon packets to account for creating of r-packets
     by k-packets in the spectrum calculation. It should only be used once the total 
     energy emitted in this way in the wavelength range in question is well known
     (calculated in the ionization cycles).


Arguments:   
       WindPtr w                   the ptr to the structure defining the wind
       PhotPtr p;                  the ptr to the structire for the photons
       double weight;              the photon weight
       int photstart, nphot;       the number of the first photon to be generated and
                                   the total number of photons to be generated by this
                                   routine

                                           

Returns:   
      When it finishes it should have generated nphot photons from k-packet elliminations.


Description:
      This routine is closely related to photo_gen_wind from which much of the code has been copied.

Notes:  
        
         
History:
          June 04 SS - coding began
	  04aug	ksl	Modified so that uses coordinate system
	  		independent routine get_random_location
			to set location of generated kpkt
          Sep  04 SS - minor bug fixed
	06may	ksl	57+ -- Updated to use plasma structure. 
			and to elimate passing entrie wind structure

************************************************************/
int
photo_gen_kpkt (p, weight, photstart, nphot)
     PhotPtr p;
     double weight;
     int photstart, nphot;
{
  int photstop;
  int icell;
  double xlum, xlumsum;
  struct photon pp;
  int nres, esc_ptr;
  int n;
  double v[3];
  double dot ();
  int wind_n_to_ij (), stuff_v (), where_in_wind (), randvec ();
  int get_random_location ();
  double rho_test, test;
  double ztest, dvds, z, tau;
  int nnscat;
  double dvwind_ds (), sobolev ();
  int nplasma;


  // Error ("Calling photo_gen_kpkt - don't want this yet.\n");
  //exit(0);

  photstop = photstart + nphot;
  Log ("photo_gen_kpkt creates nphot %5d photons from %5d to %5d \n", nphot,
       photstart, photstop);

  for (n = photstart; n < photstop; n++)
    {
      /* locate the wind_cell in which the photon bundle originates. */

      xlum = (rand () + 0.5) / (MAXRAND) * geo.f_kpkt;

      xlumsum = 0;
      icell = 0;
      while (xlumsum < xlum)
	{
	  if (wmain[icell].vol > 0.0)
	    {
	      nplasma = wmain[icell].nplasma;
	      xlumsum += plasmamain[nplasma].kpkt_emiss;
	    }
	  icell++;
	}
      icell--;			/* This is the cell in which the photon must be generated */

      /* Now generate a single photon in this cell */
      p[n].w = weight;

      pp.freq = 0.0;
      pp.grid = icell;

      /* This following block is a bad way of doing it but it'll do as
         a quick and dirty test for now. Really kpkt should be modified to 
         do what we want in a more elegant way. */

      test = pp.freq;

      while (test > em_rnge.fmax || test < em_rnge.fmin)
	{
	  kpkt (&pp, &nres, &esc_ptr);
	  if (esc_ptr == 0)
	    {
	      test = 0.0;
	    }
	  else
	    {
	      test = pp.freq;
	    }
	}


      p[n].freq = pp.freq;
      p[n].nres = nres;

      /* The photon frequency is now known. */

      /* Determine the position of the photon in the moving frame */
/*04aug -- ksl -- 52a -- Replaced next set of lines with coordinate-sytem
 * independent call.  
 */
      get_random_location (icell, p[n].x);

      /*04 dec -- The check below should not be needed.  get_random_location already
       * does this.  The actual check is in where in wind.  I've commented this out
       * and added an error check that will kill the program if this is wrong. 54b-ksl
       */
      if (geo.disk_type == 2)
	{
	  rho_test = ((p[n].x[0] * p[n].x[0]) + (p[n].x[1] * p[n].x[1]));
	  if ((fabs (p[n].x[2]) < zdisk (rho_test))
	      && (rho_test < geo.diskrad))
	    {
	      Error ("kpkt -- Photon created in disk");
	      exit (0);
	    }

	}


      if (where_in_wind (p[n].x) != 0)
	{
	  Error ("phot_gen_kpkt: generated photon not in wind\n");
	}
      if (where_in_grid (p[n].x) < 0)
	{
	  Error ("phot_gen_kpkt: generated photon not in grid\n");
	}



      p[n].grid = icell;

      nnscat = 1;
      // Determine the direction of the photon
      // Need to allow for anisotropic emission here
      if (p[n].nres < 0 || p[n].nres > NLINES || geo.scatter_mode == 0)
	{
	  /*  It was either an electron scatter so the  distribution is isotropic, or it
	     was a resonant scatter but we want isotropic scattering anyway.  */
	  randvec (p[n].lmn, 1.0);	/* The photon is emitted isotropically */
	}
      else if (geo.scatter_mode == 1)
	{			// It was a line photon and we want anisotropic scattering

	  // -1. forces a full reinitialization of the pdf for anisotropic scattering

	  randwind (&p[n], p[n].lmn, wmain[icell].lmn);

	}
      else if (geo.scatter_mode == 2)
	{			//It was a line photon and we want the thermal trapping anisotropic model

	  ztest = 1.0;
	  z = 0.0;
	  nnscat = 0;
	  while (ztest > z)
	    {
	      nnscat = nnscat + 1;
	      randvec (p[n].lmn, 1.0);	//Get a new direction
	      ztest = (rand () + 0.5) / MAXRAND;
	      dvds = dvwind_ds (&p[n]);	//Here w is entire wind
	      tau =
		sobolev (&wmain[icell], &p[n], -1.0, lin_ptr[p[n].nres],
			 dvds);
	      if (tau < 1.0e-5)
		z = 1.0;
	      else
		z = (1. - exp (-tau)) / tau;	/* probability to see if it escapes in that direction */
	    }
	}

      p[n].nnscat = nnscat;





      /* The next two lines correct the frequency to first order, but do not result in
         forward scattering of the distribution */

      vwind_xyz (&p[n], v);
      p[n].freq /= (1. - dot (v, p[n].lmn) / C);

      p[n].istat = 0;
      p[n].tau = p[n].nscat = p[n].nrscat = 0;
      p[n].origin = PTYPE_WIND;	// Call it a wind photon

    }


  return (nphot);		/* Return the number of photons generated */


  /* All done. */
}

/************************************************************
                                    Imperial College London
Synopsis:

     photo_gen_matom produces photon packets to account for creation of r-packets
     by deactivation of macro atoms in the spectrum calculation. It should only be used 
     once the total energy emitted in this way in the wavelength range in question is well known
     (calculated in the ionization cycles).


Arguments:   
       WindPtr w                   the ptr to the structure defining the wind
       PhotPtr p;                  the ptr to the structire for the photons
       double weight;              the photon weight
       int photstart, nphot;       the number of the first photon to be generated and
                                   the total number of photons to be generated by this
                                   routine

                                           

Returns:   
      When it finishes it should have generated nphot photons from macro atom deactivations.


Description:
      This routine is closely related to photo_gen_kpkt from which much of the code has been copied.

Notes:  
        
         
History:
          June 04 SS - coding began
          Aug  04 SS - modified as for gen_kpkt above by KSL
          Sep  04 SS - minor bug fixed
	06may	ksl	57+ -- Initial adaptation to plasma structure

************************************************************/
int
photo_gen_matom (p, w, weight, photstart, nphot)
     WindPtr w;			// Wind is entire wind here
     PhotPtr p;
     double weight;
     int photstart, nphot;
{
  int photstop;
  int icell;
  double xlum, xlumsum;
  struct photon pp;
  int nres;
  int n;
  double v[3];
  double dot ();
  int wind_n_to_ij (), stuff_v (), where_in_wind (), randvec ();
  int emit_matom ();
  double rho_test, test;
  int upper;
  double ztest, dvds, z, tau;
  int nnscat;
  double dvwind_ds (), sobolev ();
  int nplasma;



  photstop = photstart + nphot;
  Log ("photo_gen_matom creates nphot %5d photons from %5d to %5d \n", nphot,
       photstart, photstop);

  for (n = photstart; n < photstop; n++)
    {
      /* locate the wind_cell in which the photon bundle originates. And also decide which of the macro
         atom levels will be sampled (identify that level as "upper"). */

      xlum = (rand () + 0.5) / (MAXRAND) * geo.f_matom;

      xlumsum = 0;
      icell = 0;
      upper = 0;
      while (xlumsum < xlum)
	{
	  if (w[icell].vol > 0.0)
	    {
	      nplasma = w[icell].nplasma;
	      xlumsum += macromain[nplasma].matom_emiss[upper];
	      upper++;
	      if (upper == nlevels_macro)
		{
		  upper = 0;
		  icell++;
		}
	    }
	  else
	    {
	      icell++;
	      upper = 0;
	    }
	}

      if (upper == 0)
	{
	  /* If upper = 0 at this point it means the loop above finished with an increment of
	     icell and setting upper = 0. In such a case the process we want was the LAST process
	     in the previous cell. */
	  icell--;
	  upper = nlevels_macro;
	}

      upper--;			/* This leaves the macro atom level that deactivaties. */

      /* Now generate a single photon in this cell */
      p[n].w = weight;

      pp.freq = 0.0;
      pp.grid = icell;

      /* This following block is a bad way of doing it but it'll do as
         a quick and dirty test for now. Really kpkt should be modified to 
         do what we want in a more elegant way. */

      test = pp.freq;

      while (test > em_rnge.fmax || test < em_rnge.fmin)
	{

	  /* Call routine that will select an emission process for the
	     deactivating macro atom. If is deactivates outside the frequency
	     range of interest then ignore it and try again. SS June 04. */

	  emit_matom (w, &pp, &nres, upper);

	  test = pp.freq;
	}


      p[n].freq = pp.freq;
      p[n].nres = nres;


      /* The photon frequency is now known. */

      /* Determine the position of the photon in the moving frame */
      get_random_location (icell, p[n].x);

      /*04 dec -- The check below should not be needed.  get_random_location already
       * does this.  The actual check is in where in wind.  I've commented this out
       * and added an error check that will kill the program if this is wrong. Check
       * should be deleted once we are sure as it simply wastes time.  54b-ksl
       */
      if (geo.disk_type == 2)
	{
	  rho_test = ((p[n].x[0] * p[n].x[0]) + (p[n].x[1] * p[n].x[1]));
	  if ((fabs (p[n].x[2]) < zdisk (rho_test))
	      && (rho_test < geo.diskrad))
	    {
	      Error ("kpkt -- Photon created in disk");
	      exit (0);
	    }

	}



      if (where_in_wind (p[n].x) != 0)
	{
	  Error ("phot_gen_matom: generated photon not in wind\n");
	}
      if (where_in_grid (p[n].x) < 0)
	{
	  Error ("phot_gen_matom: generated photon not in grid\n");
	}

      p[n].grid = icell;


      // Determine the direction of the photon
      // Need to allow for anisotropic emission here
      nnscat = 1;
      if (p[n].nres < 0 || p[n].nres > NLINES || geo.scatter_mode == 0)
	{
	  /*  It was either an electron scatter so the  distribution is isotropic, or it
	     was a resonant scatter but we want isotropic scattering anyway.  */
	  randvec (p[n].lmn, 1.0);	/* The photon is emitted isotropically */
	}
      else if (geo.scatter_mode == 1)
	{			// It was a line photon and we want anisotropic scattering

	  // -1. forces a full reinitialization of the pdf for anisotropic scattering

	  randwind (&p[n], p[n].lmn, w[icell].lmn);

	}
      else if (geo.scatter_mode == 2)
	{			//It was a line photon and we want the thermal trapping anisotropic model

	  ztest = 1.0;
	  z = 0.0;
	  nnscat = 0;
	  while (ztest > z)
	    {
	      nnscat = nnscat + 1;
	      randvec (p[n].lmn, 1.0);	//Get a new direction
	      ztest = (rand () + 0.5) / MAXRAND;
	      dvds = dvwind_ds (&p[n]);	//Here w is entire wind
	      tau =
		sobolev (&w[icell], &p[n], -1.0, lin_ptr[p[n].nres], dvds);
	      if (tau < 1.0e-5)
		z = 1.0;
	      else
		z = (1. - exp (-tau)) / tau;	/* probability to see if it escapes in that direction */
	    }
	}
      p[n].nnscat = nnscat;


      /* The next two lines correct the frequency to first order, but do not result in
         forward scattering of the distribution */

      vwind_xyz (&p[n], v);
      p[n].freq /= (1. - dot (v, p[n].lmn) / C);

      p[n].istat = 0;
      p[n].tau = p[n].nscat = p[n].nrscat = 0;
      p[n].origin = PTYPE_WIND;	// Call it a wind photon

    }


  return (nphot);		/* Return the number of photons generated */


  /* All done. */
}

/************************************************************
                                    Imperial College London
Synopsis:
	emit_matom is a scaled down version of matom which deals with the emission due
        to deactivating macro atoms in the detailed spectrum part of the calculation.


Arguments:

       WindPtr w                   the ptr to the structure defining the wind
       PhotPtr p                   the packet at the point of activation
       int upper                   the upper level that we deactivate from

Returns:
       int nres                    the process by which deactivation occurs
       PhotPtr p                   the packet following deactivation
    
Description:


Notes: 

History:
        June 04  -   coding began
	06may	ksl	57+ -- Initial adaptation to plasma structure

************************************************************/

int
emit_matom (w, p, nres, upper)
     WindPtr w;
     PhotPtr p;
     int *nres;
     int upper;
{
  struct lines *line_ptr;
  struct topbase_phot *cont_ptr;
  int uplvl;
  double alpha_sp ();
  double eprbs[NBBJUMPS + NBFJUMPS];
  double penorm;
  double threshold, run_tot;
  double sp_rec_rate;
  int n, m;
  int nbbd, nbfd;
  double t_e, ne;
  double bb_cont;
  WindPtr one;
  PlasmaPtr xplasma;


  one = &w[p->grid];		//This is to identify the grip cell in which we are
  xplasma = &plasmamain[one->nplasma];

  t_e = xplasma->t_e;		//electron temperature 
  ne = xplasma->ne;		//electron number density


  /* The first step is to identify the configuration that has been excited. */

  uplvl = upper;

  /* Unlike the proper matom routine, we don't want any jumps, only deactivations here. */


  /*  The excited configuration is now known. Now compute all the probabilities of deactivation
     /jumping from this configuration. Then choose one. */


  nbbd = config[uplvl].n_bbd_jump;	//store these for easy access -- number of bb downward jumps
  nbfd = config[uplvl].n_bfd_jump;	// number of bf downared jumps from this transition


  // Start by setting everything to 0

  m = 0;			//m counts the total number of possible ways to leave the level
  penorm = 0.0;			//stores the total emission probability

  for (n = 0; n < nbbd + nbfd; n++)
    {
      eprbs[n] = 0;		//stores the individual emission probabilities SS
    }
  /* Finished zeroing. */


  /* bb */

  /* First downward jumps. */

  for (n = 0; n < nbbd; n++)
    {
      line_ptr = &line[config[uplvl].bbd_jump[n]];
      /* Since we are only interested in making an r-packet here we can (a) ignore collisional
         deactivation and (b) ignore lines outside the frequency range of interest. */
      if ((line_ptr->freq > em_rnge.fmin) && (line_ptr->freq < em_rnge.fmax))	// correct range
	{
	  bb_cont = (a21 (line_ptr) * p_escape (line_ptr, xplasma));

	  eprbs[m] = bb_cont * (config[uplvl].ex - config[line[config[uplvl].bbd_jump[n]].nconfigl].ex);	//energy difference

	  if (eprbs[m] < 0.)	//test (can be deleted eventually SS)
	    {
	      Error ("Negative probability (matom, 2). Abort.");
	      exit (0);
	    }

	  penorm += eprbs[m];
	}
      m++;
    }

  /* bf */
  for (n = 0; n < nbfd; n++)
    {
      cont_ptr = &phot_top[config[uplvl].bfd_jump[n]];	//pointer to continuum

      /* If the edge is above the frequency range we are interested in then we need not consider this
         bf process. */
      if (cont_ptr->freq[0] < em_rnge.fmax)	//means that it may contribute
	{
	  sp_rec_rate = alpha_sp (cont_ptr, xplasma, 0);
	  eprbs[m] = sp_rec_rate * xplasma->ne * (config[uplvl].ex - config[phot_top[config[uplvl].bfd_jump[n]].nlev].ex);	//energy difference
	  if (eprbs[m] < 0.)	//test (can be deleted eventually SS)
	    {
	      Error ("Negative probability (matom, 4). Abort.");
	      exit (0);
	    }
	  penorm += eprbs[m];
	}
      m++;
    }



  /* No need for upwards jumps - they don't provide emission mechanisms. */

  /* Probabilities of emission (e) are now known. The normalisation factor pnorm has
     also been recorded. the integer m now gives the total number of possibilities too. 
     now select what happens next. Start by choosing the random threshold value at which the
     event will occur. */

  threshold = ((rand () + 0.5) / MAXRAND);

  run_tot = 0;
  n = 0;
  threshold = threshold * penorm;	//normalise to total emission prob.
  while (run_tot < threshold)
    {
      run_tot += eprbs[n];
      n++;
    }
  n = n - 1;
  /* n now identifies the jump that occurs - now set nres for the return value. */
  if (n < nbbd)
    {				/* bb downwards */
      line_ptr = &line[config[uplvl].bbd_jump[n]];	//pointer for the bb transition
      *nres = line[config[uplvl].bbd_jump[n]].where_in_list;
      p->freq = line[config[uplvl].bbd_jump[n]].freq;
      /* This is the comoving frequency - changed to rest frequency by doppler */
    }
  else if (n < (nbbd + nbfd))
    {				/* bf downwards jump */
      *nres = config[uplvl].bfd_jump[n - nbbd] + NLINES + 1;
      /* continuua are indicated by nres > NLINES */
      p->freq = phot_top[config[uplvl].bfd_jump[n - nbbd]].freq[0]
	- (log (1. - (rand () + 0.5) / MAXRAND) * xplasma->t_e / H_OVER_K);
      /* Co-moving frequency - changed to rest frequency by doppler */
      /*Currently this assumed hydrogenic shape cross-section - Improve */
    }
  else
    {
      Error
	("Trying to emitt from Macro Atom but no available route (emit_matom). Abort.");
      exit (0);
    }
  return (0);
}

/* The frequency and the value of nres have been set correctly. All done. */

/********************************************************************************/


/******************************************************************************/

/* q_ioniz. This returns the collisional ionization co-efficient
Calculated following equation 5-79 of Mihalas.
*/

double
q_ioniz (cont_ptr, electron_temperature)
     struct topbase_phot *cont_ptr;
     double electron_temperature;
{
  double coeff;
  double gaunt;
  double u0;

  u0 = cont_ptr->freq[0] * H_OVER_K / electron_temperature;
  gaunt = 0.1;			//for now - from Mihalas for hydrogen

  coeff = 1.55e13 / sqrt (electron_temperature) * gaunt * cont_ptr->x[0] *
    exp (-1. * u0) / u0;

  return (coeff);
}

/******************************************************************************/

/* q_recomb. This returns the collisional recombination co-efficient
Calculated from inverse of q_ioniz.
*/

double
q_recomb (cont_ptr, electron_temperature)
     struct topbase_phot *cont_ptr;
     double electron_temperature;
{
  double coeff;
  double root_etemp;
  double q_ioniz ();

  root_etemp = sqrt (electron_temperature);
  coeff = 2.07e-16 / (root_etemp * root_etemp * root_etemp);
  coeff *= exp (cont_ptr->freq[0] * H_OVER_K / electron_temperature);
  coeff *= config[cont_ptr->nlev].g / config[cont_ptr->uplev].g;
  coeff *= q_ioniz (cont_ptr, electron_temperature);

  return (coeff);
}

/* 

!! Stuart, in 57h -- I have replaced this routine.  See bb_pdf_get_rand_limit in bb.c
This routine a number of out of range photons becasue you did not make use of
pdf->limit1 and 2 to restict the potential frequecies.  Please, check that routine,
and if you are happy delete this one.  ksl -- 060806



The following is an alternative to pdf_get_rand_limit that will sample a
   pdf assuming that the function is cubic rather than linear. Should do better
   at making a smooth continuum in the IR. 

History
	06aug	ksl	57h -- This routine has an error in it,  It does not properly
			limit x to lie between the pdf_limits, and was generating
			many photons outside the desired energy intervals.  Also
			while correct at low limit.  It would seem to be horrible
			for the high T limit, where you are being killed by the
			exponential.  
*/
//OLD57h double
//OLD57h pdf_get_rand_limit_bbir (pdf)
//OLD57h      PdfPtr pdf;
//OLD57h {
//OLD57h   double x, r;
//OLD57h   int i;
//OLD57h   double choice;

//OLD57h   r = rand () / MAXRAND;       /* r must be slightly less than 1 */
//OLD57h   r = r * pdf->limit2 + (1. - r) * pdf->limit1;

//OLD57h   i = r * NPDF;

//OLD57h   while (pdf->y[i + 1] < r && i < NPDF - 1)
//OLD57h     i++;
//OLD57h   while (pdf->y[i] > r && i > 0)
//OLD57h     i--;

  //  q = (r - pdf->y[i]) / (pdf->y[i + 1] - pdf->y[i]);

  //x = pdf->x[i] * (1. - q) + pdf->x[i + 1] * q;

//OLD57h   choice = rand () / MAXRAND;

//OLD57h   x = choice * pow (pdf->x[i], 3.) + (1. - choice) * pow (pdf->x[i + 1], 3.);
//OLD57h   x = pow (x, (1. / 3.));

//OLD57h   return (x);

//OLD57h }
