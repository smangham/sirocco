/***********************************************************
                                       Space Telescope Science Institute

 Synopsis: balance is a test program to allow one to understand how and whether various
     radiative process approach lte
   
 Arguments:		

Returns:
 
 
Description:	

	G:   An option to calculate the abundances of the light elements in a variety of situations, for comparison
	     with other routines like Cloudy
		
Notes:
	Insofar as possible the program uses python routines directly

History:
	98jul	ksl	Coding began during a visit to Oxford
	01mar	ksl	Verified program still runs basically, but note
			that a number of the underlying routines have
			changed.
	01sept	ksl	Attempted to really sort out some of the various options
			that no longer mattered and eliminate these sections
			of balance, and added options for levels.
	02jun	ksl	Added an option to set the total heating
	02jun	ksl	Added a subroutine to calculate what the derived values
			for the radiation field would be after photon passage.
 
**************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
#include "balance_templates.h"



#define LINELENGTH 132


int
main ()
{

  WindPtr w;
  PhotPtr p;
  double t_r, t_e;
  double nh, heat;
  double weight, vol, vmax;
  char choice;
  char atomic_filename[LINELENGTH];
  int lmode;
  int freq_sampling;
  int check_field ();
  int do_detail ();
  int cmode;

  NPHOT = 100000;

/* These next lines define the gird size and give the push through distance */
  ndim = 3;
  mdim = 3;
  NDIM = ndim;
  MDIM = mdim;
  NDIM2 = ndim * mdim;
  dfudge = 0.0001;
  DFUDGE = dfudge;
/* End of definition of wind array size */

/* Get the atomic data */

  strcpy (atomic_filename, "atomic/standard");
  rdstr ("Atomic_data", atomic_filename);
  get_atomic_data (atomic_filename);

/* Create the photon array */
  w = (WindPtr) calloc (sizeof (wind_dummy), NDIM2);

  if (w == NULL)
    {
      Error
	("There is a problem in allocating memory for the wind structure\n");
      exit (0);
    }

  p = (PhotPtr) calloc (sizeof (p_dummy), NPHOT);
  if (p == NULL)
    {
      Log ("There was a problem allocating p\n");
      exit (0);
    }

/* Create a set of inital choices */

  geo.line_mode = 2;
  geo.ioniz_mode = 0;
  geo.rt_mode = 0;
  geo.pl_nh = nh = 1e10;
  geo.pl_t_r = t_r = 20000;
  geo.pl_t_e = t_e = 20000;
  geo.pl_w = weight = 1;
  geo.pl_vol = vol = 1.0;
  geo.pl_vmax = vmax = 1e6;
  freq_sampling = 1;

  choice = 'l';

  pl_wind_define (w);
  levels (w, 1);
  print_levels (w);
/* Completed the initialization */

a:
  Log
    ("\nnh %8.2g t_r %8.2g t_e %8.2g weight %8.2e vol %8.2g vmax %8.2g line_mode %2d ioniz_mode %d rt_mode %d freq_sampling %d\n\n",
     nh, t_r, t_e, weight, vol, vmax, geo.line_mode, geo.ioniz_mode,
     geo.rt_mode, freq_sampling);
  Log
    ("h=Nh, r=t_r, e=t_e, w=weight, v=vol, s=vmax, p=line_mode, k=rad_trans_mode\n");
  Log ("i=ionization mode, j=calculate heating/cooling and abundances\n");
  Log ("l=lte(force t_e=t_r w=1),  a=absolute_lte(ions+levels), A=levels\n");
  Log
    ("b=frequency sampling d=dumb, f=find_te, ,g=go  H=Set heating, D -- modify to detailed rates\n");
  rdchar ("Choice", &choice);
  Log_silent ("Choice %c\n", choice);

/* Initialize the wind structure */
  w[0].w = weight;
  w[0].t_e = t_e;
  w[0].t_r = t_r;
  w[0].vol = vol;
//      w[0].dvds_ave=1.0;
  w[0].rho = nh / rho2nh;
  w[0].gain = 0.5;

  switch (choice)
    {


    case 'a':			/* Force Saha ion densities, and LTE level populations, and set t_e=t_r and weight to 1
				   This differs from the case l (lte option) in that LTE level populations are forced explicitly */

      rddoub ("LTE temp", &t_e);
      Log_silent ("LTE temp %8.2g\n", t_e);

      geo.pl_nh = nh;
      geo.pl_t_r = t_r = t_e;
      geo.pl_t_e = t_e;
      geo.pl_w = weight = 1.;
      geo.pl_vol = vol;
      geo.pl_vmax = vmax;
//              w[0].dvds_ave=1.0;
      pl_wind_define (w);

      Log
	("Calculating heating and cooling with t_e=t_r, LTE ionization & LTE exited state populations\n");
      absolute_saha (&w[0], p, nh, t_e, t_e, weight, freq_sampling);
      break;

    case 'A':
// Calculate the levels with various options
      rdint ("Level_calc(0=onthespot,1=lte)", &lmode);

      levels (w, lmode);
      print_levels (w);
      break;

    case 'b':			// freqmode
      rdint ("frequency_sampling(0=uniform)", &freq_sampling);
      if (freq_sampling == 0)
	Log ("OK -- using normal sampling with no required photon numbers\n");
      else
	{
	  freq_sampling = 1;
	  Log
	    ("OK-- setting energy bands to assure substantial photons above ionization edges\n");
	}

    case 'd':
// Increment t_e up or down to attempt to match the current value of the heating.  
// It's very dumb though because the heating is not automatically calculated.

      dumb_step (&w[0], &t_e);
      break;

    case 'D':
// Modify the current abundances to produce detail balance, and then recalculate
// the heating and cooling without modifying the abundances
//
      do_detail (&w[0]);
      xsaha (&w[0], p, nh, t_r, t_e, weight, -1, freq_sampling);
      break;
    case 'e':			// Set the electron temperature

      rddoub ("t_e", &t_e);
      Log_silent ("t_e %8.2g\n", t_e);
      break;

    case 'f':			//Adjust t_e to balance heating and cooling without changing the abundances

      Log
	("Adjusting t_e to balance heating and cooling without changing the abundances\n");
      find_te (&w[0]);
      t_e = w[0].t_e;
      break;

    case 'g':			// Adjust t_e and abundances using python's one_shot routine.  
/* The first time this routine is called after t_r or nh have been called a saha density
distribution is set up based on concentrations, and then one_shot is called.  If this
case is repeatedly executed it should eventually match the abundances and temperatures.
Note that the ionization mode is not changed by this routine.  */
      if (geo.ioniz_mode == 0)
	cmode = 3;
      if (geo.ioniz_mode == 4)
	cmode = 4;
      cycle (&w[0], p, nh, t_r, t_e, weight, cmode, freq_sampling);
      t_e = w[0].t_e;
      break;

    case 'G':			// Calculate lots of cases at once 
      if (geo.ioniz_mode == 0)
	cmode = 3;
      if (geo.ioniz_mode == 4)
	cmode = 4;
      multicycle (&w[0], p, cmode, freq_sampling);
      break;

    case 'h':			//Set the hydrogen density

      rddoub ("NH", &nh);
      Log_silent ("NH %8.2g\n", nh);
      break;

    case 'H':
      rddoub ("Total_heating", &heat);
      w[0].heat_tot = heat;
      break;

    case 'i':			//Set the mode for which ionization abundances will be calculated, e.g. saha abundances

      printf
	("0=on_the_spot,1=saha,2=fix to std vals, 3=_on_the_spot (with one t update), 4=detailed, 10=same as 3, but ground mult.\n");
      rdint ("ioniz_mode", &geo.ioniz_mode);
      geo.partition_mode = -1;

      if (geo.ioniz_mode == 10)	// Use recalc but only ground state multiplicities
	{
	  geo.ioniz_mode = 3;
	  geo.partition_mode = 0;
	}

      Log_silent ("Ionization mode %d", geo.ioniz_mode);
      break;

    case 'I':
      multi_concen (&w[0]);
      break;

    case 'j':
/* Calculate abundances and then heating and cooling with current ionization mode, etc.
but do not update basic parameters aftewards.  Multiple calls to j should give identical
answers, except for statistical changes */


      geo.pl_nh = nh;
      geo.pl_t_r = t_r;
      geo.pl_t_e = t_e;
      geo.pl_w = weight;
      geo.pl_vol = vol;
      geo.pl_vmax = vmax;
      pl_wind_define (w);



/* Assure that www is initialized before going into xsaha */
      w[0].t_e = t_e;
      w[0].t_r = t_r;
      w[0].w = weight;



/* Currently (Sept 01), python takes ioniz_mode to mean
0  on-the-spot approximation using existing t_e
1  LTE
2  Hardwired 
3  On the spot, with one_shot at updating t_e before calculating densities

*/
      xsaha (&w[0], p, nh, t_r, t_e, weight, geo.ioniz_mode, freq_sampling);


      check_field (&w[0]);

      break;

    case 'k':			//Set the radiative tranfer mode.  USED ONLY BY BALANCE 

/* Simple here means that the line shape is set to a simple square wave
   MC means we have run the MC code and so lines are only absorbed as closely to
	the way they are handled in python as possible, i.e. only if tau exceeds a
	defined value
   Simple sobolev calculates the heating using the Sobolev approx, but without the 
   MC aspects of deciding whether any part of a photon was absorbed
*/
      rdint ("Radiative transfer(0=simple,1=MC,2=simple_sobolev)",
	     &geo.rt_mode);
      Log_silent
	("Radiative transfer(0=simple,1=MC,2=simple_sobobolev): %d\n",
	 geo.rt_mode);
      break;

    case 'l':			/*Force Saha ion abundances and t_e=t_r and weight to 1 and then calculate heating and cooling.
				   In this option, the levels are calculated with the same routines as in python.  The results should be identical
				   to case a unless there is an error somewhere */
      rddoub ("LTE temp", &t_e);
      Log_silent ("LTE temp %8.2g\n", t_e);

      geo.pl_nh = nh;
      geo.pl_t_r = t_r = t_e;
      geo.pl_t_e = t_e;
      geo.pl_w = weight = 1.0;
      geo.pl_vol = vol;
      geo.pl_vmax = vmax;
      pl_wind_define (w);
//              w[0].dvds_ave=1.0;
      geo.ioniz_mode = 1;


      Log
	("Calculating heating and cooling with t_e=t_r and Saha=LTE ionization\n");
/* Assure that www is initialized before going into xsaha */
      w[0].t_e = t_e;
      w[0].t_r = t_e;
      w[0].w = weight;

      xsaha (&w[0], p, nh, t_e, t_e, weight, 1, freq_sampling);
      break;

    case 'p':			// Set the line transfer mode mode

      rdint ("line_mode", &geo.line_mode);
      Log_silent ("line_mode %d\n", geo.line_mode);
      break;

    case 'r':			//Set the radiation temperature

      rddoub ("t_r", &t_r);
      Log_silent ("t_r %8.2g\n", t_r);
      break;

    case 's':			//Set the radiation temperature

      rddoub ("vmax", &vmax);
      geo.pl_nh = nh;
      geo.pl_t_r = t_r;
      geo.pl_t_e = t_e;
      geo.pl_w = weight;
      geo.pl_vol = vol;
      geo.pl_vmax = vmax;
      pl_wind_define (w);
      Log_silent ("vmax %8.2g\n", vmax);
      break;

    case 'T':			//This is just for testing of routines
      init_bands (t_e, 0.0, 1e50, 1, &xband);
      init_freebound (1e3, 1e6, xband.f1[0], xband.f2[0]);
      break;


    case 'v':			//Set the volume

      rddoub ("volume", &vol);
      w[0].vol = vol;
      geo.pl_nh = nh;
      geo.pl_t_r = t_r;
      geo.pl_t_e = t_e;
      geo.pl_w = weight;
      geo.pl_vol = vol;
      geo.pl_vmax = vmax;
      pl_wind_define (w);
      Log_silent ("vol %8.2g\n", vol);
      break;

    case 'w':			//Set the radiative weight

      rddoub ("w", &weight);
      Log_silent ("weight %8.2g\n", weight);
      break;



    case 'q':
      goto b;
      break;

    default:
      goto a;
      break;

    }
  goto a;

b:exit (0);

}



int
multicycle (www, p, mode, freq_sampling)
     WindPtr www;
     PhotPtr p;
     int mode;
     int freq_sampling;

{
  FILE *fopen (), *fptr;
  double nh;
  double t_r;
  double t_e;
  double weight;
  double tot;

  double t_rmin, t_rmax, t_rdelt;
  double ww, w_min, w_max, w_delt;

  int n, m, first, last, nn;

  nh = 1e12;
  t_r = 30000;
  t_e = 20000;
  weight = 0.1;

  t_rmin = 2000;
  t_rmax = 100000;
  t_rdelt = 2000;

//  t_rmin = 50000;
//  t_rmax = 100000;
//  t_rdelt = 10000;

  w_max = 0;
  w_min = -6;
  w_delt = 2;

  fptr = fopen ("bal_summary.out", "w");

  for (ww = w_max; ww >= w_min; ww -= w_delt)	// Note w starts at 1
    {
      weight = pow (10, ww);
      for (t_r = t_rmin; t_r <= t_rmax; t_r += t_rdelt)
	{
	  www->gain = 0.5;
	  www->dt_e_old = 0.0;
	  www->dt_e = 0.0;
	  www->t_e = www->t_e_old = 0.9 * t_r;	//Lucy guess
//        for (n = 0; n < 5; n++)
	  for (n = 0; n < 10; n++)
	    {
	      Log ("multiCycle %d\n", n);

	      www->rho = nh / rho2nh;
	      cycle (www, p, nh, t_r, t_e, weight, mode, freq_sampling);
	      t_e = www->t_e;

	      Log ("Mo  %2d %.1f %.0f %.1f %8.0f %5.2f\n", n, log10 (nh),
		   t_r, (ww), t_e, www->heat_tot / www->lum_rad);

	    }

	  fprintf (fptr, "Mod   %.1f %.0f %.1f %8.0f %5.2f \n", log10 (nh),
		   t_r, (ww), t_e, www->heat_tot / www->lum_rad);
	  Log ("Mod    %.1f %.0f %.1f %8.0f %5.2f\n", log10 (nh), t_r, (ww),
	       t_e, www->heat_tot / www->lum_rad);

	  for (nn = 0; nn < 5; nn++)
	    {
	      first = ele[nn].firstion;
	      last = first + ele[nn].nions;
	      fprintf (fptr, "%-5s ", ele[nn].name);
	      Log ("%-5s ", ele[nn].name);
	      tot = 0;
	      for (m = first; m < last; m++)
		tot += www->density[m];
	      for (m = first; m < last; m++)
		{
		  fprintf (fptr, " %8.2e", log10 (www->density[m] / tot));
		  Log (" %8.2e", log10 (www->density[m] / tot));
		}
	      fprintf (fptr, "\n");
	      Log ("\n");
	    }
	}
    }
  fclose (fptr);
  return (0);
}


int
print_levels (w)
     WindPtr w;
{
  int n, m;
  for (n = 0; n < nions; n++)
    {
      if (ion[n].nlte > 0)
	{
	  m = 0;
	  while (ion[n].z != ele[m].z)
	    m++;
	  printf ("%-10s ", ele[m].name);

	  for (m = 0; m < ion[n].nlte; m++)
	    printf ("%8.2e ", w->levden[ion[n].nlte_first + m]);
	  printf ("\n");
	}
    }
  return (0);
}

int
multi_concen (www)
     WindPtr www;

{
  FILE *fopen (), *fptr;
  double nh;
  double t_r;
  double t_e;
  double weight;
  double tot;

  double t_rmin, t_rmax, t_rdelt;
  double ww, w_min, w_max, w_delt;

  int n, m, first, last, nn;

  nh = 1e12;
  t_r = 30000;
  t_e = 20000;
  weight = 0.1;

  t_rmin = 10000;
  t_rmax = 100000;
  t_rdelt = 2000;
//  t_rmin = 50000;
//  t_rmax = 100000;
//  t_rdelt = 10000;

  w_max = 0;
  w_min = -6;
  w_delt = 2;

  fptr = fopen ("con_summary.out", "w");

  for (ww = w_max; ww >= w_min; ww -= w_delt)	// Note w starts at 1
    {
      www->w = pow (10, ww);
      for (t_r = t_rmin; t_r <= t_rmax; t_r += t_rdelt)
	{
	  www->gain = 0.5;
	  www->dt_e_old = 0.0;
	  www->dt_e = 0.0;
	  www->t_e = www->t_e_old = 0.9 * t_r;	//Lucy guess
//        for (n = 0; n < 5; n++)
	  for (n = 0; n < 10; n++)
	    {
	      Log ("multiCycle %d\n", n);

	      www->rho = nh / rho2nh;
//            cycle (www, p, nh, t_r, t_e, weight, freq_sampling);
	      nebular_concentrations (www, 2);
	      t_e = www->t_e;

	      Log ("Mo  %2d %.1f %.0f %.1f %8.0f %5.2f\n", n, log10 (nh),
		   t_r, (ww), t_e, www->heat_tot / www->lum_rad);

	    }

	  fprintf (fptr, "Mod   %.1f %.0f %.1f %8.0f %5.2f \n", log10 (nh),
		   t_r, (ww), t_e, www->heat_tot / www->lum_rad);
	  Log ("Mod    %.1f %.0f %.1f %8.0f %5.2f\n", log10 (nh), t_r, (ww),
	       t_e, www->heat_tot / www->lum_rad);

	  for (nn = 0; nn < 5; nn++)
	    {
	      first = ele[nn].firstion;
	      last = first + ele[nn].nions;
	      fprintf (fptr, "%-5s ", ele[nn].name);
	      Log ("%-5s ", ele[nn].name);
	      tot = 0;
	      for (m = first; m < last; m++)
		tot += www->density[m];
	      for (m = first; m < last; m++)
		{
		  fprintf (fptr, " %8.2e", log10 (www->density[m] / tot));
		  Log (" %8.2e", log10 (www->density[m] / tot));
		}
	      fprintf (fptr, "\n");
	      Log ("\n");
	    }
	}
    }
  fclose (fptr);
  return (0);
}

int
check_field (one)
     WindPtr one;
{

  double trad;
  one->ave_freq /= one->j;
  one->j /= (4. * PI * one->vol);
  one->t_r = H * one->ave_freq / (BOLTZMANN * 3.832);
  one->w = PI * one->j / (STEFAN_BOLTZMANN * trad * trad * trad * trad);

  Log ("Check_field: trad %.2f  w %.2e\n", trad, one->w);
  return (0);
}

int
do_detail (one)
     WindPtr one;
{
  double newden[NIONS];
  int n, nelem;

  for (nelem = 0; nelem < nelements; nelem++)
    {
      detailed_balance (one, nelem, newden);
      wind_update_after_detailed_balance (one, nelem, newden);

    }
  for (n = 0; n < NIONS; n++)
    {
      one->density[n] = newden[n];
    }
  return (0);
}
