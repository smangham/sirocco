

/***********************************************************
                                       Space Telescope Science Institute

Synopsis:  This file contains the collected the various routines of the program
	which are associated with reading the binary versions of the kurucz files.


Arguments:



Returns:

 
Description:	
	

Notes:

	The use of these binary files is archaic and the routines ought to be removed
	from Python and full use made of the hubeny routines.  As far as I am aware,
	the only reason for these routins is in order to use the binary versions of
	the kurucz files.

	PLEASE maintain compatibility with the kurucz routines in the program cv, which
	can probably be found in ~/progs/disk and elsewhere!!!!

History:
	01aug	ksl	Added emittance_kurucz
 
**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "atomic.h"
#include "python.h"

double old_t, old_g, old_freqmin, old_freqmax;
double w_kurucz[683], f_kurucz[683];
struct Pdf pdf_kurucz;
double kjump[] = { 913.8 };



/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:

double one_kurucz(t,g,freqmin,freqmax) gets a photon frequency from the Kurucz grid 

Arguments:

	double t,g,			temperature and gravity of star in Kurucz grid
	freqmin,freqmax;		minimum and maximum frequency of interest


Returns:

 
Description:	
	

Notes:

History:
	97aug	ksl 	coded as part of python project.
	01jul	ksl	Moved into this file so all the kurucz lated 
			routines are together
 
**************************************************************/

double
one_kurucz (t, g, freqmin, freqmax)
     double t, g, freqmin, freqmax;
{
  double lambdamin, lambdamax;
  double f;
  double pdf_get_rand ();
  int kurucz_get_model ();

  if (old_t != t || old_g != g || old_freqmin != freqmin
      || old_freqmax != freqmax)
    {				/* Then we must initialize */
      kurucz_get_model (t, g, w_kurucz, f_kurucz);
      /*  Get_model returns wavelengths in Ang and flux in ergs/cm**2/Ang */
      lambdamin = C * 1e8 / freqmax;
      lambdamax = C * 1e8 / freqmin;
      if (pdf_gen_from_array
	  (&pdf_kurucz, w_kurucz, f_kurucz, 683, lambdamin, lambdamax, 1,
	   kjump) != 0)
	{
	  Error ("In one_kurucz after return from pdf_gen_from_array\n");
	}
//              pdf_sum(&pdf_kurucz,lambdamin,lambdamax,edges,nedges,w_kurucz,f_kurucz,683);
      old_t = t;
      old_g = g;
      old_freqmin = freqmin;
      old_freqmax = freqmax;
    }

  f = (C * 1.e8 / pdf_get_rand (&pdf_kurucz));
  if (f > freqmax)
    {
      Error ("one_kurucz: f too large %e\n");
      f = freqmax;
    }
  if (f < freqmin)
    {
      Error ("one_kurucz: f too small %e\n");
      f = freqmin;
    }
  return (f);
}

/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:
		kurucz(model,wavelength,flux)
		char	*model;   
		double	*wavelength; 
		double	*flux;		


  Description:	This routine reads kurucz models from the directory defined
	by KURUCZDIR.  Each model is stored in a separate file.  The files
	are listed in a file in this directory in the file kurucz_list.  The
	wavelengths are in the file "kurucz_wavelengths".  The model files
	and kurucz_wavelengths are binary files (which in this case have
	been created by kuruczwrite.

	This version of kurucz is for the new kurucz models that are on the
	tape sent to KSL in early 1991.  It cannot be used with the old models.

  Arguments:		

  Returns:
	The program returns wavelengths and fluxes.  The wavelengths
	are in angstroms, the fluxes are in ergs/cm**2/s/Angstrom.  The
	files themselves contain the fluxes in ergs/cm**2/s/Hz so a
	conversion has been made.  

  Notes:
	It may be that it would be more appropriate not to do the conversion
	described above, i.e. to leave the fluxes in ergs/cm**2/s/Hz.  However
	the earlier Kurucz tapes were in ergs/cm**2/s/Angstrom and so making
	the conversion provides the most continuity for now.

  History:
	2/4/91	ksl 	Modified from HCF's earlier subroutine of the same name
	8/8/97 	ksl	Modified to use fopen and fread (which are ANSI C routines).
	5Aug99	ksl	Removed possibility of running on MAC and modified too
			make perational on linux.  Note that the directory 
			data.91 needs to be linked symbolically in the directory 
			where the program is being run.   

 ************************************************************************/

#include <stdio.h>
#include <strings.h>
#include <string.h>
#include <math.h>
#include <fcntl.h>
#include "log.h"

#define WAVEFILE 	"data.91/kurucz_wavelengths"
#define KURUCZDIR	"data.91/"
#define	KURUCZLIST	"data.91/kurucz_index"

#define NWAVES		683	/*number of wavelengths in the Kurucz grid */
#define ARRAYSIZE 	8 * 683	/* Number of characters in each array */

#define NOK		-1
#define OK		0

#define NAMELEN		40
#define LINELEN     132
#define	DEBUG		0	/*0 if no debug, !0 if debug */

int kurucz_start = 0;

int
kurucz (model, wavelength, flux)
     char *model;
     double wavelength[];
     double flux[];
{
  FILE *fopen (), *fptr;
  char filename[80];

  /* Read in wavelength list */
  if ((fptr = fopen (WAVEFILE, "r")) == NULL)
    {
      perror ("kurucz:wavefile");
      return (NOK);
    }

//      if (fread(wavelength,1,ARRAYSIZE,fptr) != ARRAYSIZE )
  if (fread (wavelength, 8, 683, fptr) != 683)
    {
      Error ("kurucz: too few wavelength values\n");
      return (NOK);
    }
#if DEBUG
  printf ("wave_k %lf %lf %lf\n", wavelength[1], wavelength[100],
	  wavelength[500]);
#endif

  fclose (fptr);

  /* Read in flux list */
  strcpy (filename, KURUCZDIR);
  strcat (filename, model);

  if ((fptr = fopen (filename, "r")) == NULL)
    {
      perror ("kurucz:fluxfile");
      return (NOK);
    }

  if (fread (flux, 1, ARRAYSIZE, fptr) != ARRAYSIZE)
    {
      Error ("kurucz: too few flux values\n");
      return (NOK);
    }

#if DEBUG
  printf ("flux_k %f %f %f\n", flux[1], flux[100], flux[500]);
#endif

  fclose (fptr);

/* Convert from ergs/cm**2/Hz to ergs/cm**2/Angstrom */

  kurucz_change_units (wavelength, flux);

  return (OK);
}

int
kurucz_change_units (wavelength, flux)
     double *wavelength;
     double *flux;
{
  int i;

  for (i = 0; i < NWAVES; i++)
    flux[i] *= (C * 1e8 / (wavelength[i] * wavelength[i]));

  return (OK);

}

/*-************************************************************************
                    Space Telescope Science Institute
  Synopsis: 

			kurucz_get_model(t,g,wave,flux)
			double t,g,wave[],flux[];

			kurucz_which_model(t,g,no,frac)
			double t,g,abund;
			int		no[4];
			double	frac[4];

			kurucz_get_list()

  Description:	These subroutines find what Kurucz models best correspond
		to the required temperature, gravity, and abundance
		requested.  kurucz_get_model actually returns wavelengths and
		data for the best model.  kurucz_get_list needs to be called only once to
		get the list of all the models which exist.  For now
		the program returns a single model, by first finding
		the closest temperature to the one specified and then
		finding the gravity which is closest to that required.

		Therefore, the program does not really interpolate between Kurucz
		models.

  Arguments:		
			t	desired T in degrees
			g	log gravity (A change has been added to kurucz_which_model
				so that if g is large, then the assumption is that the
				g is being expressed directly, not as a log)
			no	the four models which define the region surrounding
				the desired model.  The models are sorted so that the
				model closest to the desired model will be first. At
				the present time only the closest model is returned.
			frac	the percentage contribution to apply from each
				model to the creation of the desired model
  Returns:

  Notes:
		The kurucz models are identified by a temperature and a log g.

  History:
	
	1/11/91	ksl	Began work
	2/1/91	ksl	This version revised to work with new 91 kurucz models
	3/2/96  	Modified so that when one cannot find a model of the appropriate
			temperature, then the model is scaled by the ratio of t/tmodel**4.
			Note that a better approach might be to scale the fluxes by
			the ratios one would predict from a Planck function, but this
			is not what is implemented here (at this time).
	8/8/97  	Modified so it was ANSI compatible.  Removed my line oriented versions
			of lscanf.  Made modifications so it would work on Mac.
	8/25/97 	Modified so successive calls to kurucz_get_model do not read from the disk unless
			a new model is needed.  This is a necessity if one hopes to use this program with
			Python.
	01sept	ksl	Modified so that the model is scaled by a Planck function when tbest is
			not exactly t
	02mar	ksl	Modified kurucz_which_model so that if g is large, then it is assumed
			that g has been passed in cgs units directly rather than log g.
	02may	ksl	Added a fudge in situations in which the Kurucz models have zero flux.
			Hopefull this will fix a problem which arose when I split the wavelength
			intervals in photon_gen.  The real fix is to have models which cover all
			wavelength ranges.
	02jun	ksl	Fixed potential problem that would occur if very low temperature models
			are called first (after finding the problem in my cv program).

 ************************************************************************/

char list_no[2000][NAMELEN];
float list_t[2000], list_g[2000];
int nlist, ktest = 999;
double wave[NWAVES], f[NWAVES];
char old_name[NAMELEN];


/* Get the kurucz model fluxes for a specific temperature and gravity */

int
kurucz_get_model (t, g, w, flux)
     double t, g, w[], flux[];
{
  char no[4][NAMELEN];
  double frac[4];
  double tbest, q, hnu;
  int n;
  double pow ();
  double kurucz_which_model ();
  double hnu_zero;
  int nzero;

  if (ktest == 999)
    {
      kurucz_get_list ();
      ktest = 84;
    }
  tbest = kurucz_which_model (t, g, no, frac);

  if (strcmp (old_name, no[0]) != 0)
    {
#if DEBUG
      Log ("Getting model %s\n", no[0]);
#endif
      if (kurucz (no[0], wave, f) == NOK)
	{
	  Error ("Could not read model %s\n", no[0]);
	  exit (0);
	}
      strcpy (old_name, no[0]);
    }

/* Next section corrects for the problem of fluxes being zero at certain
short wavelengths in certain cool star Kuurcz spectra.  It is a fudge!!! */

  nzero = 0;
  while (f[nzero] == 0.0)
    nzero++;
  for (n = 0; n < nzero; n++)
    {
      if (n == 0)
	{			// initialize
	  hnu_zero = HC / (wave[nzero] * ANGSTROM);
	}
      hnu = HC / (wave[n] * ANGSTROM);
      f[n] =
	f[nzero] * pow ((hnu / hnu_zero),
			5.0) * exp ((hnu_zero - hnu) / (BOLTZMANN * tbest));

    }


//  q below corrects for the possibility that one may not have a model
//  with the appropriate temperature. 
  for (n = 0; n < NWAVES; n++)
    {
      w[n] = wave[n];
      hnu = HC / (w[n] * ANGSTROM);
      q =
	(exp (hnu / (BOLTZMANN * tbest)) -
	 1.0) / (exp (hnu / (BOLTZMANN * t)) - 1.0);
      flux[n] = q * f[n];
    }
  return (0);
}


/*get the list of models */
int
kurucz_get_list ()
{
  int i;
  char line[LINELEN];
  char *fgets ();

  FILE *fptr, *fopen ();

  if ((fptr = fopen (KURUCZLIST, "r")) == NULL)
    {
      Error ("Could not open %s\n", KURUCZLIST);
      exit (-1);
    }

  i = 0;
  while ((fgets (line, LINELEN, fptr)) != NULL)
    {
      sscanf (line, "%s %f %f", list_no[i], &list_t[i], &list_g[i]);
      list_g[i] *= (0.1);	/*model list in g in tenths of log */
      i++;
    }

  nlist = i;
#if	DEBUG
  for (i = 0; i < 100; i++)
    Log ("name t g %s %g %g\n", list_no[i], list_t[i], list_g[i]);

#endif
  return (0);
}


double
kurucz_which_model (t, g, no, frac)
     double t, g;
     char no[4][NAMELEN];
     double frac[4];
{
  int i;
  double dt, dg;
  double x, xtest;
  double tbest, gbest;

/* 02mar.  The Kurucz grid is in log g and  in general
it should be called as such.  However, if g is greater 
than 100 we will assume that g has been passed in cgs units
*/

  if (g > 100)
    g = log (g);

//Log("find %g %g\n",t,g);

/* find the model that is closest to the one needed */
#if DEBUG
  Log ("No of models %d\n", nlist);
#endif
  x = 1.e20;
  for (i = 0; i < nlist; i++)
    {
      if (list_t[i] != 0.0)
	{
	  dt = list_t[i] - t;
	  xtest = dt * dt;
	  if (x > xtest)
	    {
	      x = xtest;
	      tbest = list_t[i];
	    }
	}
    }
#if DEBUG
  Log ("Best temperature %g dt %g\n", tbest, dt);
#endif
  x = 1.e50;

  for (i = 0; i < nlist; i++)
    {
      dt = list_t[i] - tbest;
      if (dt * dt < 1.e-3)
	{
	  dg = list_g[i] - g;
	  xtest = dg * dg;
	  if (x > xtest)
	    {
	      x = xtest;
	      gbest = list_g[i];
	      strcpy (no[0], list_no[i]);
#if DEBUG
	      Log ("Copying %s to no %s\n", list_no[i], no[0]);
#endif
	    }
	}
    }


  frac[1] = frac[2] = frac[3] = 0.0;
  frac[0] = 1.0;
  for (i = 1; i < 4; i++)
    strcpy (no[i], "");
#if DEBUG
  Log ("debug listname %s\n", list_no[0]);
  Log ("tbest %f gbest %f n %s\n", tbest, gbest, no[0]);
#endif

  return (tbest);

}


// Adapted from emittance_hub    ksl 01aug
int ierr_emit = 0;

double
emittance_kurucz (freqmin, freqmax, t, g)
     double freqmin, freqmax, t, g;
{
  int kurucz_get_model (), nwav, n;
  double w, x, lambdamin, lambdamax;
  double dlambda;


  lambdamin = C / (freqmax * ANGSTROM);
  lambdamax = C / (freqmin * ANGSTROM);
  kurucz_get_model (t, g, w_kurucz, f_kurucz);
  nwav = 683;
  if (lambdamax > w_kurucz[nwav - 1])
    {
      if (ierr_emit == 0)
	{
	  Error ("emittance_kurucz: Wavelengths out of range of models\n");
	  Error ("lambda %f %f  model %f %f\n", lambdamin, lambdamax,
		 w_kurucz[0], w_kurucz[nwav - 1]);
	  ierr_emit = 1;
	}
      lambdamax = w_kurucz[nwav - 1];
    }
  if (lambdamin < w_kurucz[0])
    {
      if (ierr_emit == 0)
	{
	  Error ("emittance_kurucz: Wavelengths out of range of models\n");
	  Error ("lambda %f %f  model %f %f\n", lambdamin, lambdamax,
		 w_kurucz[0], w_kurucz[nwav - 1]);
	  ierr_emit = 1;
	}
      lambdamin = w_kurucz[0];
    }
  x = 0;
  for (n = 0; n < nwav; n++)
    {
      w = w_kurucz[n];
      if (n == 0)
	{
	  dlambda = w_kurucz[1] - w_kurucz[0];
	}
      else if (n == nwav - 1)
	{
	  dlambda = w_kurucz[n] - w_kurucz[n - 1];
	}
      else
	{
	  dlambda = 0.5 * (w_kurucz[n + 1] - w_kurucz[n - 1]);
	}

      if (lambdamin < w && w < lambdamax)
	{
	  x += f_kurucz[n] * dlambda;
	}
    }
  x *= 4. * PI;
  return (x);
}
