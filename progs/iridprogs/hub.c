#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "atomic.h"
#include "python.h"

/*  The structure is defined in python.h.  Here for reference only */
//#define NPDF 200

//typedef struct Pdf {
//      double x[NPDF+1];
//      double y[NPDF+1];
//      double limit1,limit2;
//} *PdfPtr,pdf_dummy;



/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:

double one_hub(t,g,freqmin,freqmax) gets a photon frequency from a hubeny grid 

Arguments:

	double t,g,				temperature and gravity of star in Kurucz grid
	freqmin,freqmax;			minimum and maximum frequency of interest


Returns:

 
Description:	
	

Notes:

History:
	97aug	ksl 	coded as part of python project.
	01jul	ksl	Adapted from corresponding kurucz routine
 
**************************************************************/



double old_t, old_g, old_freqmin, old_freqmax;
double w_hub[4096], f_hub[4096];
struct Pdf pdf_hub;
double jump[] = { 913.8 };

double
one_hub (t, g, freqmin, freqmax)
     double t, g, freqmin, freqmax;
{
  double lambdamin, lambdamax;
  double f;
  double pdf_get_rand ();
  int get_hub (), nwav;

  if (old_t != t || old_g != g || old_freqmin != freqmin
      || old_freqmax != freqmax)
    {				/* Then we must initialize */
      nwav = get_hub (t, g, w_hub, f_hub);
      /*  Get_model returns wavelengths in Ang and flux in ergs/cm**2/Ang */
      lambdamin = C * 1e8 / freqmax;
      lambdamax = C * 1e8 / freqmin;
      if (pdf_gen_from_array
	  (&pdf_hub, w_hub, f_hub, nwav, lambdamin, lambdamax, 1, jump) != 0)
	{
	  Error ("In one_hub after return from pdf_gen_from_array\n");
	}
//              pdf_sum(&pdf_hub,lambdamin,lambdamax,edges,nedges,w_hub,f_hub,nwav);
      old_t = t;
      old_g = g;
      old_freqmin = freqmin;
      old_freqmax = freqmax;
    }

//      f=(C*1.e8/rand_pdf(&pdf_hub));
  f = (C * 1.e8 / pdf_get_rand (&pdf_hub));
  if (f > freqmax)
    {
      Error ("one_hub: f too large %e\n");
      f = freqmax;
    }
  if (f < freqmin)
    {
      Error ("one_hub: f too small %e\n");
      f = freqmin;
    }
  return (f);
}



/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:

  Description:	These subroutines get a interpolate between a
  series of hubeny's models in t and usually log g.  

  Arguments:		

	int get_hub(t,g,w,flux,files)

		double 	t			temperature
		double	g			log g
		double	w[],flux[]	wavelengths and fluxes to be returned
		char	files[]     full name of file containing list of spectra


  Returns:

	get_hub--returns the number of wavelengths in the spectrum

  Notes:

  The program interpolates between various models in two dimensions
  (nominally T and g).  It takes a temperature and a gravity
  and searches all the models for the set which has T just below
  the desired value and just above the desired value.  The grid need
  not be equally spaced, but the program first interpolates
  between the set of models with T just below the desired value of
  T and then interpolates between the set of models with T just
  above the desired value.  It then linearly interpolates the result.

  The logic of the program is such that you need to have at
  least one model at the low and hi T edges of the grid.  If you
  are off the grid on the hi side, you interpolate between
  available models at the max T.  Similarly if you are off on
  the low side you interpolate between the models at the min. T.


  An example of part of a file is as follows:

  5000. 115000.                      	!Min and max value of T 
  7 9									!Min and max value of G 
  /pisces/data1/long/hubeny/db/zero.0920.11 10000 8.
  /pisces/data1/long/hubeny/db/test/hhe158lt.0920.11 15000 8.
  /pisces/data1/long/hubeny/db/test/hhe158lt.0920.11 15000 6.
  /pisces/data1/long/hubeny/db/test/hhe258lt.0920.11 25000 8.
  /pisces/data1/long/hubeny/db/db.final/hhe358lt.0920.11 35000 8.

  If NWAVES < the number of wavelengths in the files then only
  the first nwave wavelengths will be returned.
  
  Note that in this example if you asked for t=12500 g=7 model.  It
  would weight the 10000 8 model by 0.5 and the two 15000 degree models
  by 0.25

  History:

  10/92   ksl  Coded
  12/92   ksl  Updated to read all models when first attempt to
			   read a model is made
  4/93   ksl  Updated to so that the grid need not be equally spaced
               in either t or g
 99aug5	ksl   Updated for linux.  The main change is to remove lscanf and
	      to alter the way filenaming is handled.
 00feb20 ksl	Increased maximum number of wavelengths to 4096.  Also cleaned
		out all -Wall cautions.  Note that this routine has more
		external variables than I would program it with today.
 01jul	ksl	Modified for use in Python
 ************************************************************************/


#include <stdio.h>
#include <math.h>

#define LINELEN	132
#define NMODELS 60
#define NWAVES 4096
#define DEBUG  0
#define EPST   0.1
#define  EPSG   .01

int igethublist = -1;

int nmod, nwav;
double hub_tmin, hub_tmax, hub_tstep, waves[NWAVES];
double hub_gmin, hub_gmax, hub_gstep;
struct model_desc
{
  char name[80];
  double t, g;
  double flux[NWAVES];
}
model[NMODELS];


int
get_hub (t, g, w, flux)
     double t, g, w[], flux[];

{
  double fluxtlo[NWAVES], fluxthi[NWAVES];
  extern double hub_tmin, hub_tmax, waves[NWAVES];
  extern int igethublist, nmod, nwav;
  double dt, thi, tlo, dthi, dtlo, wtlo, wthi;
  double dg, ghi, glo, dghi, dglo;
  double whi, wlo;
  int i, iglo, ighi;
  int get_allmodels ();

//printf("find t, g %f %f\n",t,g);

/* If this is the first time the subroutine is called, it reads
all the models */

  if (igethublist < 0)
    {
      get_allmodels (hubeny_list);
      printf ("hub Read all input models \n");
      printf ("hub t %f %f\n", hub_tmin, hub_tmax);
      printf ("hub g %f %f\n", hub_gmin, hub_gmax);
    }

  igethublist = 1;

/* Force t and g to be within the grid */

  thi = tlo = 0;
  if (t > hub_tmax)
    {
      t = hub_tmax;
    }
  if (t < hub_tmin)
    {
      t = hub_tmin;
    }
  if (g < hub_gmin)
    g = hub_gmin;
  if (g > hub_gmax)
    g = hub_gmax;


/* Now you must find which models to use and
with what weights */

/* First establish the boundaries in t */
  thi = dthi = tlo = dtlo = 1e33;
  for (i = 0; i < nmod; i++)
    {
      dt = model[i].t - t;
#if DEBUG
      printf ("dt mod.t t %.2g %.2g %.2g\n", dt, model[i].t, t);
#endif
      if (dt >= 0 && dt < dthi)
	{
	  thi = model[i].t;
	  dthi = dt;
	}
      if (dt < 0 && (-dt) < dtlo)
	{
	  tlo = model[i].t;
	  dtlo = (-dt);
	}
    }
  wtlo = (thi - t) / (thi - tlo);
  wthi = (1. - wtlo);

#if DEBUG
  printf ("thi tlo %.2g %.2g %.2g %.2g\n", thi, tlo, wtlo, wthi);
#endif

/* So now we have defined thi and tlo.  We now
need to get the weights for different gravities 
at tlow and thi*/

/* First work on the low T end */

  ghi = dghi = glo = dglo = 1e33;
  iglo = ighi = 0;

  for (i = 0; i < nmod; i++)
    {
      if (fabs (model[i].t - tlo) < EPST)
	{

	  dg = model[i].g - g;
#if DEBUG
	  printf ("match tlo: %d %f %f %f\n", i, model[i].t, model[i].g, dg);
#endif
	  if (dg >= 0 && dg < dghi)
	    {
	      ghi = model[i].g;
	      dghi = dg;
	      ighi = i;
	    }
	  if (dg < 0 && (-dg) < dglo)
	    {
	      glo = model[i].g;
	      dglo = (-dg);
	      iglo = i;
	    }
	}
    }

#if DEBUG
  printf ("low T glo ghi iglo ighi %.2lg %.2lg %d %d %.2lg %.2lg\n", glo, ghi,
	  iglo, ighi, dghi, dglo);
#endif

  if (ghi == glo)
    {
      for (i = 0; i < nwav; i++)
	{
	  fluxtlo[i] = model[ighi].flux[i];
	}
    }
  else
    {
      wlo = (ghi - g) / (ghi - glo);
      whi = (1. - wlo);
      for (i = 0; i < nwav; i++)
	{
	  fluxtlo[i] = whi * model[ighi].flux[i] + wlo * model[iglo].flux[i];
	}
#if DEBUG
      printf ("weight glo: ig lo-hi w lo-hi %d %d %f %f\n", iglo, ighi, wlo,
	      whi);
#endif

    }

/* Now do the hi T end*/
  ghi = dghi = glo = dglo = 1e33;
  iglo = ighi = 0;
  for (i = 0; i < nmod; i++)
    {
      if (fabs (model[i].t - thi) < EPST)
	{
	  dg = model[i].g - g;

#if DEBUG
	  printf ("match thi: %d %f %f %f\n", i, model[i].t, model[i].g, dg);
#endif
	  if (dg >= 0 && dg < dghi)
	    {
	      ghi = model[i].g;
	      dghi = dg;
	      ighi = i;
	    }
	  if (dg < 0 && (-dg) < dglo)
	    {
	      glo = model[i].g;
	      dglo = (-dg);
	      iglo = i;
	    }
	}
    }
#if DEBUG
  printf ("hi T glo ghi iglo ighi %.2lg %.2lf %d %d %.2lg %.2lg\n", glo, ghi,
	  iglo, ighi, dghi, dglo);
#endif
  if (ghi == glo)
    {
      for (i = 0; i < nwav; i++)
	{
	  fluxthi[i] = model[ighi].flux[i];
	}
    }
  else
    {
      wlo = (ghi - g) / (ghi - glo);
      whi = (1. - wlo);
      for (i = 0; i < nwav; i++)
	{
	  fluxthi[i] = whi * model[ighi].flux[i] + wlo * model[iglo].flux[i];
	}
#if DEBUG
      printf ("weight hi T ig lo-hi w lo-hi %d %d %f %f\n", iglo, ighi, wlo,
	      whi);
#endif

    }

/* Now combine the high and the low temperature */

  for (i = 0; i < nwav; i++)
    {
      w[i] = waves[i];		/*note waves is external */
      flux[i] = wthi * fluxthi[i] + wtlo * fluxtlo[i];
    }


  return (nwav);

}


/* This routine reads the "filelist" which contains all the
models.  It is called only one time */

int
get_allmodels (filename)
     char filename[];
{
  FILE *iptr, *jptr, *fopen ();
  int kkk, nw;
  double g, t, f[NWAVES], w[NWAVES];
  char name[80];
  char fname[LINELEN];
  char line[LINELEN];
  extern double hub_tmin, hub_tmax, waves[NWAVES];
  extern int nmod, nwav;


/* Create the complete name to the file */
  printf ("Index of models being used in %s  \n", filename);

/* Now read in the models which exist */

  if ((iptr = fopen (filename, "r")) == NULL)
    {
      fprintf (stderr, "Error opening file %s\n", filename);
      exit (1);
    }

  if ((fgets (line, LINELEN, iptr)) != NULL)
    {
      sscanf (line, "%lf %lf", &hub_tmin, &hub_tmax);
    }
  else
    {
      fprintf (stderr, "help: Unexpected EOF\n");
      exit (0);
    }


  if ((fgets (line, LINELEN, iptr)) != NULL)
    {
      sscanf (line, "%lf %lf", &hub_gmin, &hub_gmax);
    }
  else
    {
      fprintf (stderr, "help: Unexpected EOF\n");
      exit (0);
    }


  nmod = 0;
  while ((fgets (line, LINELEN, iptr)) != NULL)
    {
      sscanf (line, "%s %lf %lf", fname, &t, &g);
      sprintf (name, "%s", fname);

      if (t < hub_tmin || t > hub_tmax || g < hub_gmin || g > hub_gmax)
	{
	  fprintf (stderr, "Model %s %g %g out of range and ignored\n", fname,
		   t, g);
	}
      else
	{
	  model[nmod].t = t;
	  model[nmod].g = g;
	  strcpy (model[nmod].name, name);


	  printf ("Reading model %d %s %f %f\n",
		  nmod, model[nmod].name, model[nmod].t, g);

/* now read the model */


	  if ((jptr = fopen (name, "r")) == NULL)
	    {
	      fprintf (stderr, "Error trying to open %s\n", model[nmod].name);
	      exit (1);
	    }

	  else
	    {			/* read the file */
	      nw = 0;

	      /* if nw=NWAVES fgets will not be executed */

	      while (nw < NWAVES && fgets (line, LINELEN, jptr) != NULL)
		{
		  sscanf (line, "%lf %lf", &w[nw], &f[nw]);
		  nw++;
#if DEBUG
		  printf ("ret from fgets %lf %lf\n", w[nw], f[nw]);
#endif
		}

	      fclose (jptr);
	    }


	  if (nmod == 0)
	    {
	      nwav = nw;
	      for (kkk = 0; kkk < nwav; kkk++)
		waves[kkk] = w[kkk];
	    }
	  if (nmod > 0 && nwav != nw)
	    {
	      fprintf (stderr,
		       "Error: File %s did had %d waves instead of %d\n",
		       name, nw, nwav);
	      exit (1);
	    }
	  for (kkk = 0; kkk < nwav; kkk++)
	    model[nmod].flux[kkk] = f[kkk];
#if DEBUG
	  for (kkk = 0; kkk < nwav; kkk++)
	    printf ("%lf %lf\n", w[kkk], f[kkk]);
#endif
	  nmod++;
	}

    }

  fclose (iptr);

  return (0);
}

double
emittance_hub (freqmin, freqmax, t, g)
     double freqmin, freqmax, t, g;
{
  int get_hub (), nwav, n;
  double w, x, lambdamin, lambdamax;
  double dlambda;

  lambdamin = C / (freqmax * ANGSTROM);
  lambdamax = C / (freqmin * ANGSTROM);
  nwav = get_hub (t, g, w_hub, f_hub);
  if (lambdamax > w_hub[nwav - 1] || lambdamin < w_hub[0])
    {
      Error ("emittance_hub: Wavelengths out of range of models\n");
      Error ("lambda %f %f  model %f %f\n", lambdamin, lambdamax, w_hub[0],
	     w_hub[nwav - 1]);
      exit (0);
    }
  x = 0;
  for (n = 0; n < nwav; n++)
    {
      w = w_hub[n];
      if (n == 0)
	{
	  dlambda = w_hub[1] - w_hub[0];
	}
      else if (n == nwav - 1)
	{
	  dlambda = w_hub[n] - w_hub[n - 1];
	}
      else
	{
	  dlambda = 0.5 * (w_hub[n + 1] - w_hub[n - 1]);
	}
      if (lambdamin < w && w < lambdamax)
	{
	  x += f_hub[n] * dlambda;
	}
    }
  x *= 4. * PI;
  return (x);
}
