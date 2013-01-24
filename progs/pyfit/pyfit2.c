


/**************************************************************************
                    Center for Astrophysical Sciences
                        Johns Hopkins University


  Synopsis:  

  Description:	

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:

  This is roughly what the Makefile should look like. Uncomment the mv 
  statement to store binary in ~/bin

CC = gcc
CFLAGS = -c -g -I$$HOME/include -Wall
LDFLAGS= -L$$HOME/lib -lm -lkpar
BIN = $$HOME/bin

pyfit: pyfit.o
	gcc pyfit.o  $(LDFLAGS) -o pyfit
#	mv $@ $(BIN)

 

  History:
2003	ksl	Adapted from kslfit


 ************************************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


#include "log.h"



#define LINELENGTH 132
#define NWAVES      8500
#define MAX_MOD    10
#define BIG        1.e32
#define MAXTOT     10000
#define NREG	10

struct spect
{
  char name[LINELENGTH];
  float ang;
  float w[NWAVES], w_edge[NWAVES];
  float f[NWAVES];
  float e[NWAVES];
  int good[NWAVES];
  int nw, ng;
}
spec, mod_dummy[MAX_MOD], mod[MAX_MOD], mod_best;


struct res
{
  char name[LINELENGTH];
  float ang;
  float norm;
  float chinu;
}
results[MAXTOT];

int nres;

float xx_bestnorm, xx_normchi;
float normchi_best;

#include "templates.h"

int
main (argc, argv)
     int argc;
     char *argv[];
{

  FILE *fptr, *optr, *fopen ();
  char parfile[LINELENGTH];
  char source[LINELENGTH];
  char outname[LINELENGTH];
  char outputfile[LINELENGTH];
  char output_specfile[LINELENGTH];
  char model_list[LINELENGTH];
  char line[LINELENGTH];
  char modelname[LINELENGTH];
  char waverange[NREG][LINELENGTH];
  int i, m, mspec;
  int n, nreg;
  float x;
  float incmin, incmax;
  int nmods;

/* Parse the command line */
  if (argc == 1)
    strcpy (parfile, "pyfit.pf");
  else if (argc == 2)
    {
      strcpy (parfile, argv[1]);
      if ((strstr (parfile, ".pf")) == NULL)
	strcat (parfile, ".pf");
    }
  else
    {
      printf ("Usage: pyfit [file.pf]\n");
      exit (0);
    }

/* Initialize variables that will be read by rdpar */

  strcpy (source, "sscyg_ghrs1");
  strcpy (model_list, "kwd.ls");
  for (n = 0; n < NREG; n++)
    {
      strcpy (waverange[n], "1525-1575");
    }
  strcpy (outputfile, "pyfit.out");

  normchi_best = BIG;
  nres = 0;
  incmin = 0.0;
  incmax = 90.;

/* Get data from a parameter file */


  opar (parfile);


  rdstr ("source_file", source);
  rdstr ("model_list", model_list);
  rdflo ("minimum.inclination.angle", &incmin);
  rdflo ("maximum.inclination.angle", &incmax);
  rdint ("number.of.independent.regions", &nreg);
  for (n = 0; n < nreg; n++)
    {
      rdstr ("wavelengths.to.fit", waverange[n]);
    }
  rdstr ("output_file_root", outname);

  strcpy (outputfile, outname);
  strcat (outputfile, ".results");
  strcpy (output_specfile, outname);
  strcat (output_specfile, ".out");

  cpar (parfile);

/* End of input section */

  if ((optr = fopen (outputfile, "w")) == NULL)
    {
      printf ("Failed to open file %s\n", outputfile);
      exit (0);
    }



  /* Get the data */

  get_data (source);
  /* Open the input file.  Exit if it is not opened successfully */
  if ((fptr = fopen (model_list, "r")) == NULL)
    {
      printf ("Failed to open file %s\n", model_list);
      exit (0);
    }



/* Read and print the input file */
  nmods = 0;
  while (fgets (line, LINELENGTH, fptr) != NULL)
    {
      sscanf (line, "%s", modelname);
      mspec = get_python_model (modelname);
      if ((nmods % 100) == 0)
	printf ("%s\n", modelname);
      nmods++;

/* Now we have the data and the models and we can begin a various comparisons
*/

      for (m = 0; m < mspec; m++)
	{
	  if (incmin <= mod[m].ang && mod[m].ang <= incmax)
	    {
	      strcpy (results[nres].name, modelname);
	      results[nres].ang = mod[m].ang;
	      x = 0;
	      for (n = 0; n < nreg; n++)
		{
		  parse_waverange (waverange[n]);
		  chisq_pure (&mod[m]);
		  results[nres].norm = xx_bestnorm;
		  results[nres].chinu = xx_normchi / spec.ng;
		  x += results[nres].chinu;
		}
	      if (x < normchi_best)
		{
		  strcpy (mod_best.name, mod[m].name);
		  mod_best.ang = mod[m].ang;
		  normchi_best = x;
		  for (i = 0; i < spec.nw; i++)
		    {
		      mod_best.w[i] = mod[m].w[i];
		      mod_best.f[i] = xx_bestnorm * mod[m].f[i];

		    }
		}
	      nres++;

	    }

	}
    }

  /* OK now write out the best fitting model and sum up */
  for (i = 0; i < nres; i++)
    {
      fprintf (optr, "%30s %5.1f %10.3e %10.3e  %10.2f\n", results[i].name,
	       results[i].ang, results[i].norm, results[i].chinu,
	       results[i].chinu / (normchi_best / spec.ng));
    }

  fclose (optr);

  if ((optr = fopen (output_specfile, "w")) == NULL)
    {
      printf ("Failed to open file %s\n", output_specfile);
      exit (0);
    }

  fprintf (optr, "# %s %.1f %10.3e\n", mod_best.name, mod_best.ang,
	   normchi_best / spec.ng);
  for (i = 0; i < spec.nw; i++)
    {
      fprintf (optr, "%10.3f %10.3e %10.3e %10.3e %2d\n", spec.w[i],
	       spec.f[i], spec.e[i], mod_best.f[i], spec.good[i]);

    }



  printf ("# %s %.1f %10.3e\n", mod_best.name, mod_best.ang,
	  normchi_best / spec.ng);
/* OK now end program and go on to the next important task. */
  exit (0);

}


int
get_data (filename)
     char filename[];
{
  FILE *ptr;
  char dummy[LINELENGTH];
  int n;

  if ((ptr = fopen (filename, "r")) == NULL)
    {
      Error ("Could not open filename %s\n", filename);
      exit (0);
    }
  n = 0;
  while (n < NWAVES && (fgets (dummy, LINELENGTH, ptr)) != NULL)
    {
      if ((dummy[0] != '#'))
	{
	  sscanf (dummy, "%e %e %e", &spec.w[n], &spec.f[n], &spec.e[n]);
	  n++;
	}
    }
  spec.nw = n;

/* We assume that the wavelenghts in the grid represent the midpoints of
   the data bins. Now we need to define the edges */


  spec.w_edge[0] = spec.w[0] - 0.5 * (spec.w[1] - spec.w[0]);
  for (n = 1; n < spec.nw; n++)
    {
      spec.w_edge[n] = 0.5 * (spec.w[n - 1] + spec.w[n]);
    }
  spec.w_edge[spec.nw] =
    spec.w_edge[spec.nw - 1] + 0.5 * (spec.w_edge[spec.nw - 2] +
				      spec.w_edge[spec.nw - 1]);

  fclose (ptr);
}


int
get_python_model (filename)
     char filename[];
{
  FILE *ptr;
  char dummy[LINELENGTH];
  int n, m, nscanned;
  float ang[MAX_MOD], wavelength, f[MAX_MOD];
  int mspec, nang, nwaves;

  if ((ptr = fopen (filename, "r")) == NULL)
    {
      Error ("Could not open filename %s\n", filename);
      exit (0);
    }


// Get the title line

  for (m = 0; m < 6; m++)
    ang[m] = 999;
  mspec =
    fscanf (ptr,
	    "#Freq.        Lambda  Emitted  Star+BL  Disk     Wind     Scattered  Absorbed  Ang:%f  Ang:%f  Ang:%f  Ang:%f  Ang:%f Ang:%f",
	    &ang[0], &ang[1], &ang[2], &ang[3], &ang[4], &ang[5]);

  n = 0;
  while (n < NWAVES && (fgets (dummy, LINELENGTH, ptr)) != NULL)
    {
      if ((dummy[0] != '#'))
	{
	  nscanned =
	    sscanf (dummy, "%*e %f %*e %*e %*e %*e %*e %*e %e %e %e %e %e %e",
		    &wavelength, &f[0], &f[1], &f[2], &f[3], &f[4], &f[5]);
	  if (nscanned != mspec + 1)
	    Error ("Problem in file %s with line %s\n", filename, dummy);
	  for (m = 0; m < mspec; m++)
	    {
	      mod[m].w[n] = wavelength;
	      mod[m].f[n] = f[m];
	    }
	  n++;
	}
    }

  nwaves = n;
  for (m = 0; m < mspec; m++)
    {
      mod[m].ang = ang[m];
      strcpy (mod[m].name, filename);
      mod_dummy[m].nw = mod[m].nw = nwaves;
    }

/* OK, now invert the order of the indices.*/

  for (m = 0; m < mspec; m++)
    {
      for (n = 0; n < nwaves; n++)
	{
	  mod_dummy[m].w[n] = mod[m].w[nwaves - 1 - n];
	  mod_dummy[m].f[n] = mod[m].f[nwaves - 1 - n];
	}

      mod_dummy[m].w_edge[0] =
	mod_dummy[m].w[0] - 0.5 * (mod_dummy[m].w[1] - mod_dummy[m].w[0]);
      for (n = 1; n < nwaves; n++)
	{
	  mod_dummy[m].w_edge[n] =
	    0.5 * (mod_dummy[m].w[n - 1] + mod_dummy[m].w[n]);
	}

      mod_dummy[m].w_edge[nwaves] =
	mod_dummy[m].w_edge[nwaves - 1] +
	0.5 * (mod_dummy[m].w_edge[nwaves - 2] +
	       mod_dummy[m].w_edge[nwaves - 1]);

      /* OK, now put them on the same wavelength scale as the spectrum */

      regrid (&mod_dummy[m], &mod[m]);

    }


  fclose (ptr);

  return (mspec);
}



/* Regridding is non-trivial and depends on assumption about the input models.  We assume that
   the fluxes in the input models represent an average value within a wavelength interval. We
   also assume that wavelengths increase through the arrays
 */

regrid (modin, modout)
     struct spect *modin;	// input model with any old wavelength scale

     struct spect *modout;	// regridded model on the scale of the data

{
  int i, imin, imax, n;
  float x, y, wbot, wtop;

  imin = 0;
  for (n = 0; n < spec.nw; n++)
    {
      modout->f[n] = 0;
// Identify by imin the first input interval which contributes to the flux  
      while (modin->w_edge[imin] <= spec.w_edge[n] && imin < modin->nw)
	imin++;
      imin--;
// Identify by imax the last input bin which contributes to this output interval
      imax = imin;
      while (modin->w_edge[imax] < spec.w_edge[n + 1] && imax < modin->nw)
	imax++;
      imax--;
// Now calculate the average flux in the output wavelength enterval
      x = y = 0;
//      Log("regrid imin %d imax %d\n",imin,imax);
      for (i = imin; i <= imax; i++)
	{
	  wbot = modin->w_edge[i];
	  if (wbot < spec.w_edge[n])
	    wbot = spec.w_edge[n];
	  wtop = modin->w_edge[i + 1];
	  if (wtop > spec.w_edge[n + 1])
	    wtop = spec.w_edge[n + 1];
	  x += (modin->f[i]) * (wtop - wbot);
	  y += (wtop - wbot);
	}
      if (y > 0)
	modout->f[n] = x / y;
    }
// This next is for debugging: Check average fluxes before/after regridding
  x = y = 0;
  for (i = 0; i < modin->nw; i++)
    x += modin->f[i];
  for (i = 0; i < spec.nw; i++)
    y += modout->f[i];
//  Log ("regrid: in %e out %e\n", x / modin->nw, y / spec.nw);

}


float
chisq_pure (mod)
     struct spect *mod;
{
  double xx, xz, zz, ee, eee;	// These are double to increase exponent range

  float c;
  int n;

  xx = xz = zz = 0.0;

  for (n = 0; n < spec.nw; n++)
    {
      if (spec.e[n] > 0.0 && spec.good[n])
	{
	  eee = spec.e[n];
	  ee = 1. / (eee * eee);
	  xx += spec.f[n] * spec.f[n] * ee;
	  xz += spec.f[n] * mod->f[n] * ee;
	  zz += mod->f[n] * mod->f[n] * ee;
	}
    }
  xx_bestnorm = xz / zz;
  xx_normchi =
    xx - 2. * (xx_bestnorm) * xz + (xx_bestnorm) * (xx_bestnorm) * zz;
  c = xx - 2. * xz + zz;
  if (c > BIG)
    c = BIG;
//  Log ("chisq:xx %8.3g xz %8.3g zz %8.3g bestnorm %8.3g normchi %8.3g chisq %8.3g \n",
  //       xx, xz, zz, xx_bestnorm, xx_normchi, c);
  return (c);

}

int
parse_waverange (waverange)
     char waverange[];
{
  int i, n;
  int nintervals;
  float www, w1[10], w2[10];
  char *xxx;

/* parse the string */
  xxx = waverange;
  nintervals = 0;
  while (xxx != NULL && nintervals < 10)
    {
      sscanf (xxx, "%f-%f", &w1[nintervals], &w2[nintervals]);
//      Log ("found interval n %d w1,w2 %f %f\n", nintervals + 1,
//         w1[nintervals], w2[nintervals]);
      nintervals++;
      if ((xxx = index (xxx, ',')) == NULL)
	break;
      xxx++;
    }


  for (n = 0; n < spec.nw; n++)
    spec.good[n] = 0;

  spec.ng = 0;
  for (n = 0; n < spec.nw; n++)
    {
      www = spec.w[n];
      for (i = 0; i < nintervals; i++)
	{
	  if (w1[i] <= www && www <= w2[i])
	    {
	      spec.good[n] = 1;
	      spec.ng++;
	    }
	}
    }

//  Log ("parse_wavelength_range: %d datapoints in %d wavelength intervals\n",
//       spec.ng, nintervals);
  if (spec.ng < 1)
    {
      Error ("parse_wavelength_range:  No data points to fit\n");
      exit (0);
    }

  return (0);
}
