


/**************************************************************************
                        Space Telescope Science Institute


  Synopsis:  py_compare 

  Description:	

  	This routine compare two python models with one another 

	Usage 

	py_compare [-m firstmodel nspec -m secondmodel nspec] pffile  

	Old

  	This is a routine to fit a spectrum to a series of python models.
	The intent is to fit a series of regions (lines) rather than
	the ovrall spectral shape.  An overall spectral shape is also
	fit from the sum of all the accepted regions, but this does not contribute to 
	the "best spectrum".  pyfit3 treats each inclination angle in the python 
	output files separately; it does not interpolate to produce a best fit for
	and intermediate inclination.

  Arguments:  (Input via .pf file)		

  The parameter file looks something like this 

first_model                             sscyg_fuse4    !spectrum to fit
first_spectrum
second_model
second_spectrum
number.of.independent.regions                    4     !number of regions to fit
wavelengths.to.fit                         920-960     
chi2.renorm                                       1    !Weight to give this region
wavelengths.to.fit                        960-1000
chi2.renorm                                       1
wavelengths.to.fit             1000-1025,1026-1050
chi2.renorm                                       1
wavelengths.to.fit                       1140-1185
chi2.renorm                                       1
output_file_root                             sscyg_fuse



  Returns:

The beginning of the output file (.results) looks something like this:

#  ratio  ave_chinu    ratio_ind...                      ratio_comb ave_norm angle params                                                                                          filename
   1.500  1.29e+03      2.089    1.983    3.411    1.754    1.945  1.15e+00  15.0  5.00e-10  5.00e-11  4.00e+00  0.00e+00  1.00e+00
 1.00e+00  1.00e+10 sscyg_kgrid/sscyg_00_00_00_00_00_00_00.spec


The program also produces the spectra in the file output_file_root.out.  This contaians the 
best overall fit as well as the individual spectra.

# sscyg_kgrid/sscyg_04_03_01_02_00_02_02.spec 30.0  3.439e+03
# sscyg_kgrid/sscyg_04_02_05_00_02_00_02.spec 75.0  7.062e+02 920-960
# sscyg_kgrid/sscyg_02_01_03_00_01_00_02.spec 75.0  4.465e+02 960-1000
# sscyg_kgrid/sscyg_04_03_01_02_00_02_02.spec 30.0  7.475e+02 1000-1025,1026-1050
# sscyg_kgrid/sscyg_04_04_00_01_02_02_00.spec 60.0  1.413e+02 1140-1185
# sscyg_kgrid/sscyg_01_00_01_01_00_01_00.spec 15.0  9.020e+02 920-960,960-1000,1000-1025,1026-1050,1140-1185
   904.500  2.433e-12  1.569e-11  2.379e-11  0    1.198e-11  0  1.797e-11  0  2.537e-11  0  1.477e-11  0  1.335e-11  0
   904.600  1.222e-12  6.032e-12  2.144e-11  0    1.108e-11  0  1.749e-11  0  2.287e-11  0  1.378e-11  0  1.354e-11  0
   904.700  7.269e-13  4.312e-12  2.085e-11  0    1.156e-11  0  1.735e-11  0  2.224e-11  0  1.462e-11  0  1.394e-11  0

The comment lines describe the various fits, and then each fit (and the region fitted) correspond to the 
various spectra

  Notes:


pyfit3 simply calculates the best normilization in each region and then
the chi**2 for each spectrum.  In terms of mattchin the model spectra to
the data, it uses the standard kslfit approach of binning the model spectra
into the binning of the data spectra.  

There is no allowance for reddening or anything else, although this would
not in principle be difficult to add.

It is not entirely obvious that pure chi**2 is the best apporach, but that
is the one used currently.

!! It would be desirable to print out the parameters of the various models
in the results.out file, but this would require storing this information
as the fits improved or finding the models in the list.  

  History:
090408	ksl	Adapted from pyfit3

 ************************************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


#include "log.h"



#define LINELENGTH 200
#define NCOMMENTS   10
#define NWAVES      14000
#define MAX_MOD    10
#define BIG        1.e32
#define MAXTOT     100000
#define NREG	10
#define MAXVAR	10



char waverange[NREG][LINELENGTH];

// frac_weight is reg_weight renormalized to 1
float reg_weight[NREG],
  frac_weight[NREG];

float overall, results[NREG];


#include "py_compare.h"
#include "templates.h"

int
main (argc, argv)
     int argc;
     char *argv[];
{

  FILE *optr, *fopen ();
  char parfile[LINELENGTH];
  char diagfile[LINELENGTH];

  char first_model[LINELENGTH];
  char second_model[LINELENGTH];
  int first_spec, second_spec;
  struct Spec first_spectrum, second_spectrum;

  char outname[LINELENGTH];
  char outputfile[LINELENGTH];
  char output_specfile[LINELENGTH];
  int i;
  int n, nreg;
  float x;


/* Initialize variables that will be read by rdpar.  This needs to be
 * some of these. Note that some of the names are important, because
 * they are used to determine whether the user has provided entries on
 * the command line*/

  for (n = 0; n < NREG; n++)
    {
      strcpy (waverange[n], "1525-1575");
      reg_weight[n] = 1.0;
    }
  strcpy (outname, "py_compare");


  strcpy (first_model, "Unknown");
  strcpy (second_model, "Unknown");
  first_spec = 0;
  second_spec = 0;



/* Parse the command line */
  if (argc == 1)
    {
      strcpy (parfile, "py_compare.pf");
    }
  else
    {
      for (i = 1; i < argc; i++)
	{
	  if (strcmp (argv[i], "-h") == 0)
	    {
	      printf ("Usage: py_compare [file.pf]\n");
	      exit (0);
	    }
	  else if (strcmp (argv[i], "-m1") == 0)
	    {
	      if (sscanf (argv[i + 1], "%s", first_model) != 1)
		{
		  Error
		    ("py_compare: Expected first model after -m switch\n");
		  exit (0);
		}
	      i++;
	      if (sscanf (argv[i + 1], "%d", &first_spec) != 1)
		{
		  Error
		    ("py_compare: Expected first spectrum after modelname\n");
		  exit (0);
		}
	      i++;

	    }
	  else if (strcmp (argv[i], "-m2") == 0)
	    {
	      if (sscanf (argv[i + 1], "%s", second_model) != 1)
		{
		  Error
		    ("py_compare: Expected first model after -m switch\n");
		  exit (0);
		}
	      i++;
	      if (sscanf (argv[i + 1], "%d", &second_spec) != 1)
		{
		  Error
		    ("py_compare: Expected first spectrum after modelname\n");
		  exit (0);
		}
	      i++;

	    }
	  else if (strcmp (argv[i], "-o") == 0)
	    {
	      if (sscanf (argv[i + 1], "%s", outname) != 1)
		{
		  Error
		    ("py_compare: Expected root for output files  after -o switch\n");
		  exit (0);
		}
	      i++;

	    }
	  else if (strncmp (argv[i], "-", 1) == 0)
	    {
	      Error ("python: Unknown switch %s\n", argv[i]);
	      exit (0);
	    }

	}


      /* The last argument is always the pf file unless there were no arguments */

      strcpy (parfile, argv[argc - 1]);
      if ((strstr (parfile, ".pf")) == NULL)
	strcat (parfile, ".pf");

    }

  /* Construct an outputfile name if it is not given but model list or source comes from command line */



/* Open a log file */

  get_root (diagfile, parfile);
  strcat (diagfile, ".log");
  Log_init (diagfile);

/* Get data from a parameter file */


  opar (parfile);


  if (strcmp (first_model, "Unknown") == 0)
    {
      rdstr ("First_model", first_model);
      rdint ("First_spectrum", &first_spec);
    }
  if (strcmp (second_model, "Unknown") == 0)
    {
      rdstr ("Second_model", second_model);
      rdint ("Second_spectrum", &second_spec);
    }

  rdint ("number.of.independent.regions", &nreg);
  for (n = 0; n < nreg; n++)
    {
      rdstr ("wavelengths.to.fit", waverange[n]);
      rdflo ("chi2.renorm", &reg_weight[n]);
    }

  if (strcmp (outname, "py_compare") == 0)
    {
      get_root (outname, parfile);

    }

  strcpy (outputfile, outname);
  strcat (outputfile, ".out");
  strcpy (output_specfile, outname);
  strcat (output_specfile, ".spec");

  cpar (parfile);

/* End of input section */

  /* Compute the fractional weights */
  x=0.0;
  for (n=0;n<nreg;n++){
	  x+=reg_weight[n];
  }
  for (n=0;n<nreg;n++){
	  frac_weight[n]=reg_weight[n]/x;
  }

  /* Read the models */

  if (!get_python_spec (first_model, first_spec, &first_spectrum))
    {
      Error ("Could not read %s spectrum %d \n", first_model, first_spectrum);
      exit (0);
    }


  if (!get_python_spec (second_model, second_spec, &second_spectrum))
    {
      Error ("Could not read %s spectrum %d \n", second_model,
	     second_spectrum);
      exit (0);
    }

  overall = 0;

  for (n = 0; n < nreg; n++)
    {
      parse_waverange (&first_spectrum, waverange[n]);

      /* Now do the fit */

      x = metric (&first_spectrum, &second_spectrum);

      /* Now report the results somehow */

      results[n] = x;

      overall += frac_weight[n] * x;



    }


  /* Report the results */

  optr=fopen(outputfile,"a");


  printf ("%20s %2d  %20s %2d %8.3f", first_model, first_spec, second_model,
	  second_spec, overall);
  for (n = 0; n < nreg; n++)
    {
      printf (" %8.3f ", results[n]);
    }

  printf ("\n");



  fprintf (optr,"%20s %2d  %20s %2d %8.3f", first_model, first_spec, second_model,
	  second_spec, overall);
  for (n = 0; n < nreg; n++)
    {
      fprintf (optr," %8.3f ", results[n]);
    }

  fprintf (optr,"\n");


  fclose(optr);


  exit (0);

}





/*
 * This version adapted to change what is passed to
 * the routine 
 */

int
parse_waverange (spec, waverange)
     SpecPtr spec;
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
      nintervals++;
      if ((xxx = index (xxx, ',')) == NULL)
	break;
      xxx++;
    }


  for (n = 0; n < spec->nw; n++)
    spec->good[n] = 0;

  spec->ng = 0;
  for (n = 0; n < spec->nw; n++)
    {
      www = spec->w[n];
      for (i = 0; i < nintervals; i++)
	{
	  // Next line modified to allow overlapping waveranges
	  if (w1[i] <= www && www <= w2[i] && spec->good[n] == 0)
	    {
	      spec->good[n] = 1;
	      spec->ng++;
	    }
	}
    }

  if (spec->ng < 1)
    {
      Error
	("parse_wavelength_range:  No data points in waverange %s\n",
	 waverange);
      exit (0);
    }

  return (0);

}


float
metric (one, two)
     SpecPtr one, two;
{
  int n;
  float x, y;

  x = 0;
  y = 0;

  for (n = 0; n < one->nw; n++)
    {
      if (one->good[n])
	{
	  x += fabs (one->f[n] - two->f[n]);
	  y += 0.5 * (one->f[n] + two->f[n]);
	}
    }

  if (y > 0.0)
    {
      return (x / y);
    }
  return (1.0);

}
