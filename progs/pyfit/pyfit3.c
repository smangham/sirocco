


/**************************************************************************
                        Space Telescope Science Institute


  Synopsis:  pyfit3

  Description:	
  	This is a routine to fit a spectrum to a series of python models.
	The intent is to fit a series of regions (lines) rather than
	the overall spectral shape.  An overall spectral shape is also
	fit from the sum of all the accepted regions, but this does not contribute to 
	the "best spectrum".  pyfit3 treats each inclination angle in the python 
	output files separately; it does not interpolate to produce a best fit for
	and intermediate inclination.

  Arguments:  (Input via .pf file or switches from command line (See help
  below)		

  The parameter file looks something like this 

source_file                             sscyg_fuse4    !spectrum to fit
model_list                      sscyg_kgrid/models.ls  !list of input models/params
absorption_file                          absorp.dat    !Optional absorption file
minimum.inclination.angle                        0     
maximum.inclination.angle                        90
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


  The input file ls should looks something like this, where the first variable is the 
  location of the python spec (or spec_smo) file.

sscyg_kgrid/sscyg_00_00_00_00_00_00_00.spec 5.00e-10 5.00e-11 4 0.00 1.00 1.0 1.00e+10
sscyg_kgrid/sscyg_00_00_00_00_00_00_01.spec 5.00e-10 5.00e-11 4 0.00 1.00 1.0 3.00e+10
sscyg_kgrid/sscyg_00_00_00_00_00_00_02.spec 5.00e-10 5.00e-11 4 0.00 1.00 1.0 7.00e+10


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
the chi**2 for each spectrum.  In terms of matching the model spectra to
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
2003	ksl	Adapted from kslfit
03Oct29	ksl	Modify so that can keep track of the indvidual bestfits
		as well as the global one.  This is intended to speed 
		thing up, but the programming style is hindered by the 
		fact that I cannot get everything into memory at once.
03nov	ksl	Modified parse waverange so that it would be able to 
		handle overlapping wavelength ranges properly
		Added error check to quit if too many models were read in
		Added fit to sum total of all regions
04dec	ksl	Modified what to the results file, so that instead of 
		the pure chi**2 one was writing the ration of the particular
		fit to the best fit for each element
04dec	ksl	Modified to accept the results of py_smooth, and new naming
		conventions in py54
06aug	ksl	Added the possibility of absorption for a second time!  I
		had apparently made this modification on perch earlier this
		year but somehow overwrote it.
070707	ksl	Added a basic capability to keep track of what the columns 
		mean in the input files.  This program is a kluge!
070804	ksl	Modified to handle the format of the .ls file for both
		old-style versions of the spectra, e.g sscyg_kgrid, and
		newer ones.  It's still a mess.
070908	ksl	Modified so only handles the newest versions of outputs
		produced by gen_grid.py.  This version has comment lines
		in the .ls file that use the keyword Variable.  Older
		versions need to be convered to this.    I also tried
		to clean up the way the columns are commented
070917	ksl	Made additional changes so that best fit models are
		summarized at the beginning of outputs file
090114	ksL	Added check to see that there were not more variables than
		allowed for by MAXVAR, and changed the way comments wer
		searched for.
090120	ksl	Minor change to record outputs in a log file
090303	ksl	Added changes to allow program to be run from command line.
090608	ksl	Added help as in python
090609	ksl	Added normalized absolute value as a metric for the fit


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

struct spect
{
  char name[LINELENGTH];
  float ang;
  float w[NWAVES], w_edge[NWAVES];
  float f[NWAVES];
  float e[NWAVES];
  int good[NWAVES];
  int nw, ng, icol;
  int nres;			// An index for the best models into res structure
}
spec, mod_dummy[MAX_MOD], mod[MAX_MOD], reg_mod_best[NREG], mod_best, absorp;


/* The results of the individual fits are stored here, note
 * that this structure also contains all the information 
 * on the input files */
struct res
{
  char name[LINELENGTH];
  float xpar[MAXVAR];
  float ang;
  int icol;
  float norm, reg_norm[NREG];
  float chinu, reg_chinu[NREG];
}
results[MAXTOT];

/* The next two arrays contain the values of the variables that were varied
 * and the number of variations of each variable 
 */

float unique[MAXVAR][10];
int numvar[10];
int init_unique = 1;


char waverange[NREG][LINELENGTH];

int nres, npar;

float xx_bestnorm, xx_normchi;
float normchi_best, reg_normchi_best[NREG], reg_weight[NREG],
  frac_weight[NREG];
// frac_weight is reg_weight renormalized to 1

#include "pyfit3_templates.h"

int
main (argc, argv)
     int argc;
     char *argv[];
{

  FILE *fptr, *optr, *fopen ();
  char parfile[LINELENGTH];
  char diagfile[LINELENGTH];
  char source[LINELENGTH];
  char outname[LINELENGTH];
  char outputfile[LINELENGTH];
  char output_specfile[LINELENGTH];
  char model_list[LINELENGTH];
  char line[LINELENGTH];
  char modelname[LINELENGTH];
  char absorption_file[LINELENGTH];
  char comments[NCOMMENTS][LINELENGTH], word[LINELENGTH];
  char command[LINELENGTH];
  int ncomments;
  int i, m, mspec;
  int n, nreg;
  int iabs;
  float x;
  float incmin, incmax;
  int nmods;
  int nnn, ii;
  float xx[MAXVAR], xdummy;
  int fit_type;
  double timer ();


/* Initialize variables that will be read by rdpar.  This needs to be
 * before the parsing of the command line because I may want to reset
 * some of these. Note that some of the names are important, because
 * they are used to determine whether the user has provided entries on
 * the command line*/

  strcpy (source, "source_file.dat");
  strcpy (model_list, "model_list.ls");
  strcpy (absorption_file, "none");
  for (n = 0; n < NREG; n++)
    {
      strcpy (waverange[n], "1525-1575");
      reg_normchi_best[n] = BIG;
      reg_weight[n] = 1.0;
    }
  strcpy (outname, "output_file_root");

  normchi_best = BIG;
  nres = 0;
  incmin = 0.0;
  incmax = 90.;
  fit_type = 0;


/* Parse the command line */
  if (argc == 1)
    {
      strcpy (parfile, "pyfit.pf");
    }
  else
    {
      for (i = 1; i < argc; i++)
	{
	  if (strcmp (argv[i], "-h") == 0)
	    {
	      help ();
	      exit (0);
	    }
	  else if (strcmp (argv[i], "-m") == 0)
	    {
	      if (sscanf (argv[i + 1], "%s", model_list) != 1)
		{
		  Error ("pyfit3: Expected model list  after -m switch\n");
		  exit (0);
		}
	      if (strcmp (outname, "output_file_root") == 0)
		{
		  strcpy (outname, "construct");
		}
	      i++;

	    }
	  else if (strcmp (argv[i], "-s") == 0)
	    {
	      if (sscanf (argv[i + 1], "%s", source) != 1)
		{
		  Error
		    ("pyfit3: Expected source spectrum  after -s switch\n");
		  exit (0);
		}
	      if (strcmp (outname, "output_file_root") == 0)
		{
		  strcpy (outname, "construct");
		}
	      i++;

	    }
	  else if (strcmp (argv[i], "-o") == 0)
	    {
	      if (sscanf (argv[i + 1], "%s", outname) != 1)
		{
		  Error
		    ("pyfit3: Expected root for output files  after -o switch\n");
		  help ();
		  exit (0);
		}
	      i++;

	    }
	  else if (strcmp (argv[i], "-f") == 0)
	    {
	      if (sscanf (argv[i + 1], "%d", &fit_type) != 1)
		{
		  Error
		    ("pyfit3: Expected integer for fit_type after -f switch\n");
		  help ();
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


  if (strcmp (source, "source_file.dat") == 0)
    {
      rdstr ("source_file", source);
    }
  if (strcmp (model_list, "model_list.ls") == 0)
    {
      rdstr ("model_list", model_list);
    }
  else
    {
      Log ("Using model_list %s from command line\n", model_list);
    }
  rdstr ("absorption_file(none=no_abs)", absorption_file);
  rdflo ("minimum.inclination.angle", &incmin);
  rdflo ("maximum.inclination.angle", &incmax);
  rdint ("number.of.independent.regions", &nreg);
  for (n = 0; n < nreg; n++)
    {
      rdstr ("wavelengths.to.fit", waverange[n]);
      rdflo ("chi2.renorm", &reg_weight[n]);
    }

  if (strcmp (outname, "output_file_root") == 0)
    {
      rdstr ("output_file_root", outname);
    }
  else if (strcmp (outname, "construct") == 0)
    {
      get_root (outname, parfile);
      strcat (outname, ".");
      get_root (word, source);
      strcat (outname, word);
      strcat (outname, ".");
      strcpy (word, "");
      strncpy (word, model_list, i = strcspn (model_list, "/"));
      word[i] = '\0';


      strcat (outname, word);
      Log ("Created output file root %s\n", outname);

    }

  strcpy (outputfile, outname);
  strcat (outputfile, ".results");
  strcpy (output_specfile, outname);
  strcat (output_specfile, ".out");

  cpar (parfile);

/* Initiate a timer */

  Log ("Starting timer: %f\n", timer ());

  /* Create symbolic links to the output files for the last time this is run
   */


  system ("rm last.results");
  system ("rm last_spec.out");

  sprintf (command, "ln -s %s  last.results\n", outputfile);
  system (command);
  sprintf (command, "ln -s %s last_spec.out \n", output_specfile);
  system (command);

/* End of input section */

  if ((optr = fopen (outputfile, "w")) == NULL)
    {
      printf ("Failed to open file %s\n", outputfile);
      exit (0);
    }



  /* Get the data */
  get_data (source);

  if (strncmp (absorption_file, "none", 4) != 0)
    {
      get_abs (absorption_file);
      iabs = 1;			//TRUE
    }
  else
    {
      iabs = 0;			// FALSE
    }


  /* Open the file containing the list of models.  Exit if it is not opened successfully */
  if ((fptr = fopen (model_list, "r")) == NULL)
    {
      printf ("Failed to open file %s\n", model_list);
      exit (0);
    }


/* Now create a renormalized set of weights for the model fitting*/
  x = 0;
  for (n = 0; n < nreg; n++)
    x += reg_weight[n];
  for (n = 0; n < nreg; n++)
    {
      frac_weight[n] = reg_weight[n] / x;
      Log ("frac_weight %d %f\n", n, frac_weight[n]);
    }


/* Now automatically create one more set of regions to fit consisting of the
 * concatenation of all of the regions */
  strcpy (waverange[nreg], waverange[0]);
  for (n = 1; n < nreg; n++)
    {
      strcat (waverange[nreg], ",");
      strcat (waverange[nreg], waverange[n]);
    }
/* Next line is important; it prevents the combined fit from being included
 * in the weighting; it is possible that one should use an input variable
 * for this
 */

  frac_weight[nreg] = reg_weight[nreg] = 0.0;

  for (n = 0; n <= nreg; n++)
    {
      parse_waverange (waverange[n]);
      reg_mod_best[n].ng = spec.ng;
      for (i = 0; i < spec.nw; i++)
	{
	  reg_mod_best[n].w[i] = spec.w[i];
	  reg_mod_best[n].good[i] = spec.good[i];
	}
    }

/* Read and proccess the individual models 

   If the number of parameters in the list exceeds that allowed in the scanf line
   modifications need to be made to the format statement.  If the number of parameters
   is less than what is shown, this is properly handled.
*/

  nmods = 0;
  ncomments = 0;
a:while (fgets (line, LINELENGTH, fptr) != NULL)
    {
      sscanf (line, "%s", word);
      if (word[0] == '#')
	{			// Trap comment lines 
	  if (ncomments < NCOMMENTS)
	    strcpy (comments[ncomments], line);
	  printf ("Comment: %s", line);
	  ncomments++;

	  goto a;
	}


      nnn =
	sscanf (line, "%s %e %e %e %e %e %e %e %e %e %e %e", modelname,
		&xx[0], &xx[1], &xx[2], &xx[3], &xx[4], &xx[5], &xx[6],
		&xx[7], &xx[8], &xx[9], &xdummy);
      if (nnn == MAXVAR + 2)
	{
	  Error
	    ("Too many variables in input line.  Increase size of MAXVAR, and rerun\n");
	  printf ("%s\n", line);
	  exit (0);

	}

      if (nmods == 0)
	npar = nnn;
      else if (npar != nnn)
	{
	  Error ("Not all models have same number of variables\n");
	  exit (0);
	}

      /* get_python_models returns the number of spectra which were
       * read.  If the file called modelname does not exist thne
       * mspec will be 0
       */

      mspec = get_python_model (modelname);

      if ((nmods % 100) == 0)
	printf ("%s\n", modelname);
      nmods++;

      /* Identify the allowed valus of variables */

      update_unique (xx, nnn);

/* Now we have the data and the models and we can begin a various comparisons
*/

      for (m = 0; m < mspec; m++)
	{
	  if (incmin <= mod[m].ang && mod[m].ang <= incmax)
	    {
	      strcpy (results[nres].name, modelname);
	      for (ii = 0; ii < nnn; ii++)
		results[nres].xpar[ii] = xx[ii];

	      results[nres].ang = mod[m].ang;
	      results[nres].icol = mod[m].icol;
	      results[nres].chinu = 0.0;
	      results[nres].norm = 0.0;
	      if (iabs)
		model_multiply (&mod[m]);
	      x = 0;
	      // Here are the individual fits which include the overall fit
	      for (n = 0; n <= nreg; n++)
		{
		  // Instead of parsing the waveranges every time just copy
		  for (i = 0; i < spec.nw; i++)
		    spec.good[i] = reg_mod_best[n].good[i];
		  spec.ng = reg_mod_best[n].ng;

		  if (fit_type == 1)
		    {
		      metric_abs (&mod[m]);
		    }
		  else
		    {

		      chisq_pure (&mod[m]);	// Actually fit the model to the data 
		      xx_normchi /= spec.ng;	// Convert to chi**2/dof
		    }

		  results[nres].reg_chinu[n] = xx_normchi;
		  results[nres].reg_norm[n] = xx_bestnorm;

		  /* Collect the information that is needed for the best overall fit */
		  results[nres].chinu += reg_weight[n] * xx_normchi;
		  results[nres].norm += frac_weight[n] * xx_bestnorm;

		  /* Now trap and store individual best spectra  */
		  if (xx_normchi < reg_normchi_best[n])
		    {
		      strcpy (reg_mod_best[n].name, mod[m].name);
		      reg_mod_best[n].ang = mod[m].ang;
		      reg_normchi_best[n] = xx_normchi;
		      reg_mod_best[n].nres = nres;	// Associate the proper results with the model
		      for (i = 0; i < spec.nw; i++)
			{
			  reg_mod_best[n].f[i] = xx_bestnorm * mod[m].f[i];

			}

		    }
		}

	      // Now do the sum that will allow one to trap the overall best
	      // Next section is for handling the best fit overall
	      // Note that although we use the weighting to control the relationship
	      // between the fits, a horrible chi**2 in one of the regions could
	      // still dominate the others.  It may be that there are better
	      // approaches for this.  04dec--ksl

	      if (results[nres].chinu < normchi_best)
		{
		  strcpy (mod_best.name, mod[m].name);
		  mod_best.ang = mod[m].ang;
		  normchi_best = results[nres].chinu;
		  mod_best.nres = nres;	// Associate the proper results with the model
		  for (i = 0; i < spec.nw; i++)
		    {
		      mod_best.f[i] = results[nres].norm * mod[m].f[i];

		    }
		}
	      nres++;
	      if (nres > MAXTOT)
		{
		  Error ("MAXTOT (number of spectra x angles) exceeded\n");
		  exit (0);
		}

	    }

	}
    }

  /* Now write out the results of all the fits to the .results file */


  fprintf (optr, "# Source_spectrum %s\n", source);
  fprintf (optr, "# Models          %s\n", model_list);

  /* First summarize the very best model */

  fprintf (optr,
	   "# Col: %2d  %50s %s %.1f %10.3e   # The-best-weighted-compromise\n",
	   1, waverange[nreg], mod_best.name, mod_best.ang, normchi_best);



  /* Now sumarize the results for the indivicual fits, including the combined region */

  for (n = 0; n <= nreg; n++)
    {
      fprintf (optr, "# Col: %2d  %50s %s %.1f %10.3e \n", n + 2,
	       waverange[n], reg_mod_best[n].name, reg_mod_best[n].ang,
	       reg_normchi_best[n]);


    }

/* Now process the comments */

  n = 0;
  for (i = 0; i < ncomments; i++)
    {

      nnn = sscanf (comments[i], "%*s Variable   %s ", word);
      if (nnn > 0)
	{			// We have caught the column names
	  fprintf (optr, "# Col: %2d  %s\n", nreg + 3 + n, word);
	  write_unique (optr, n);	// Write out the variations of this variable
	  n = n + 1;
	}
      else
	{
	  Log ("This comment that was not interpreted: %s", comments[i]);
	}

    }

  if (n + 1 != npar)
    {
      Error ("Ncols(%d) != npar  (%d)\n", n + 1, npar);
    }
  if (n == 0)
    {
      Error ("This was an old style grid model - Please fix .ls file\n");
    }


  fprintf (optr, "# Col: %2d  %s\n", nreg + npar + 2, "Inclination");


/* Finished processing the comments */
  fprintf (optr, "# Col: %2d  weighted chi_nu of this model\n",
	   nreg + npar + 3);
  fprintf (optr,
	   "# Col: %2d  normalization factor for weighted best fit for this model\n",
	   nreg + npar + 4);
  fprintf (optr, "# Col: %2d  Modelname \n", nreg + npar + 5);
  fprintf (optr, "# Col: %2d  Spectrum  \n", nreg + npar + 6);


  fprintf (optr, "# Notes\n");
  fprintf (optr,
	   "# First %d columns are not chi**2 but ratios of chi**2 for this model compared to very best \n",
	   nreg + 2);
  fprintf (optr,
	   "# Col 1 is weighted best fit; Col %d-%d are best individ; %d is unweighted best \n",
	   2, nreg + 1, nreg + 2);
  fprintf (optr,
	   "# Columns  %d - %d describe the  models; %d is inclination\n",
	   nreg + 3, nreg + npar + 2, nreg + npar + 2);

  /* Now write the best models as comments */

  fprintf (optr, "# Show results for the best fit models: \n");

  fprintf (optr, "#");
  write_record (optr, mod_best.nres, nreg);

  for (n = 0; n <= nreg; n++)
    {

      fprintf (optr, "#");
      write_record (optr, reg_mod_best[n].nres, nreg);
    }

  /* Now write the individual records to the open file */
  fprintf (optr, "#  Begin results for individual fits:\n");
  for (i = 0; i < nres; i++)
    {
      write_record (optr, i, nreg);
    }


  fclose (optr);

  /* Now write out the best fitting models to the .out file. Include
   * some comment lines that describe the models  */


  write_output_spectrum (output_specfile, nreg);

  /* Finished */
  Log ("Finished at: %f\n", timer ());
/* OK now end program and go on to the next important task. */
  exit (0);

}

/***********************************************************
Space Telescope Science Institute

Synopsis:
print out some basic help on how to run the program
Arguments:              

Returns:
 
Description:    
                
Notes:

The easiest way to create the message, at least initially, is simply to to type
out what you want to appear on the screen and then as \n\ to all of the lines, including
the ones with nothing in them

History:

090608	ksl	Added help, mainly as a way to understand the usage

**************************************************************/

int
help ()
{

  char *some_help;

  some_help = "\
\n\
This is a routine to fit a spectrum to a series of python models.\n\
The intent is to fit a series of regions (lines) rather than\n\
the overall spectral shape.  An overall spectral shape is also\n\
\n\
Usage: pyfit [-h] -m model_list -s source_spectrm  -o output_root [file.pf]\n\
\n\
where:\n\
  	-h 			produces this help file\n\
	-m model_list  		the name of file containing all the models and the variables they contain.\n\
	-s source_spectrum	Input spectrum in standard format (wave,flux,err)\n\
	-o output_root 		replaces the output_file with a new name; normally the ouput file is \n\
				computed from the model_list and spectrum name or given through the\n\
				pf file\n\
	-f fit_type		Type of fit.  The choices are:\n\
       					1 	normalized absolute value\n\
					other	normalize chi2\n\
\n\
";				// End of string to provide one with help

  printf ("%s\n", some_help);

  exit (0);
}

/* This routine is intended to collect the various values that
 * have been varied to make the models that are being fit to 
 * the data
 */

int
update_unique (xx, nnn)
     float xx[];
     int nnn;
{
  int i, n;
  int inew;
  if (init_unique)
    {
      for (n = 0; i < MAXVAR; i++)
	{
	  numvar[n] = 0;
	}
      init_unique = 0;
    }

  for (i = 0; i < nnn; i++)
    {
      inew = 1;
      for (n = 0; n < numvar[i]; n++)
	{
	  if (xx[i] == unique[i][n])
	    inew = 0;
	}
      if (inew)
	{
	  unique[i][numvar[i]] = xx[i];
	  numvar[i]++;
	}
    }
  return (0);
}

/* This writes out one of the colected variables to a comment line */
int
write_unique (optr, n)
     FILE *optr;
     int n;
{
  int i;
  fprintf (optr, "# Values in grid: ");
  for (i = 0; i < numvar[n]; i++)
    {
      fprintf (optr, "%8.2g", unique[n][i]);
    }
  fprintf (optr, "\n");
  return (0);
}

int
write_record (optr, nrec, nreg)
     FILE *optr;
     int nrec, nreg;
{
  int n;

  /* Col 1 containts the ratio of the weighted fit for this model to the best fit */
  fprintf (optr, "%8.3f  ", results[nrec].chinu / (normchi_best));
  /*  Cols 2 through nreg+1 contain the ratio of the weighted fit for this region to the best fit for this reg */
  for (n = 0; n <= nreg; n++)
    fprintf (optr, " %8.3f",
	     results[nrec].reg_chinu[n] / reg_normchi_best[n]);
  fprintf (optr, "   ");	/* space to spearate fits from parameters */
  /* Columns nreg+2 thorough nreg+npar+1 contain the parameters used to create this model */
  for (n = 0; n < npar - 1; n++)
    fprintf (optr, " %9.2e", results[nrec].xpar[n]);
  /* Column nreg+npar+2 contains the angle for this particular spectrum
     Column nreg+npar+3 contains the average chinu
     Column nreg+npar+4 contains the normalization for this model (based on the fit to the regions (separately)
   */
  fprintf (optr, "   %5.1f    %9.2e %9.2e ", results[nrec].ang,
	   results[nrec].chinu / nreg, results[nrec].norm / nreg);
  /* Column nreg+npar+5 contains the model name */
  fprintf (optr, " %30s %3d\n", results[nrec].name, results[nrec].icol);

  return (0);
}

int
write_output_spectrum (filename, nreg)
     char filename[];
     int nreg;
{
  FILE *optr, *fopen ();
  int n, i;

  if ((optr = fopen (filename, "w")) == NULL)
    {
      printf ("Failed to open file %s\n", filename);
      exit (0);
    }

  fprintf (optr, "# %s %.1f %10.3e\n", mod_best.name, mod_best.ang,
	   normchi_best);

  printf ("# %s %.1f %10.3e\n", mod_best.name, mod_best.ang, normchi_best);

  for (n = 0; n <= nreg; n++)
    {
      fprintf (optr, "# %s %.1f %10.3e %s\n", reg_mod_best[n].name,
	       reg_mod_best[n].ang, reg_normchi_best[n], waverange[n]);
      printf ("# %s %.1f %10.3e %s\n", reg_mod_best[n].name,
	      reg_mod_best[n].ang, reg_normchi_best[n], waverange[n]);


    }

  for (i = 0; i < spec.nw; i++)
    {
      fprintf (optr, "%10.3f %10.3e %10.3e %10.3e %2d  ", spec.w[i],
	       spec.f[i], spec.e[i], mod_best.f[i], spec.good[i]);
      for (n = 0; n <= nreg; n++)
	fprintf (optr, " %10.3e %2d", reg_mod_best[n].f[i],
		 reg_mod_best[n].good[i]);
      fprintf (optr, "\n");


    }


  return (0);
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
  return (0);
}





/**************************************************************************
                        Space Telescope Science Institute


  Synopsis:  
	get_abs (filename) simplye gets a wavelenght value file and 
	intepolates it onto the spectrum grid.

  Description:	

  Arguments:  (Input via .pf file)		

  Returns:


  Notes:


  History:
06aug	ksl	Added the possibility of absorption for a second time!


 ************************************************************************/


int
get_abs (filename)
     char filename[];
{
  FILE *ptr;
  char dummy[LINELENGTH];
  int n, nwaves;
  struct spect abs_dummy;


  if ((ptr = fopen (filename, "r")) == NULL)
    {
      Error ("Could not open filename %s\n", filename);
      exit (0);
    }

  nwaves = 0;
  while (nwaves < NWAVES && (fgets (dummy, LINELENGTH, ptr)) != NULL)
    {
      if ((dummy[0] != '#'))
	{
	  sscanf (dummy, "%e %e", &abs_dummy.w[nwaves], &abs_dummy.f[nwaves]);
	  nwaves++;
	}
    }
  abs_dummy.nw = nwaves;


  fclose (ptr);

  abs_dummy.w_edge[0] =
    abs_dummy.w[0] - 0.5 * (abs_dummy.w[1] - abs_dummy.w[0]);
  for (n = 1; n < nwaves; n++)
    {
      abs_dummy.w_edge[n] = 0.5 * (abs_dummy.w[n - 1] + abs_dummy.w[n]);
    }

  abs_dummy.w_edge[nwaves] =
    abs_dummy.w_edge[nwaves - 1] +
    0.5 * (abs_dummy.w_edge[nwaves - 2] + abs_dummy.w_edge[nwaves - 1]);

  /* OK, now put them on the same wavelength scale as the spectrum */

  regrid (&abs_dummy, &absorp);

  //Debug steps to verify that absorption file looks correct
//DEBUG  ptr = fopen ("abs_check.out", "w");
//DEBUG
//DEBUG  for (n = 0; n < absorp.nw; n++)
//DEBUG    {
//DEBUG      fprintf (ptr, "%f %f\n", absorp.w[n], absorp.f[n]);
//DEBUG    }
//DEBUG  fclose (ptr);

  return (0);
}



/**************************************************************************
                        Space Telescope Science Institute


  Synopsis:  

  Description:	

  Arguments:  (Input via .pf file)		

  Returns:


  Notes:


  History:
06aug	ksl	Added the possibility of absorption for a second time!
07sep	ksl	Modified so that missing files do not crash the program.
		This was because in some cases one does not always get
		every model calculated on Royal.


 ************************************************************************/

int
model_multiply (model)
     struct spect *model;
{
  int n;

  for (n = 0; n < model->nw; n++)
    {
      model->f[n] *= absorp.f[n];
    }
  return (0);
}

int
get_python_model (filename)
     char filename[];
{
  FILE *ptr;
  char dummy[LINELENGTH];
  int n, m, nscanned;
  float ang[MAX_MOD], wavelength, f[MAX_MOD];
  int mspec, nwaves;

  if ((ptr = fopen (filename, "r")) == NULL)
    {
      Error ("Could not open filename %s\n", filename);
      return (0);
//OLD      exit (0);
    }


// Get the title line

  for (m = 0; m < 6; m++)
    ang[m] = 999;

  /* Kluge to handle both old and new style spectral lines.  Originally the #Freq line was the first line.  
   * Today 07, we use that line tor record the versio of Python run.  */
  fgets (dummy, LINELENGTH, ptr);
  mspec =
    sscanf (dummy,
	    "#Freq.        Lambda  Emitted  Star+BL  Disk     %*s       %*s       %*s       Ang:%f  Ang:%f  Ang:%f  Ang:%f  Ang:%f Ang:%f",
	    &ang[0], &ang[1], &ang[2], &ang[3], &ang[4], &ang[5]);

  if (mspec == 0)
    {				// It's a new style line
      fgets (dummy, LINELENGTH, ptr);
      mspec = sscanf (dummy,
		      "# Freq.        %*s %*s %*s %*s %*s       %*s       %*s       Ang:%f  Ang:%f  Ang:%f  Ang:%f  Ang:%f Ang:%f",
		      &ang[0], &ang[1], &ang[2], &ang[3], &ang[4], &ang[5]);

      if (mspec == 0)
	{
	  Error ("Could not scan header:%s", dummy);
	  exit (0);
	}


    }

  /* Now read the rest of the file */

  n = 0;
  while (n < NWAVES && (fgets (dummy, LINELENGTH, ptr)) != NULL)
    {
      if ((dummy[0] != '#'))
	{
	  nscanned =
	    sscanf (dummy, "%*e %e %*e %*e %*e %*e %*e %*e %e %e %e %e %e %e",
		    &wavelength, &f[0], &f[1], &f[2], &f[3], &f[4], &f[5]);
	  if (nscanned != mspec + 1)
	    Error ("Problem in file %s (%d scanned/%d expected): %s",
		   filename, nscanned, mspec, dummy);
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
      mod[m].icol = m;
      strcpy (mod[m].name, filename);
      mod_dummy[m].nw = mod[m].nw = nwaves;
    }

/* OK, now invert the order of the indices.*/

  for (m = 0; m < mspec; m++)
    {
      if (mod[m].w[0] > mod[m].w[1])	// Then the wavelenght order needs fixing 
	{
	  for (n = 0; n < nwaves; n++)
	    {
	      mod_dummy[m].w[n] = mod[m].w[nwaves - 1 - n];
	      mod_dummy[m].f[n] = mod[m].f[nwaves - 1 - n];
	    }
	}
      else
	{
	  for (n = 0; n < nwaves; n++)	// They are increasing wavelength order
	    {
	      mod_dummy[m].w[n] = mod[m].w[n];
	      mod_dummy[m].f[n] = mod[m].f[n];
	    }
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

int
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

  modout->nw = spec.nw;
  return (0);
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

  /* c is pure chi**2, whereas xx_normchi is the normalized chi**2 
   * Although the routine passes pure chi**2 back that is not what
   * is being used.  This is because generally that is what we want.
   * The exact approach probably reflects a similar routine in kslfit
   * but it is surely awkward
   *
   * 090609 ksl  */
  c = xx - 2. * xz + zz;
  if (c > BIG)
    c = BIG;
  return (c);

}


/* Various portions of this routine, and in particular the way it
 * passes values back through xx_normchi, and xx_bestnorm are
 * for backwards compatiblity
 */

float
metric_abs (mod)
     struct spect *mod;
{
  double x, y, xav, yav;
  double ratio;
  double xtest;
  int n;

  x = 0;
  y = 0;
  xav = 0;
  yav = 0;
  xtest=0;

  for (n = 0; n < spec.nw; n++)
    {
      if (spec.e[n] > 0.0 && spec.good[n])
	{
	  xav += spec.f[n];
	  yav += mod->f[n];
	}
    }

  /* Pass bestnorm back so we can use it */
  xx_bestnorm = ratio = xav / yav;

  for (n = 0; n < spec.nw; n++)
    {
      if (spec.e[n] > 0.0 && spec.good[n])
	{
	  x += fabs (spec.f[n] - ratio * mod->f[n]);

	  xtest +=  (spec.f[n] - ratio * mod->f[n]);
	}


    }

  if (xav > 0.0)
    {
      xx_normchi = x / xav;
      return (xx_normchi);
    }

  xx_normchi = 1.;
  return (1.0);



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
	  // Next line modified to allow overlapping waveranges
	  if (w1[i] <= www && www <= w2[i] && spec.good[n] == 0)
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
      Error
	("parse_wavelength_range:  No data points to fit for waverange %s\n",
	 waverange);
      exit (0);
    }

  return (0);
}
