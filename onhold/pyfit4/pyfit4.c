


/**************************************************************************
                        Space Telescope Science Institute


  Synopsis:  pyfit4

  Description:	
  	This is a routine to fit a spectrum to a series of python models.
	It is essentially a test of my new model reading routines

  Arguments:  (Input via .pf file)		

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
080925	ksl	This is a heavily revised version of pyfit3, which takes
		advantage of a new set of model reading routines.  Much
		of the documentation in the comments above is now 
		inaccurate concerning the format of the input and output
		files.  The current version can read a considerably wider
		range of models than just python models.

Bugs -- There is currently no way to limit the models being fit in a
	file.  A more generalized procedure is now needed.

	The output files do not feed my old progrms for interpretation
	since the headers to the output file, particularly the results
	file are not constructed properly

 ************************************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


#include "log.h"
#include "models.h"



#define MAX_MOD    10
#define BIG        1.e32
#define NREG	10

struct Spectrum mod_dummy[MAX_MOD], mod[MAX_MOD], reg_mod_best[NREG],
  mod_best;

struct Model absorp;


/* The results of the individual fits are stored here, note
 * that this structure also contains all the information 
 * on the input files */

struct res
{
  char name[LINELEN];
  float par[NPARS];
  int npars;			/* Really not needed as should be the
				   same for all of the models, but just
				   for convenience
				 */
  float norm, reg_norm[NREG];   /* norm is a weighted normalization for various spectra
				   reg_norm is the normalization for the individual regions */
  float chinu, reg_chinu[NREG]; /* chinu is the weighted normalization; reg_chinu, the
				   individual chi_nu */
  int nreg;                     /* Again not really needed, but this structure doesn't 
				   take much space
				   */
}
results[NMODS];




/* array where the variious regions are stored */
char waverange[NREG][LINELEN];

// int nres, npar;

float xx_bestnorm, xx_normchi;
float normchi_best, reg_normchi_best[NREG], reg_weight[NREG],
  frac_weight[NREG];

/* frac_weight is reg_weight renormalized to 1
 */


int
main (argc, argv)
     int argc;
     char *argv[];
{

  FILE *optr, *fopen ();
  char parfile[LINELEN];
  char source[LINELEN];
  char outname[LINELEN];
  char outputfile[LINELEN];
  char output_specfile[LINELEN];
  char model_list[LINELEN];
  char absorption_file[LINELEN];
  int i, m;
  int n, nreg;
  int iabs;
  float x;
  float incmin, incmax;
  int nmod, modstop;
  CompPtr one_comp;
  ModelPtr one_model;
  double timer();
  float chisq_pure ();
  int nres;

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
  strcpy (absorption_file, "none");
  for (n = 0; n < NREG; n++)
    {
      strcpy (waverange[n], "1525-1575");
      reg_normchi_best[n] = BIG;
      reg_weight[n] = 1.0;
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
  rdstr ("absorption_file(none=no_abs)", absorption_file);
  rdflo ("minimum.inclination.angle", &incmin);
  rdflo ("maximum.inclination.angle", &incmax);
  rdint ("number.of.independent.regions", &nreg);
  for (n = 0; n < nreg; n++)
    {
      rdstr ("wavelengths.to.fit", waverange[n]);
      rdflo ("chi2.renorm", &reg_weight[n]);
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


  /* Initiate a timer */

  Log ("Starting timer: %f\n",timer());



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


  /* Now get the models */

  get_models (model_list, "disk", 1);	// We have to actually get the data here, because we are not
  // going through one_model
  one_comp = get_component ("disk");	// Not actually necessary; we know it's comp 0


/* Record the fact that we hare read all the models */

  Log("Finished reading models %f\n",timer());

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

  modstop = one_comp->modstop;

  for (nmod = one_comp->modstart; nmod < modstop; nmod++)
    {
      one_model = &mods[nmod];

      /* Set up the index for storing the results of the fits in results */
      if (nmod == one_comp->modstart)
	{
	  nres = 0;
	}
      else
	{
	  nres++;
	}

      /* The first thing we need to do is multiply each model by the absorption */
      if (iabs)
	{
	  model_multiply (one_model);
	}


      /* Now populate the results record results[nres] with the model 
       * information
       */

      strcpy (results[nres].name, one_model->name);
      results[nres].npars = one_comp->npars;
      results[nres].nreg = nreg;

      for (i = 0; i < results[nres].npars; i++)
	{
	  results[nres].par[i] = one_model->par[i];
	}

      /* Now do the model fits for each region in turn */

      for (n = 0; n <= nreg; n++)
	{
	  // Instead of parsing the waveranges every time just copy

	  for (i = 0; i < spec.nw; i++)
	    spec.good[i] = reg_mod_best[n].good[i];
	  spec.ng = reg_mod_best[n].ng;


	  chisq_pure (one_model);	// Actually fit the model to the data 
	  xx_normchi /= spec.ng;	// Convert to chi**2/dof

	  /* Store the reults of the individual fits */

	  results[nres].reg_chinu[n] = xx_normchi;
	  results[nres].reg_norm[n] = xx_bestnorm;

	  /* Collect the information that is needed for the best overall fit.  
	   * Note that although it appears that we are including the fit that
	   * contains the concatenated set of regions, its weight is zero. */

	  results[nres].chinu += reg_weight[n] * xx_normchi/ nreg;
	  results[nres].norm += frac_weight[n] * xx_bestnorm;

	  if (xx_normchi < 0.0)
	    {
	      Error ("pyfit4: xx_normchi<0.0 %f %s\n", xx_normchi,
		     one_model->name);
	    }

	  /* Now trap and store individual best spectra  */

	  if (xx_normchi < reg_normchi_best[n])
	    {
	      strcpy (reg_mod_best[n].name, one_model->name);
	      reg_mod_best[n].npars = one_comp->npars;
	      for (i = 0; i < one_comp->npars; i++)
		{
		  reg_mod_best[n].par[i] = one_model->par[i];
		}
	      for (i = 0; i < spec.nw; i++)
		{
		  reg_mod_best[n].f[i] = xx_bestnorm * one_model->f[i];

		}
	      reg_normchi_best[n] = xx_normchi;

	    }
	}


      /* Now trap the overall best fit */

      if (results[nres].chinu < normchi_best)
	{
	  strcpy (mod_best.name, one_model->name);
	  normchi_best = results[nres].chinu;
	  for (i = 0; i < one_comp->npars; i++)
	    {
	      mod_best.par[i] = one_model->par[i];
	    }
	  for (i = 0; i < spec.nw; i++)
	    {
	      mod_best.f[i] = results[nres].norm * one_model->f[i];

	    }
	}


    }

  /* All of the work is done at this point, and stored in various structures
   * The remainder of the program writes the results to files */


  fprintf (optr,
	   "# Col: %2d  %50s %s %10.3e   # The-best-weighted-compromise\n",
	   1, waverange[nreg], mod_best.name, normchi_best);



  /* Now sumarize the results for the indivicual fits, including the combined region */

  for (n = 0; n <= nreg; n++)
    {
      fprintf (optr, "# Col: %2d  %50s %s %10.3e \n", n + 2,
	       waverange[n], reg_mod_best[n].name, reg_normchi_best[n]);


    }


  for (i = 0; i < one_comp->nmods; i++)
    {
      write_record (optr, i);
    }

    /* Finally, write out the best fitting models to the .out file. Include
     *    some comment lines that describe the models  */


    write_output_spectrum (output_specfile, nreg);

    Log ("Finished at: %f\n",timer());

  exit (0);
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
  char dummy[LINELEN];
  int nwaves;
  float w[NWAVES], f[NWAVES], f_out[NWAVES];
  int nn;

  struct Model w_in, w_out;


  /* Get the data */
  if ((ptr = fopen (filename, "r")) == NULL)
    {
      Error ("get_abs: Could not open filename %s\n", filename);
      exit (0);
    }

  nwaves = 0;
  while (nwaves < NWAVES && (fgets (dummy, LINELEN, ptr)) != NULL)
    {
      if ((dummy[0] != '#'))
	{
	  sscanf (dummy, "%e %e", &w[nwaves], &f[nwaves]);
	  nwaves++;
	}
    }

  fclose (ptr);


  /* Now regrid this model */


  w_in.w = w;
  w_in.f = f;
  w_in.nwaves = nwaves;

  w_out.w = spec.w;
  w_out.nwaves = spec.nw;
  w_out.f = f_out;

  regrid (&w_in, &w_out);

  /* So now we have everything gridded the way we want it, but we need to make it
   * permanent
   */


  assign_wavelengths (&absorp, w_out.nwaves, w_out.w);
  absorp.f = calloc (w_out.nwaves, sizeof (float));

  nn=0;
  while (nn < w_out.nwaves)
    {
      absorp.f[nn] = w_out.f[nn];
      nn++;
    }

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
0809	ksl	Changed this to model 


 ************************************************************************/

int
model_multiply (model)
     ModelPtr model;
{
  int n;

  for (n = 0; n < model->nwaves; n++)
    {
      model->f[n] *= absorp.f[n];
    }
  return (0);
}



float
chisq_pure (mod)
     ModelPtr mod;
{
  double xx, xz, zz, ee, eee;	// These are double to increase exponent range
  double xrenorm;

  float c;
  int n;

  xx = xz = zz = 0.0;

  /* Get the average value of the spectra and the model ignoring the errors */
  for (n = 0; n < spec.nw; n++)
    {
      if (spec.e[n] > 0.0 && spec.good[n])
	{
	  xz += spec.f[n] * mod->f[n];
	  zz += mod->f[n] * mod->f[n];
	}
    }

  xrenorm = xz / zz;		// This is approximately the ratio of the data to the model

  for (n = 0; n < spec.nw; n++)
    {
      if (spec.e[n] > 0.0 && spec.good[n])
	{
	  eee = spec.e[n];
	  ee = 1. / (eee * eee);
	  xx += spec.f[n] * spec.f[n] * ee;
	  xz += spec.f[n] * xrenorm * mod->f[n] * ee;
	  zz += xrenorm * mod->f[n] * xrenorm * mod->f[n] * ee;
	}
    }
  /* Now calculate the normalized chi**2 */

  xx_bestnorm = xz / zz;
  xx_normchi =
    xx - 2. * (xx_bestnorm) * xz + (xx_bestnorm) * (xx_bestnorm) * zz;

  /* Now correct the normalization */

  xx_bestnorm *= xrenorm;


  /* Finally calculate the straight chi**2, which is not likely to be that useful */

  c = xx - 2. * xz + zz;
  if (c > BIG)
    c = BIG;

  /* Finally just check that we are not ending up with a normalization which is negative */

  if (xx_normchi < 0)
    {
      Error
	("chisq:xx %8.3g xz %8.3g zz %8.3g bestnorm %8.3g normchi %8.3g chisq %8.3g \n",
	 xx, xz, zz, xx_bestnorm, xx_normchi, c);

    }
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

int
write_record (optr, nrec)
     FILE *optr;
     int nrec;
{
  int n;
  /* Col 1 containts the ratio of the weighted fit for this model to the best fit */
  fprintf (optr, "%8.3f  ", results[nrec].chinu / (normchi_best));

  /*  Cols 2 through nreg+1 contain the ratio of the weighted fit for this region to the best fit for this reg */

  for (n = 0; n <= results[nrec].nreg; n++)
    fprintf (optr, " %8.3f",
	     results[nrec].reg_chinu[n] / reg_normchi_best[n]);
  fprintf (optr, "   ");	/* space to spearate fits from parameters */


  /* Columns nreg+2 thorough nreg+npar+1 contain the parameters used to create this model */

  for (n = 0; n < results[nrec].npars; n++)
    fprintf (optr, " %9.2e", results[nrec].par[n]);

  /* 
   Column nreg+npar+2 contains the average chinu, the weighted chinu
   Column nreg+npar+3 contains the weighted normalization for this model (based on the fit to the regions (separately)
   */
  fprintf (optr, "	%9.2e %9.2e ",
	   results[nrec].chinu, results[nrec].norm);

  /* Column nreg+npar+4 contains the model name */
  fprintf (optr, " %30s\n", results[nrec].name);

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

  fprintf (optr, "# %s %10.3e\n", mod_best.name, normchi_best);

  printf ("# %s %10.3e\n", mod_best.name, normchi_best);

  for (n = 0; n <= nreg; n++)
    {
      fprintf (optr, "# %s %10.3e %s\n", reg_mod_best[n].name,
	       reg_normchi_best[n], waverange[n]);
      printf ("# %s %10.3e %s\n", reg_mod_best[n].name,
	      reg_normchi_best[n], waverange[n]);


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
