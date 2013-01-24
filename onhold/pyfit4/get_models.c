


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:

  Description:	These are a generic set of routines to read a grid or grids
	of models  The grid can be any reasonable number of dimension

  Arguments:		

  Returns:

  Notes:
	DO NOT INCLUDE ANYTHING THAT IS NOT GENERIC HERE!  

	These routines were written for kslfit, but should be quite generic.
	The python version does not regrid to a predifined wavelength
	scale.

	For 1d models we assume the models are read in increasing order of
	the parameter of interest.  

	2 or more parameter models can be read in any order.  Bilinear interpolation
	is done for the 2 parmeter models (as long as the grid is fully filled).
	To deal with cases in which the grid is not fully filled, we follow 
	the same approach as for get_hub, i.e. for parameters t and g (i.e.
	the first parameter being t, and the second g., we first search the grid
	for all models with t just below (above) the desired value.  We then 
	interpolate on g for the model just below, and on g for the model just
	above, and finally interpolate between these two.  
	The algorithm for interpolating between models 

	080928 - This is the kslfit version of these routines, and this version
       has been updated very significantly. The goal is to end up with a single
       version that can be used everywhere but this may be difficult.


  History:
04jul	ksl	Adapted for use in python
04dec	ksl	54a -- miniscule change to eliminate warning with 03
0809	ksl	This is a modification of the routine used in python.  The
		goal is ulitmately to end up with a single routine that 
		can be generically used to handle models that neeed to be
		read in and interpolated.  
	

 ************************************************************************/
#include	<math.h>
#include	<stdio.h>
#include        <stdlib.h>
#include	<strings.h>
#include	<string.h>
// #include     "atomic.h"
// #include     "python.h"      //This needs to come before modles.h so that what is in models.h is used
#include	"log.h"

#include         "models.h"
#include	"model_prototypes.h"
#define    	BIG 1e32

#define 	LOCATE_ERR 	-99

/* Get all the models of one type and regrid them onto the wavelength grid of the data */
int get_models_init = 0;	// needed so can initialize nmods_tot the first time this routine is called

int
get_models (modellist, modelname, ichoice)
     char modellist[];		// filename containing location and associated parameters of models
     char modelname[];		// Descriptor for the model
     int ichoice;		// 0 means defer reading the models until later; 1 means read the models now


{
  char *dummy;
  int n, m, nxpar;
  float xpar[NPARS], xmin[NPARS], xmax[NPARS];
  int nw, nwaves;
  int npars;
  int model_locate ();
  char word[LINELEN], filetype[LINELEN];
  int ix, iy;
  int nline;
  int check_comment ();


  nwaves = 0;
  /* Assume we are dealing with a simple set of spectral files */
  strcpy (filetype, "Simple");
  ix = 1;
  iy = 2;

  /* Initialize the structures to known values if this is the first time get_models is called */

  if (get_models_init == 0)
    {
      init_models ();
      get_models_init = 1;
    }


  if ((n = model_locate (modelname)) != LOCATE_ERR)
    {
      Log ("No need to read %s again\n", modelname);
      return (n);		// This set of models has been read in previously
    }
  else
    {
      Log ("Creating %s\n", modelname);
    }

  ascii_read (modellist);

  if (check_comment ("File_type", filetype))
    {
      strcpy (filetype, "Simple");
    }

  if (strcmp (filetype, "Python") == 0)
    {
      n=get_python_models (modellist, modelname, ichoice);
      return (n);
    }

  Log
    ("get_models: Assuming %s with filetype %s is trying to read simple models\n",
     modellist, filetype);


  n = ncomps;			// We want to add a new component to the end of the components


/* Now initialize the model summary structure */
  strcpy (comp[ncomps].name, modellist);
  strcpy (comp[ncomps].var, modelname);
  comp[ncomps].modstart = nmods_tot;
  comp[ncomps].filetype = SIMPLE;

  npars = -1;			// Don't yet know how many pars there are for this model type
//  comp[ncomps].xpdf.limit1 = -99.;
//  comp[ncomps].xpdf.limit2 = -99.;

  n = nmods_tot;		// This is the starting point since may have read models in before
/* Establish initial limits on xmin and xmax so that they can be properly populated */
  for (m = 0; m < NPARS; m++)
    {
      xmin[m] = (BIG);
      xmax[m] = (-BIG);
      comp[ncomps].xmod.par[m] = -99;
    }

  nw = -1;			// Initiallize nw
  nline = 0;
  while (n < NMODS && nline < ascii_nlines)
    {
      dummy = ascii_store[nline];

      if (dummy[0] == '#' || dummy[0] == '!')
	{
	  if (sscanf (dummy, "%*s Filetype %s", word) == 1)
	    {
	      strcpy (filetype, word);
	    }
	}			//skip comment lines in models
      else
	{

		/* In order to allow for a normalization, as for use with kslfit, the
		 * parameters start with xpar[1], and not xpar[0]
		 * ??? This will be need to reconciled in the versions in pyfit3 and python
		 */

	  nxpar =
	    sscanf (dummy, "%s %f %f %f %f %f %f %f %f %f ",
		    mods[n].name, &xpar[1], &xpar[2], &xpar[3], &xpar[4],
		    &xpar[5], &xpar[6], &xpar[7], &xpar[8], &xpar[9]);

	  xpar[0]=1.0; // Assure xpar[0] is not undefined



	/* As defined here npars includes the normalization */
	  if (npars < 0)
	    {
	      comp[ncomps].npars = npars = nxpar;	// Initialize to value in first good model

	    }

	  if (nxpar < npars)
	    {
	      Error ("get_models: nxpar (%d) < npars (%d) in line %s\n",
		     nxpar, npars, dummy);
	      exit (0);
	    }

	  for (m = 0; m < npars; m++)
	    {
	      mods[n].par[m] = xpar[m];
	      if (xpar[m] > xmax[m])
		xmax[m] = xpar[m];
	      if (xpar[m] < xmin[m])
		xmin[m] = xpar[m];
	    }
	  for (m = m; m < NPARS; m++)
	    mods[n].par[m] = -99;

	  mods[n].ix = ix;
	  mods[n].iy = iy;
	  mods[n].ncomp = ncomps;	// Associate the model with a component


	  if ((n % 100) == 0)
	    printf ("Model n %d %s\n", n, mods[n].name);

	  /* It is necessary to read at least one model to be able to create the wavelength and flux grid in comp.
	   * Or we would need to check when we tried to populate if the first time.  However unless ichoice is 0
	   * we also read all of the models at this point  */

	  if (n == 0 || ichoice)
	    nwaves = get_one_model (mods[n].name);

	  n++;
	  nmods_tot = n;
	  nline++;
	}
    }

  if (n == NMODS)
    {
      Error
	("get_models: Reached maximum number of models %d. Please increase NMODS in .h file \n",
	 n);
      exit (0);
    }
/* Now complete the initialization of the modsum structure for this set of models*/
  comp[ncomps].modstop = nmods_tot;
  comp[ncomps].nmods = comp[ncomps].modstop - comp[ncomps].modstart;

/* Comp contains a model structure and so we need to allocate space for the 
 * w and f in this structue at this point 
 */

  /* ??? get_models Seems like one needs to set the parameters here.  Also note that there a partial
   * ??? initialization of the comp in init_models.  Seems like it should all be done in one
   * ??? place
   */


  comp[ncomps].nwaves = nwaves;
  comp[ncomps].xmod.w = calloc (nwaves, sizeof (float));
  comp[ncomps].xmod.f = calloc (nwaves, sizeof (float));


  for (n = 0; n < nwaves; n++)
    {
      comp[ncomps].xmod.w[n] = mods[comp[ncomps].modstart].w[n];
    }

  if (comp[ncomps].nmods == 0)
    {
      Error
	("get_models: No models from %s were read. Please check list of models!\n",
	 comp[ncomps].name);
      exit (0);
    }

  /* The next 3 lines set a normalization that is used by kslfit.  They are mostly
   * not relevant to python, where comp[ncomp[.min[0] refers to a normalization for 
   * a model.  I've kept this initialization for now */
  comp[ncomps].min[0] = 0;
  comp[ncomps].max[0] = 1000;

  for (m = 0; m < npars; m++)
    {
      comp[ncomps].min[m] = xmin[m];
      comp[ncomps].max[m] = xmax[m];
    }

  ncomps++;
  return (ncomps-1);		/* Return the component number, which is one less
 				than the current numer of components  */
}

/* Determine whether a model component of this varname exists and if so
 * return the component number
 */


int
model_locate (modelname)
     char *modelname;
{
  int n;
  n = 0;
  Log ("model_locate %d\n", ncomps);
  while (n < ncomps)
    {
      printf ("ModelLocate %s var %s\n", modelname, comp[n].var);
      if (strcmp (modelname, comp[n].var) == 0)
	{
	  Log ("Matched\n");
	  return (n);
	}
      n++;
    }
  return (-99);			// There is no model of this name
}

/* Get a single model model 
   This routine simple reads a model file from disk and puts the result into the structure
   onemod.  This file need not have the same wavelengths as the data nor other models.

080915	ksl	Added error to catch the case where the model being read in has more
		wavelengths than allowed, ie. more than NWAVES
080916	ksl	Changed variables to floats to save space, and made changes
		to reflect the fact that w,f in the structure are only 
		allocated when the data are read in
080917	ksl	Modified call to match with get_python_model
 */

float *last_model_w;
int last_model_nwaves = 0;

int
get_one_model (filename)
     char filename[];
{
  FILE *ptr;
  char dummy[LINELEN];
  int n, nn, nmod;
  float w[NWAVES], f[NWAVES], f_out[NWAVES];
  struct Model w_in, w_out;

  /* First find the model */

  nmod = 0;
  while (nmod <= nmods_tot && strcmp (mods[nmod].name, filename) != 0)
    nmod++;

  if (nmod > nmods_tot)
    {
      Error
	("get_one_model: Could not locate filename %s in Model struc with %d models\n",
	 filename, nmods_tot);
      exit (0);
    }
  else if (mods[nmod].nwaves > 0) {
	  /* We have already read this model*/
	  return(mods[nmod].nwaves);
  }


  /* Check if this a PYTHON model and if so read that type of file */
  if (comp[mods[nmod].ncomp].filetype == PYTHON)
    {
      get_python_model (filename);
      return (0);

    }
  /* OK we know which model we are dealing with and so we can consume the dat */

  if ((ptr = fopen (filename, "r")) == NULL)
    {
      Error ("get_one_model: Could not open filename -%s-\n", filename);
      exit (0);
    }




  n = 0;
  while (n < NWAVES && (fgets (dummy, LINELEN, ptr)) != NULL)
    {
      if ((dummy[0] != '#'))
	{
	  sscanf (dummy, "%e %e", &w[n], &f[n]);
	}
      n++;
    }


  fclose (ptr);

  if (n >= NWAVES)
    {
      Error ("get_one_model: model %s has more than %d wavelengths\n",
	     filename, NWAVES);
    }

  /* This is the point at which the model should be 
   * regridded if that is required.  
   */

  if (spec.regrid_spectra == 1)
    {
      w_in.w = w;
      w_in.f = f;
      w_in.nwaves = n;
      w_out.w = spec.w;
      w_out.nwaves = spec.nw;
      w_out.f = f_out;
      regrid (&w_in, &w_out);

      assign_wavelengths (&mods[nmod], w_out.nwaves, w_out.w);

      mods[nmod].f = calloc (w_out.nwaves, sizeof (float));
      nn = 0;
      while (nn < w_out.nwaves)
	{
	  mods[nmod].f[nn] = w_out.f[nn];
	  nn++;
	}

      return (w_out.nwaves);	/* Want to populate comp.nwaves properly */

    }

  else
    {

      /* Now allocate space for the flux and models */

      assign_wavelengths (&mods[nmod], n, w);


      /* Now we can finish up */

      mods[nmod].f = calloc (n, sizeof (float));
      nn = 0;
      while (nn < n)
	{
	  mods[nmod].f[nn] = f[nn];
	  nn++;
	}

    }

  n_xflux++;			// A diagnostic telling one how many flux arrays have been created
  fclose (ptr);
  return (n);
}

/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  model interpolates between models and places the results in 
  	xcomp->xmod  (It is important to note that the model is un-renormalized
	so any renormalization has to be done outside of this routine.)

  Description:	

  Arguments:		

  Returns:
  	//OLD - Returns nwaves in the the spectrum if a new spectrum was calculated
	0 if it was the old spectrum.

	Returns a pointer to the model in the comp structure.

  Notes:

  It one has a dense grid of models covering all of the parameter space on needs
  it's fairly straightforward to figure out how to interpolate.  If the grid is not
  dense however, or there are missing models (as can happen for example with hubeny
  stellar atmospheres where attempts to generate a complete grid fail because the
  hubeny atmospheres fail to converge) the results can be problematic.

  The solution here is to assume that the first parameter is the most important, the 
  second the next most important, etc. So if you have a grid with t and g. If
  you are sant to interpolate a model that lies half-way between t1 and t2, and g1 and g2, and assume
  the following models exist in the grid, (t1,g1) (t2,g1), (t2,g2), then (t1,g1) gets 
  a weight of 1/2 since it was the only model with a temperature of t, while t2,g1 and t2,g2
  both get weights of 1/4.
  
  The solution is "practical" because it always gives a result, and
  because it gives the an appropriate answer when the grid is filled.  It is also
  easily extensible to more than two dimensions.  
  
  But it should be noted, the approach is also dangerous, since it really
  does depend on the user having put the variables in in a sensible order.

  History:
04aug	ksl	Adapted from similar routine used by kslfit.  In this
		case the results are placed in xcomp->mod
06jul	ksl	57g -- fixed bug when trying to retrieve the best model
		when the first parameter (T) was exactly at a model value
		but this modeld would be excluded by the gravity calculation.
		Also did a little clean-up of comments
080928	ksl	Until now this routine has been named model, which was 
		extremely confusing
080929	ksl	There is a difference in the way the current version of these
		functions operates in pythona and pyfit3, which ??? needs to
		be fixed.  For kslfit, we need to preserve parameter 0 for 
		a normalization, but we really want to models to be interpolated
		without the normallization.  ??? Needs to be reconciled i
		universal version ???

 ************************************************************************/
int nmodel_error = 0;

ModelPtr
make_model (CompPtr xcomp)
{
  int j, n;
  int good_models[NMODS];	// Used to establish which models are to be included in creating output model
  double xmin[NPARS], xmax[NPARS];	// The vertices of a completely filled grid
  double weight[NMODS];		// The weights assigned to the models
  double hi, lo, delta, wtot;
  int ngood;
  double f;
  int nwaves;
  double flux[NWAVES];


/* First determine whether we already have interpolated this
exact model previously */

  n = 0;
  while (n < xcomp->npars && xcomp->xmod.par[n] == xcomp->new[n])
    {
      n++;
    }
  if (n == xcomp->npars)
    {
      return (&xcomp->xmod);	// This was the model stored in comp already
    }


  /* So we need to create a new model */
  /* First identify the models of interest */

  n = 0;
  while (n < xcomp->modstart)
    {
      weight[n] = good_models[n] = 0;
      n++;
    }
  while (n < xcomp->modstop)
    {
      weight[n] = good_models[n] = 1;
      n++;
    }
  while (n < nmods_tot)
    {
      weight[n] = good_models[n] = 0;
      n++;
    }

  /* Start with j=1 because 0 is reserved for normalization
   * but we intend to return an unrenormalized model. 
   * ??? The sum here is up to < npars because in this version
   * of the code we have inclued the normalization in the count
   * In the old kslfit version of the code we counted differently.
   * */

  for (j = 1; j < xcomp->npars; j++)
    {
      xmax[j] = xcomp->max[j];
      xmin[j] = xcomp->min[j];
      hi = BIG;
      lo = -BIG;
      for (n = xcomp->modstart; n < xcomp->modstop; n++)
	{
	  if (good_models[n])
	    {
	      delta = mods[n].par[j] - xcomp->new[j];
	      if (delta > 0.0 && delta < hi)
		{
		  xmax[j] = mods[n].par[j];
		  hi = delta;
		}
	      if (delta <= 0.0 && delta >= lo)
		{
		  xmin[j] = mods[n].par[j];
		  lo = delta;
		}
	    }
	}
      // So at this point we know what xmin[j] and xmax[j] and we 
      // need to prune good_models
      for (n = xcomp->modstart; n < xcomp->modstop; n++)
	{
	  // Next lines excludes the models which are out of range.
	  if (mods[n].par[j] > xmax[j] || mods[n].par[j] < xmin[j])
	    good_models[n] = 0;
	  // Next line modifies the weight of this model assuming a regular grid
	  // If xmax==xmin, then par[j] was outside of the range of the models and
	  // so we need to weight the remaining models fully.
	  if (good_models[n] && xmax[j] > xmin[j])
	    {
	      f = (xcomp->new[j] - xmin[j]) / (xmax[j] - xmin[j]);
	      if (mods[n].par[j] == xmax[j])
		{
		  // Then the model is at the maximum for this parameter
		  weight[n] *= f;
		}
	      else
		weight[n] *= (1. - f);

/* 57g -- If the weight given to a model is going to be zero, it needs to be
excluded from furthur consideration -- 07jul ksl */
	      if (weight[n] == 0.0)
		good_models[n] = 0;

	    }
	}
    }

  // At this point, we should have all the input models we want to include in the
  // final output weighting, as well as the relative weighting of the models.  

  wtot = 0;
  ngood = 0;
  for (n = xcomp->modstart; n < xcomp->modstop; n++)
    {
      if (good_models[n])
	{
	  wtot += weight[n];
	  ngood++;
	}
    }
  if (wtot == 0)
    {
      Error
	("model: Wtot must be greater than 0 or something is badly wrong\n");
      exit (0);
    }
  for (n = xcomp->modstart; n < xcomp->modstop; n++)
    {
      if (good_models[n])
	weight[n] /= wtot;
    }

// So now we know the absolute weighting.

  if (ngood == 0)
    {
      Error ("model: No models from %s survived pruning\n", xcomp->name);
      exit (0);
    }
  else if (ngood == 1 && nmodel_error < 20)
    {
      Error
	("model: Only one model after pruning for parameters, consider larger model grid\n");
      for (n = xcomp->modstart; n < xcomp->modstop; n++)
	{
	  if (good_models[n])
	    {
	      Error ("model: %s %8.2f %8.2f\n", mods[n].name, xcomp->new[1],
		     xcomp->new[2]);
	    }
	}
      nmodel_error++;
    }

  nwaves = xcomp->nwaves;

// Now create the spectrum 
  for (j = 0; j < nwaves; j++)
    {
      flux[j] = 0;
    }


  for (n = xcomp->modstart; n < xcomp->modstop; n++)
    {
      if (good_models[n])
	{
	  get_one_model (mods[n].name);
	  for (j = 0; j < nwaves; j++)
	    {
	      flux[j] += weight[n] * mods[n].f[j];
	    }
	}
    }



  for (j = 0; j < nwaves; j++)
    {
      xcomp->xmod.f[j] = flux[j];
    }
  for (j = 0; j < xcomp->npars; j++)
    {
      xcomp->xmod.par[j] = xcomp->new[j];
    }

  return (&xcomp->xmod);	// This is a pointer to the new model
}


/* 
 * Initialize the Model and component structures
 */

int
init_models ()
{
  int n;

  n = 0;
  while (n < NMODS)
    {
      strcpy (mods[n].name, "");
      mods[n].w = NULL;
      mods[n].f = NULL;
      mods[n].nwaves = -1;
      n++;
    }

  n_xwav = 0;
  n_xflux = 0;
  nmods_tot = 0;
  ncomps = 0;			// The number of different sets of models that have been read in

  n = 0;
  while (n < NCOMPS)
    {
      strcpy (comp[n].name, "");
      comp[n].xmod.w = NULL;
      comp[n].xmod.f = NULL;
      n++;
    }

  return (0);
}


/* This routine populates the wavelength structure if necessary
 * and assigns onemod->w
 */
int
assign_wavelengths (onemod, nwaves, waves)
     struct Model *onemod;
     int nwaves;
     float waves[];
{
  int k;
  int nn;




  k = 0;
  while (k < n_xwav)
    {
      if (nwaves == xwav[k].nwaves)
	{
	  // it's probably the same as a set of old ones but we need to check 

	  nn = 0;
	  while (nn < nwaves && xwav[k].w[nn] == waves[nn])
	    {
	      nn++;
	    }
	  if (nn == nwaves)
	    {
	      break;
	    }
	}

      k++;
    }
  if (k < n_xwav)
    {
      // The is no need to allocate new wavelengths

      onemod->w = xwav[k].w;
      onemod->nwaves = nwaves;
      return (0);
    }


  /* If we do need a new set of wavelengths we create the array and populate it */

  xwav[n_xwav].w = calloc (nwaves, sizeof (float));
  nn = 0;
  while (nn < nwaves)
    {
      xwav[n_xwav].w[nn] = waves[nn];
      nn++;
    }
  xwav[n_xwav].nwaves = nwaves;
  Log ("Allocating new wavelength array for %s\n", onemod->name);

  onemod->w = xwav[n_xwav].w;
  onemod->nwaves = nwaves;
  n_xwav++;

  /* Now we can finish up */
  return (0);
}



/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  
	parse_xy (line, xcol, ycol, w, f)

  Description:	

  	parse a character string line into discrete words and attempt
	to convert the words corresponding the the xcolth and ycolth
	word into a float

  Arguments:  (Input via .pf file)		


  Returns:

  Returns 0 on success, and a postive integer (see code) if
  one, or both of the conversions have failed.

  Notes:


  History:
  0809	ksl	Coded as part of the effort to get better model 
  		reading routines

 ************************************************************************/


int
parse_xy (line, xcol, ycol, w, f)
     char *line;
     int xcol, ycol;
     float *w, *f;
{
  char delims[] = " 	";
  char *result = NULL;
  int i, ii, xerr, yerr;
  char sentence[132];

  strcpy (sentence, line);

  *w = 0.0;
  *f = 0.0;
  xerr = 1;
  yerr = 10;

  result = strtok (sentence, delims);
  i = 1;
  while (result != NULL)
    {
      // printf ("result is %s\n", result);
      if (i == xcol)
	{
	  if ((ii = sscanf (result, "%e", w)) == 1)
	    {
	      xerr = 0;
	    }
	}
      if (i == ycol)
	{
	  if ((ii = sscanf (result, "%e", f)) == 1)
	    {
	      yerr = 0;
	    }
	}
      i++;
      result = strtok (NULL, delims);
    }

  return (xerr + yerr);
}
