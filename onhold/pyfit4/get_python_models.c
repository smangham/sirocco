/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:
  	get_python_models contains a set of routine to read grids of
	models that were made with the radiative transfer program
	python.  

  Description:	These are a generic set of routines to read a grid or grids
	of models  The grid can be any reasonable number of dimension

  Arguments:		

  Returns:

  Notes:
	DO NOT INCLUDE ANYTHING THAT IS NOT GENERIC HERE!  This includes
	anything specific to python.  This is important because this
	set of routines needs to be usable outside of the envirnment
	of python.


  History:
0809	ksl	Coded.  These routines were partially derived from those
		in the current version of pyfit3
	

 ************************************************************************/
#include	<math.h>
#include	<stdio.h>
#include        <stdlib.h>
#include	<strings.h>
#include	<string.h>
#include	"log.h"
#include         "models.h"
#include	"model_prototypes.h"
#define    	BIG 1e32

#define 	LOCATE_ERR 	-99
/* Read a standard python output file */

#define MAX_ANG	6

/* 
 
Notes - This routine is only intended to be called after get_models and so some of
the initialization is not necessary.  In particular the file should already be
read.

080919	ksl Modfified to include resampling
 */
int
get_python_models (modellist, modelname, ichoice)
     char modellist[];		// filename containing location and associated parameters of models
     char modelname[];		// Descriptor for the model
     int ichoice;		// 0 means defer reading the models until later; 1 means read the models now
{
  char line[LINELEN];
  float ang[MAX_ANG];
  int nang, nangs;
//  int ncomp;

  char *dummy;
  int n, m, nxpar;
  float xpar[NPARS], xmin[NPARS], xmax[NPARS];
  int nw, nwaves;
  int npars;
  char word[LINELEN], filetype[LINELEN];
  char word2[LINELEN];
  int ix, iy;
  int nline;
  int ang_in_deg;


  Log ("get_python_models: Assuming %s is a specidal Python file\n",
       modellist);

  n = ncomps;			// We want to add a new component to the end of the components


  /* Now initialize the model summary structure */
  strcpy (comp[ncomps].name, modellist);
  strcpy (comp[ncomps].var, modelname);
  comp[ncomps].modstart = nmods_tot;
  comp[ncomps].filetype = PYTHON;

  npars = -1;			// Don'd yet know how many pars there are for this model type

/*Now retrieve the line that gives the angles */
  get_comment_line ("Angles", line);
  nangs = sscanf (line, "Angles %e %e %e %e %e %e",
		  &ang[0], &ang[1], &ang[2], &ang[3], &ang[4], &ang[5]);
  if (nangs == 0)
    {
      Error
	("get_python_models: Could not find Angles line in Python file %s\n",
	 modellist);
      exit (0);
    }
  /* We are now ready to populate it all. */
/* Beginning of copied section */



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

	  nxpar =
	    sscanf (dummy, "%s %f %f %f %f %f %f %f %f %f ",
		    word, &xpar[1], &xpar[2], &xpar[3], &xpar[4],
		    &xpar[5], &xpar[6], &xpar[7], &xpar[8], &xpar[9]);

	  xpar[0]=1.0; // Assure Xpar[0] is not undefined



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

	  /* We have the information we need at this point but for python files we are going to create 
	   * more than one line of spectra
	   */

	  nang = 0;

	  while (nang < nangs)
	    {
	      strcpy (mods[n].name, word);	// There are multiple spectra and we need to copy this each time
	      /* Make the filename unique since we search on this to find a model */

	      ang_in_deg = ang[nang];

	      sprintf (word2, ".Ang%02d", ang_in_deg);
	      strcat (mods[n].name, word2);

	      /* Add in the angle to the parameter file */

	      for (m = 0; m < npars; m++)
		{
		  mods[n].par[m] = xpar[m];
		}

	      mods[n].par[npars] = xpar[npars] = ang[nang];
	      comp[ncomps].npars = npars + 1;;

	      for (m = 0; m < comp[ncomps].npars; m++)
		{
		  if (xpar[m] > xmax[m])
		    xmax[m] = xpar[m];
		  if (xpar[m] < xmin[m])
		    xmin[m] = xpar[m];
		}

	      for (m = comp[ncomps].npars; m < NPARS; m++)
		mods[n].par[m] = -99;

	      mods[n].ix = ix;
	      mods[n].iy = iy;
	      mods[n].ncomp = ncomps;	// Associate the model with a component


	      if ((n % 1000) == 0)
		printf ("Model n %d %s\n", n, mods[n].name);

	      n++;
	      nang++;
	    }

	  /* It is necessary to read at least one model to be able to create the wavelength and flux grid in comp.
	   * Or we would need to check when we tried to populate if the first time.  However unless ichoice is 0
	   * we also read all of the models at this point  */

	  nmods_tot = n;
	  if (n == nang || ichoice)
	    nwaves = get_python_model (mods[n - 1].name);

	}

      nline++;
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

  for (m = 0; m < comp[ncomps].npars; m++)
    {
      comp[ncomps].min[m] = xmin[m];
      comp[ncomps].max[m] = xmax[m];
    }

  ncomps++;
  return (ncomps-1);		/* Return the component number */

  /* End of copied section */
}





int
get_python_model (modelname)
     char modelname[];
{
  FILE *ptr;
  char dummy[LINELEN];
  char filename[LINELEN];
  int i, n, m, nscanned;
  float ang[MAX_ANG], wavelength[NWAVES], f[MAX_ANG][NWAVES];
  float wav[NWAVES], flux[NWAVES];	// properly ordered wavelengths and fluxes
  float f_out[NWAVES];
  struct Model w_in, w_out;


  int mspec, nwaves;

  /* The model name may have a Ang designator attached, we need
   * to search on this and strip it off if neccessary
   */
  get_python_filename (modelname, filename);

  if ((ptr = fopen (filename, "r")) == NULL)
    {
      Error ("get_python_model: Could not open filename %s\n", filename);
      exit (0);
    }


// Get the title line

  for (m = 0; m < MAX_ANG; m++)
    ang[m] = 999;

  /* Kluge to handle both old and new style spectral lines.  Originally the #Freq line was the first line.  
   * Today 07, we use that line tor record the versio of Python run.  
   * 080917 - There would be several better ways to do this.  One would be simply to look for comment liens
   * at the beginning and parse them and the second would be to actuall construct the scan string on the
   * fly given MAX_ANG*/
  fgets (dummy, LINELEN, ptr);
  mspec =
    sscanf (dummy,
	    "#Freq.        Lambda  Emitted  Star+BL  Disk     %*s       %*s       %*s       Ang:%f  Ang:%f  Ang:%f  Ang:%f  Ang:%f Ang:%f",
	    &ang[0], &ang[1], &ang[2], &ang[3], &ang[4], &ang[5]);

  if (mspec == 0)
    {				// It's a new style line
      fgets (dummy, LINELEN, ptr);
      mspec = sscanf (dummy,
		      "# Freq.        %*s %*s %*s %*s %*s       %*s       %*s       Ang:%f  Ang:%f  Ang:%f  Ang:%f  Ang:%f Ang:%f",
		      &ang[0], &ang[1], &ang[2], &ang[3], &ang[4], &ang[5]);

      if (mspec == 0)
	{
	  Error ("get_python_model: Could not scan header:%s", dummy);
	  exit (0);
	}


    }

  /* Now read the rest of the file */

  n = 0;
  while (n < NWAVES && (fgets (dummy, LINELEN, ptr)) != NULL)
    {
      if ((dummy[0] != '#'))
	{
	  nscanned =
	    sscanf (dummy, "%*e %e %*e %*e %*e %*e %*e %*e %e %e %e %e %e %e",
		    &wavelength[n], &f[0][n], &f[1][n], &f[2][n], &f[3][n],
		    &f[4][n], &f[5][n]);
	  if (nscanned != mspec + 1)
	    Error ("Problem in file %s (%d scanned/%d expected): %s",
		   filename, nscanned, mspec, dummy);
	}
      n++;
    }


  fclose (ptr);

  /* At this point, we have all of the data and it is internally held.  The wavelength order is however
   * backwards
   */


  nwaves = n;

  /* Get an array in the correct wavelength order.  Note that one could have been clever and done the reversal
   * during the sscanf stage above, but this should not be dominating our performance.
   */

  for (n = 0; n < nwaves; n++)
    {
      wav[n] = wavelength[nwaves - 1 - n];
    }

  i = 0;
  for (m = 0; m < mspec; m++)
    {
      /* The next step is slightly tricky because at this point 
       * the models have unique names, with .Ang70 or whatever
       * attached.  So must just compare the "filename" portion
       */

      while (i < nmods_tot
	     && strncmp (mods[i].name, filename, strlen (filename)) != 0)
	i++;
      if (i == nmods_tot)
	{
	  Error
	    ("get_python_model: Could not locate filename %s and spectrum %d in Model struc\n",
	     filename, m);
	  exit (0);
	}


      /* Reverse the order of the wavelengths.  Note that one could be more clever
       * about this if there was no possiblity of having to regrid, which requres
       * both the wavelengths and the flux in the right order. */

      for (n = 0; n < nwaves; n++)
	{
	  flux[n] = f[m][nwaves - 1 - n];
	}

//

      if (spec.regrid_spectra == 1)
	{
	  w_in.w = wav;
	  w_in.f = flux;
	  w_in.nwaves = nwaves;
	  w_out.w = spec.w;
	  w_out.nwaves = spec.nw;
	  w_out.f = f_out;
	  regrid (&w_in, &w_out);

	  assign_wavelengths (&mods[i], w_out.nwaves, w_out.w);

	  mods[i].f = calloc (w_out.nwaves, sizeof (float));

	  n = 0;
	  while (n < w_out.nwaves)
	    {
	      mods[i].f[n] = w_out.f[n];
	      n++;
	    }


	}

      else
	{



	  //
	  assign_wavelengths (&mods[i], nwaves, wav);

	  mods[i].f = calloc (nwaves, sizeof (float));

	  n = 0;
	  while (n < nwaves)
	    {
	      mods[i].f[n] = flux[n];
	      n++;
	    }
	}

      i++;			// Since we don't want to continue filling in the same model.

    }




  fclose (ptr);

  if (spec.regrid_spectra == 1)
    {
      return (w_out.nwaves);
    }

  return (nwaves);		/* This is needed to properly populate the comp value of nwaves */
}






/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

     get_python_filename (modelname,filename) 

  Description:	

     The get_models routines use the filename as a way
     to lacate an specfic model.  Hence the names must
     be unique.  However, for python files, there are
     multiple spectra in each actual file.  Sometimes
     we defer actually populating a model spectrum
     until it is actually needed.

     This routine parses the unique model name and
     recovers the filename which is to be read.

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:


  History:
080925	ksl	Part of effort to update model reading routines

 ************************************************************************/

int
get_python_filename (modelname, filename)
     char *modelname, *filename;
{
  char *tail;
  int k;

  tail = strstr (modelname, ".Ang");
  if (tail == NULL)
    {
      strcpy (filename, modelname);
      printf ("No extension\n");
      return (0);
    }

  k = strlen (modelname) - strlen (tail);

  strcpy (filename, "");
  strncat (filename, modelname, k);
  return (1);
}
