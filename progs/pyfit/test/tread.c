


/**************************************************************************
                        Space Telescope Science Institute


  Synopsis:  model_compare 

  Description:	

  	This routine compare two python models with one another 

	Usage 

	model_compare [-m firstmodel nspec -m secondmodel nspec] pffile  

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

#include <stdio.h>
#include <string.h>

#define MAXWORD 80

int line2words_nwords = 0;
char line2words_words[MAXWORD][LINELENGTH];


typedef struct Spec
{
  char name[LINELENGTH];
  float ang;
  float w[NWAVES];
  float f[NWAVES];
  int good[NWAVES];
  int nw, ng;
} spec_dummy, *SpecPtr;




int
main (argc, argv)
     int argc;
     char *argv[];
{
  char filename[LINELENGTH];
  FILE *fopen (), *optr;
  int n, nspec;
  struct Spec model;
  int get_python_model ();

  strcpy (filename, "polar.spec");
  nspec = 9;


  printf ("reading %s", filename);

  n=get_python_model (filename, nspec, &model);

  printf ("spectrum %d wavelengths %d\n", n,model.nw);

  if (n==0){
	  printf("Quitting\n");
	  exit(0);
  }


  optr = fopen ("foo.out", "w");

  for (n = 0; n < model.nw; n++)
    {
      fprintf (optr, "%8.2f %8.3e\n", model.w[n], model.f[n]);
    }



  exit (0);
}









/**************************************************************************
                        Space Telescope Science Institute


  Synopsis:   Read a specific column of a python model

  Description:	

  Arguments:  

  Returns:


  Notes:


  History:
09may	ksl	Adapted from similar routin in pyfit3, but this returns
		a single spectrum


 ************************************************************************/


int
get_python_model (filename, nspec, one_spec)
     char filename[];
     int nspec;
     SpecPtr one_spec;
{
  FILE *ptr;
  char dummy[LINELENGTH];
  int mspec, nwaves;

  if ((ptr = fopen (filename, "r")) == NULL)
    {
      Error ("Could not open filename %s\n", filename);
      return (0);
    }


  nwaves = 0;

  mspec = 5;			// This is the offset for reading a model

  while (fgets (dummy, LINELENGTH, ptr) != NULL)
    {
      line2words (dummy);

      /* Process the comment lines */
      if (dummy[0] == '#')
	{
	  if (strncmp (line2words_words[1], "Freq", 4) == 0)
	    {
	      mspec = 2;
	      while (mspec < line2words_nwords)
		{
		  if (strncmp (line2words_words[mspec], "Ang", 3) == 0)
		    {
		      break;
		    }

		  mspec = mspec + 1;
		}
	    }
	  if (mspec == line2words_nwords)
	    {
	      Error ("get_python_model: Could not parse Freq line\n");
	      exit (0);
	    }
	}
      else
	      /* Process the data lines */
	{
		if (mspec + nspec-2 >= line2words_nwords){
	      		Error ("get_python_model: spectrum %d does not exist\n",nspec);
			return(0);
		}

	  sscanf (line2words_words[1], "%f", &one_spec->w[nwaves]);
	  sscanf (line2words_words[mspec + nspec-2], "%e",
		  &one_spec->f[nwaves]);
	  nwaves++;

	}
    }

  one_spec->nw = nwaves;

  fclose (ptr);

  return (mspec+nspec-2);
}




#define DELIM   " "
int
line2words (string)
     char *string;
{
  int ntokens = 0;


  char *sep = strtok (string, DELIM);

  while (sep != NULL)
    {
      strncpy (line2words_words[ntokens++], sep, LINELENGTH);
      sep = strtok (NULL, DELIM);
    }

  line2words_nwords = ntokens;

  return (ntokens);

}
