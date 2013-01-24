#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


#include "log.h"




#include "py_compare.h"
#include "templates.h"

int line2words_nwords;
char line2words_words[MAXWORD][LINELENGTH];



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
get_python_spec (filename, nspec, one_spec)
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
