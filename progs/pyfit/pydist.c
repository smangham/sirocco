


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  
  	pydist tries to give one a sense of how many models in a grid
	approximate the data 

  Description:	

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:

	060825--This might have been more straightforward to have written
	in python.  Reading the file over and over again is rather odd.

 

  History:
2004	ksl	Coded as better.c
06aug	ksl	Tried to see what this program did and improve it
070708	ksl	Updated to read information about the columns, if
		it exists.
090120	ksl	Small updates to correct for new column format, and
		to add a logfile

 ************************************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


#include "log.h"

#define LINELENGTH 256
#define NMAX  100
#define NCOLS 30

int
main (argc, argv)
     int argc;
     char *argv[];
{

  FILE *fptr, *optr, *fopen ();
  char parfile[LINELENGTH];
  char infile[LINELENGTH];
  char outfile[LINELENGTH];
  char diagfile[LINELENGTH];
  char line[LINELENGTH];
  char format[LINELENGTH];
  char root[LINELENGTH];
  float xx[NMAX], yy[NMAX];
  int ncol, nx;
  int i, nlines, n;
  float ratio, x, cutoff;
  int col_chi, colmin, colmax;
  int ncols, colno;
  char colname[NCOLS][LINELENGTH];
/* Parse the command line */
  if (argc == 1)
    strcpy (parfile, "pydist.pf");
  else if (argc == 2)
    {
      strcpy (parfile, argv[1]);
      if ((strstr (parfile, ".pf")) == NULL)
	strcat (parfile, ".pf");
    }
  else
    {
      printf ("Usage: pydist [file.pf]\n");
      exit (0);
    }

/* Initialize variables that will be read by rdpar */

  strcpy (root, "rootname.of.results.file");
  ncol = 10;
  col_chi = 1;
  colmin = 9;
  colmax = 17;
  cutoff = 1.3;

/* Get data from a parameter file */


  opar (parfile);


  rdstr ("input_file_root", root);

/* Begin section to get an idea of what is in the input files*/

  strcpy (infile, root);
  strcat (infile, ".results");
  strcpy (outfile, root);
  strcat (outfile, ".best");

  /* Open a log file */

  strcpy(diagfile,root);
  strcat(diagfile,".log");
        Log_init (diagfile);

/* End of creating filenames */
  if ((optr = fopen (outfile, "w")) == NULL)
    {
      Log ("Failed to open file %s\n", outfile);
      exit (0);
    }

/* Process the comments, including getting colmin and colmax */
  /* Open the input file.  Exit if it is not opened successfully */
  if ((fptr = fopen (infile, "r")) == NULL)
    {
      Log ("Failed to open file %s\n", infile);
      exit (0);
    }

  /* 070707 - this section has been updated to read more information 
   * about what the .results file contains
   */

  for (i = 0; i < NCOLS; i++)
    {
      strcpy (colname[i], "Unknown");
    }

  Log ("Extracting column information from %s\n", infile);
  ncols = 1;
  while (fgets (line, LINELENGTH, fptr) != NULL)
    {

      if (sscanf (line, "%*s Col: %d %s", &colno, colname[ncols]) == 2)
	{
	  if (ncols != colno)
	    {
	      Error ("Col no problem.  Check file!");
	      exit (0);
	    }
	  ncols++;
	}
      /* Next lines if for an older version of pyfit3 */
      else if (sscanf (line, "%*s colmin %d colmax %d", &colmin, &colmax) ==
	       2)
	break;
      /* This is for the current 0901 version of pyfit3 */
      else if (sscanf (line, "%*s Columns %d - %d",&colmin,&colmax)==2)
	      break;
//      else  {
//            printf("%s\n",line);
//      }
    }
  fclose (fptr);

/* Now describe what is being done for the Reader */

  if (ncols > 1)
    {
      for (i = 1; i < ncols; i++)
	{
	  Log ("Col: %d %s\n", i, colname[i]);
	}
    }
  else
    {
	    Error("Apparently .results file contained no column info! Continuing.\n");
    }

/* Now get the rest of the information you need to run the program */

  rdint ("col_chi(1.is.overall.best)", &col_chi);	// column containing the normalized chi*2 you want
  rdflo ("cutoff.ratio", &cutoff);	// Maximum value to keep

  cpar (parfile);

/* Finished getting inputs */

  Log ("Using chi2 from col %2d : %s\n", col_chi, colname[col_chi]);
  Log ("Keeeping modesl within %f of mimimum\n", cutoff);

/* OK now process everything by repeatedly reading the file */
  for (ncol = colmin; ncol <= colmax; ncol++)
    {

      /* Create the format line */
      strcpy (format, "");
      for (i = 1; i < col_chi; i++)
	{
	  strcat (format, "%*s ");
	}
      strcat (format, "%f ");
      for (i = col_chi; i < ncol - 1; i++)
	strcat (format, "%*s ");
      strcat (format, "%e");
      Log ("Col %4d : Parameter %4d  %s\n", ncol, ncol - colmin + 1,
	   colname[ncol]);
      for (n = 0; n < NMAX; n++)
	yy[n] = 0;
      /* Open the input file.  Exit if it is not opened successfully */
      if ((fptr = fopen (infile, "r")) == NULL)
	{
	  Log ("Failed to open file %s\n", infile);
	  exit (0);
	}

/* Read and the input file and accumulate statistics*/
      nlines = 0;
      nx = 0;
      while (fgets (line, LINELENGTH, fptr) != NULL)
	{
	  if (line[0] != '#')
	    {
	      sscanf (line, format, &ratio, &x);
	      n = 0;
	      while (n < nx && x != xx[n])
		n++;
	      if (n == nx)
		{		//It's a new x
		  xx[n] = x;
		  nx++;
		}
	      if (ratio < cutoff)
		{
		  yy[n]++;
		  /* Only write the output file once */
		  if (ncol == colmin)
		    fprintf (optr, "%s", line);
		}
	      nlines++;
	    }
	}


/* OK now simple print the statistics */

      for (n = 0; n < nx; n++)
	{
	  Log ("%10.2e  %3.0f\n", xx[n], yy[n]);
	}
      fclose (fptr);
    }

  exit (0);
}
