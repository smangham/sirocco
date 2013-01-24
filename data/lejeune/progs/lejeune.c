


/**************************************************************************
                    Center for Astrophysical Sciences
                        Johns Hopkins University


  Synopsis:   Read a set of lejeune models and write them out in
	in individual files.  Also write out a file that gives the
        temperatures and gravities of the models

  Description:	

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:

Look at http://vizier.cfa.harvard.edu/viz-bin/Cat?J/A%2bAS/125/229#sRM2.1
for an explanation

  This is roughly what the Makefile should look like. Uncomment the mv 
  statement to store binary in ~/bin

CC = gcc
CFLAGS = -c -g -I$$HOME/include -Wall
LDFLAGS= -L$$HOME/lib -lm -lkpar
BIN = $$HOME/bin

lejeune: lejeune.o
	gcc lejeune.o  $(LDFLAGS) -o lejeune
#	mv $@ $(BIN)

 

  History:
1986	ksl	Coded as better.c

 ************************************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#define NWAVES 1221


#include "log.h"

#define LINELENGTH 132

int
main (argc, argv)
     int argc;
     char *argv[];
{

  FILE *fptr, *fopen (), *optr, *qptr;
  char parfile[LINELENGTH];
  char infile[LINELENGTH];
  char outfile[LINELENGTH];
  char filelist[LINELENGTH];
  char line[LINELENGTH];
  char firstword[LINELENGTH];
  float wave[NWAVES], f[NWAVES];
  int i, j, nwords;
  float x[8],z;


/* Parse the command line */
  if (argc == 1)
    strcpy (parfile, "lejeune.pf");
  else if (argc == 2)
    {
      strcpy (parfile, argv[1]);
      if ((strstr (parfile, ".pf")) == NULL)
	strcat (parfile, ".pf");
    }
  else
    {
      printf ("Usage: lejeune [file.pf]\n");
      exit (0);
    }

/* Initialize variables that will be read by rdpar */

  strcpy (infile, "file.ls");

/* Get data from a parameter file */


  opar (parfile);


  rdstr ("filename", infile);

  cpar (parfile);

/* End of input section */


  if ((qptr = fopen ("lejeune.ls", "w")) == NULL)
    {
      printf ("Failed to open file %s\n", infile);
      exit (0);
    }
  /* Open the input file.  Exit if it is not opened successfully */
  if ((fptr = fopen (infile, "r")) == NULL)
    {
      printf ("Failed to open file %s\n", infile);
      exit (0);
    }

/* Read and print the input file */
  i = 0;
  while (i < 1221)
    {
      if (fgets (line, LINELENGTH, fptr) == NULL)
	{
	  Error ("Problem reading wavelengths\n");
	}
      nwords =
	sscanf (line, "%e %e %e %e %e %e %e %e", &x[0], &x[1], &x[2], &x[3],
		&x[4], &x[5], &x[6], &x[7]);
      for (j = 0; j < nwords; j++)
	{
	  wave[i + j] = x[j];
	}
      i += j;
    }

a:if (fgets (line, LINELENGTH, fptr) == NULL)
    {
      Error ("Problem reading model name or end of file\n");
    }
  nwords =
    sscanf (line, "%e %e %e %e %e %e %e %e", &x[0], &x[1], &x[2], &x[3],
	    &x[4], &x[5], &x[6], &x[7]);

  if (nwords != 6)
    {
      Error ("Failed to find t and g properly\n");
      exit (0);
    }
  sprintf (outfile, "lej_t%05.0fg%04.1f.txt", x[1], x[2]);
  fprintf (qptr, "%-30s %5.0f %5.1f\n", outfile, x[1], x[2]);

  optr=fopen(outfile,"w");

/* Now read the fluxes */

  i = 0;
  while (i < 1221)
    {
      if (fgets (line, LINELENGTH, fptr) == NULL)
	{
	  Error ("Problem reading fluxes\n");
	}
      nwords =
	sscanf (line, "%e %e %e %e %e %e %e %e", &x[0], &x[1], &x[2], &x[3],
		&x[4], &x[5], &x[6], &x[7]);
      for (j = 0; j < nwords; j++)
	{
	  f[i + j] = x[j];
	}
      i += j;
    }
/* Now write everything out */

// See the doucmentation for the conversion */

   for(j=0;j<i;j++){
	z=(0.4*f[j]*2.997925e17/(wave[j]*wave[j]));
  fprintf (optr, "%10.1f %10.4e\n", (wave[j]*10.0),z);
} 

fclose(optr);

  goto a;

  exit (0);



}
