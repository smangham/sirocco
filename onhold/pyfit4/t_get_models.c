


/**************************************************************************
                    Space Telescope Science Institute


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

t_get_models: t_get_models.o
	gcc t_get_models.o  $(LDFLAGS) -o t_get_models
#	mv $@ $(BIN)

 

  History:
2004	ksl	Coded as better.c

 ************************************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


#include "log.h"
#include "models.h"

#define LINELENGTH 132

int
main (argc, argv)
     int argc;
     char *argv[];
{

  char parfile[LINELENGTH];
  char infile[LINELENGTH];
  int spectype;
  double xpar[10];
  char fileout[LINELENGTH];
  int n, m,nn;


  strcpy (fileout, "foo.dat");

/* Parse the command line */
  if (argc == 1)
    strcpy (parfile, "t_get_models.pf");
  else if (argc == 2)
    {
      strcpy (parfile, argv[1]);
      if ((strstr (parfile, ".pf")) == NULL)
	strcat (parfile, ".pf");
    }
  else
    {
      printf ("Usage: t_get_models [file.pf]\n");
      exit (0);
    }

/* Initialize variables that will be read by rdpar */

  strcpy (infile, "disk06/python.ls");

/* Get data from a parameter file */
 
  get_data("spectrum.dat");

  opar (parfile);


  rdstr ("filename", infile);

  cpar (parfile);

/* End of input section */

  spectype = get_models (infile, "disk06", 1);

  printf ("returned spectype %d\n", spectype);

/* This should not really do anything because we have already read these 
 * models the first time
 */

  spectype = get_models (infile, "disk06", 1);

  printf ("returned spectype %d\n", spectype);


  /* Now print some details of this model */

  printf ("There are %d components\n", ncomps);

  for (n = 0; n < ncomps; n++)
    {
      printf ("var       %s\n", comp[n].var);
      printf ("name      %s\n", comp[n].name);
      printf ("filetype  %d\n", comp[n].filetype);
      printf ("nmods     %d    modstart  %d modstop %d\n", comp[n].nmods,
	      comp[n].modstart, comp[n].modstop);
      printf ("npars     %d\n", comp[n].npars);
      for (m = 0; m < comp[n].npars;m++)
	{
	  printf ("min    %8.2g     max   %8.2g\n", comp[n].min[m], comp[n].max[m]);
	}
    }

/* Now get a model */

  nn=0;

      for (m = 0; m < comp[0].npars;m++)
	{
	  xpar[m]=0.5*(comp[0].min[m]+ comp[0].max[m]);
	}

  rdint ("comp.to.calculate", &nn);

  while (xpar[0] > 0.0)
    {

      for (m = 0; m < comp[0].npars;m++)
	{
	  printf ("min %f max %f\n", comp[n].min[m], comp[n].max[m]);
	  rddoub ("xpar", &xpar[m]);
	}

      model ("disk06", xpar);

      put_model (fileout, spectype);

      printf ("wave arrays instantiated  %d\n", n_xwav);
      printf ("flux arrays instantiated  %d\n", n_xflux);
    }


  model ("OH.NO", xpar);
/* OK now end program and go on to the next important task. */
  exit (0);

}


int
put_model (filename, spectype)
     char *filename;
     int spectype;
{
  FILE *fopen (), *fptr;
  int n;

  /* Open the input file.  Exit if it is not opened successfully */

  if ((fptr = fopen (filename, "w")) == NULL)
    {
      printf ("Failed to open file %s\n", filename);
      exit (0);
    }

  n = 0;
  while (n < comp[spectype].nwaves)
    {
      fprintf (fptr, "%10.3f %10.4e\n", comp[spectype].xmod.w[n],
	       comp[spectype].xmod.f[n]);
      n++;
    }

  fclose (fptr);
  return (0);

}
