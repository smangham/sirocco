


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  This is a program that prepares kurucz files for use by
	python and other programs.  The basic purpose is to assure that
	all kurucz models have the same set of wavelengths and that the
        that the units are ergs/cm**2/s/A.  	

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

ascii_expand.c: ascii_expand.o
	gcc ascii_expand.c.o  $(LDFLAGS) -o ascii_expand
#	mv $@ $(BIN)

 

  History:
04aug	ksl	Coded as part of python52 effort

 ************************************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


#include "log.h"

#define LINELENGTH 132
#define ARRAYSIZ	10000

int
main (argc, argv)
     int argc;
     char *argv[];
{

  FILE *fptr, *fopen ();
  char parfile[LINELENGTH];
  char infile[LINELENGTH];
  char wavefile[LINELENGTH];
  char line[LINELENGTH];
  char firstword[LINELENGTH];
  double wo[ARRAYSIZ], w[ARRAYSIZ];
  double flux[ARRAYSIZ];
  int no[ARRAYSIZ], n[ARRAYSIZ];
  int no_pts, n_pts;
  int raw_kurucz(),writekurucz();


/* Parse the command line */
  if (argc == 1)
    strcpy (parfile, "ascii_expand.pf");
  else if (argc == 2)
    {
      strcpy (parfile, argv[1]);
      if ((strstr (parfile, ".pf")) == NULL)
	strcat (parfile, ".pf");
    }
  else
    {
      printf ("Usage: ascii_expand.c [file.pf]\n");
      exit (0);
    }

/* Initialize variables that will be read by rdpar */

  strcpy (wavefile, "fp00t50000g50k2.dat");
  strcpy (infile, "files.ls");

/* Get data from a parameter file */


  opar (parfile);


  rdstr ("wavelength.file", wavefile);
  rdstr ("list.of.files.to.expand", infile);

  cpar (parfile);

/* End of input section */


/* Open the file that contais all the wavelenghts and read it */

/*Read a file with all the wavelengths */
  no_pts = raw_kurucz ("fp00t50000g50k2.dat", no, wo, flux);
  printf ("No of points: %d\n", no_pts);

  /* Open the input file.  Exit if it is not opened successfully */
  if ((fptr = fopen (infile, "r")) == NULL)
    {
      printf ("Failed to open file %s\n", infile);
      exit (0);
    }

/* Read and print the input file */
  while (fgets (line, LINELENGTH, fptr) != NULL)
    {
      sscanf (line, "%s", firstword);
      printf ("%s\n", firstword);
      n_pts = raw_kurucz (firstword, n, w, flux);
      printf ("No of points: %d\n", n_pts);

      writekurucz (firstword, no_pts, no, wo, n_pts, n, w, flux);

    }


/* OK now end program and go on to the next important task. */
  exit (0);

}


#define NSKIP 3
#define WAVEMAX 750000 // 750,000 nm == 750 microns
// Note that the read statement needs improving to to longward of this

int
raw_kurucz (filename, n, wave, flux)
     char *filename;
     int n[];
     double wave[], flux[];

{

  FILE *fopen (), *ptr;
  char *fgets ();

  char name[LINELENGTH];
  char line[LINELENGTH];
  int i;

  strcpy (name, filename);

  if ((ptr = fopen (name, "r")) == NULL)
    {
      printf ("Failed to open file %s\n", name);
      exit (0);
    }

  for (i = 0; i < NSKIP; i++)
    fgets (line, LINELENGTH, ptr);
  i = 0;
  while (fgets (line, LINELENGTH, ptr) != NULL)
    {
      sscanf (line, "%*s %4d %lf %*f %lf", &n[i], &wave[i], &flux[i]);
      wave[i] *= 10;
      i++;
      if (wave[i - 1] > WAVEMAX)
	return (i);
    }

  return (i);
}



int writekurucz (filename, no_pts, no, wo, n_pts, n, wave, flux)
char filename[];
int no_pts, n_pts;
int no[], n[];
double wo[];
double wave[], flux[];
{
  FILE *optr, *fopen ();
  char outfile[LINELENGTH];
  double f[ARRAYSIZ];
  int i,j;
  int kurucz_change_units();

/* Initialize f to zero and then stuff it with
the appropriate variables from flux */
  for (i = 0; i < no_pts; i++)
    {
      f[i] = 0.0;
    }

  for (i = 0; i < n_pts; i++)
    {
      j = 0;
      while (no[j] != n[i])
	j++;
      f[j] = flux[i];
    }

  i=strlen(filename);
  printf("%d %s\n",i,filename);

  strcpy(outfile,"");
  strncat (outfile, filename, i - 4);	// Assumes original extensions are .dat
  strncat (outfile, ".txt",4);  // New extensions are .txt
  printf("%s\n",outfile);

  kurucz_change_units(wo,f,no_pts);

  optr = fopen (outfile, "w");
  for (i = 0; i < no_pts; i++)
    {
      fprintf (optr,"%10.2f %10.3e\n", wo[i], f[i]);
    }


return(0);
}


/* Convert from ergs/cm**2/Hz to ergs/cm**2/Angstrom */


#define C                               2.997925e10

int
kurucz_change_units (wavelength, flux,npts)
     double *wavelength;
     double *flux;
{
  int i;

  for (i = 0; i < npts; i++)
    flux[i] *= (C * 1e8 / (wavelength[i] * wavelength[i]));

  return (0);

}

