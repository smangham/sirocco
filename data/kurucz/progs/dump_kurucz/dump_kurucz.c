
/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:

  Description:	This program takes the spectrum files on Bob Kurucz's web
        site and breaks it up into individual ascii files which can then
        be used in model fitting.  

  Arguments:		

  Returns:

  Notes:
	I suspect that one needs to check for each file one downloads 
	what the number of wavelengths is.

  History:
	00dec	ksl	Coded
	

 ************************************************************************/
#include	<math.h>
#include	<stdio.h>
#include	<strings.h>
#define LINELENGTH 132
#define SPECMAX    1221		// Number of wavelengths in the grid //

main (argc, argv)
     int argc;
     char *argv[];
{

  FILE *fptr, *optr, *sptr, *fopen ();
  char infile[80], outfile[80];
  char line[LINELENGTH];
  char firstword[LINELENGTH];
  float x[8];
  float wave[1221], flux[1221], cont[1221];
  float teff, gravity;

  int stage;
  int nwords;
  int ii, n;


  stage = 0;
/* Parse the command line */
  if (argc == 2)
    {
      strcpy (infile, argv[1]);
    }
  else
    {
      printf ("Usage: PROGRAM filename\n");
      exit (0);
    }

  sptr = fopen ("kuruczlist", "w");

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
      if (strncmp (firstword, "END", 3) == 0)
	{
	  stage = 1;
	  ii = 0;
	}
      else if (strncmp (firstword, "TEFF", 4) == 0)
	{
	  stage = 2;
//        printf ("%s", line);
	  sscanf (line, "%*s %f %*s %f", &teff, &gravity);
	  printf ("TEFF %f %f\n", teff, gravity);
	  strcpy (outfile, "");
	  sprintf (outfile, "fp00k2t%05.0fg%03.1f.dat", teff, gravity);
	  optr = fopen (outfile, "w");
	  fprintf (sptr, "%-20s  %10.0f %5.1f\n", outfile, teff, gravity);
	  ii = 0;
	}
      else if (stage == 1)
	{			// Then this is a wavelength record 
//        printf ("%s", line);
	  nwords =
	    sscanf (line, "%10f%10f%10f%10f%10f%10f%10f%10f", &x[0], &x[1],
		    &x[2], &x[3], &x[4], &x[5], &x[6], &x[7]);
//        printf ("%d Wavelength %f %f %f %f %f %f %f %f\n", nwords, x[0],
//                x[1], x[2], x[3], x[4], x[5], x[6], x[7]);
	  for (n = 0; n < nwords; n++)
	    {
	      wave[ii] = x[n]*10.;  // Convert from nanometers to Angstroms
	      ii++;
	    }
	}
      else if (stage == 2)
	{			// This is a flux record
//        printf ("%s", line);
	  nwords =
	    sscanf (line, "%10e%10e%10e%10e%10e%10e%10e%10e", &x[0], &x[1],
		    &x[2], &x[3], &x[4], &x[5], &x[6], &x[7]);
//        printf ("%d Flux %e %e %e %e %e %e %e %e\n", nwords, x[0], x[1],
//                x[2], x[3], x[4], x[5], x[6], x[7]);
	  for (n = 0; n < nwords; n++)
	    {
	      flux[ii] = x[n];
	      ii++;
	    }
	  if (ii == SPECMAX)
	    {
	      ii = 0;
	      stage = 3;

	    }
	}
      else if (stage == 3)
	{			// this is a continuum record
//        printf ("%s", line);
	  nwords =
	    sscanf (line, "%10e%10e%10e%10e%10e%10e%10e%10e", &x[0], &x[1],
		    &x[2], &x[3], &x[4], &x[5], &x[6], &x[7]);
//       printf ("%d Flux %e %e %e %e %e %e %e %e\n", nwords, x[0], x[1],
//                x[2], x[3], x[4], x[5], x[6], x[7]);
	  for (n = 0; n < nwords; n++)
	    {
	      cont[ii] = x[n];
	      ii++;
	    }
	  if (ii == SPECMAX)
	    {
	      for (n = 0; n < SPECMAX; n++)
		{
		  flux[n] *= 2.997e18 / (wave[n] * wave[n]);	// Convert to f_lambda
		  cont[n] *= 2.997e18 / (wave[n] * wave[n]);	// Convert to f_lambda
		  //printf ("%f %9.3e %9.3e\n", wave[n], flux[n], cont[n]);
		  //fprintf(optr,"%.1f %9.3e %9.3e\n", wave[n], flux[n], cont[n]);
		  fprintf (optr, "%.1f %9.3e\n", wave[n], flux[n]);
		}
	      fclose (optr);
	    }
	}

    }


  fclose (sptr);
}
