

/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:

  Description:	This program creates an sm file.  It
    draws a solid WD and a solid disk as viewed at an
    arbitrary inclination angle theta.  It then reads
    a file of the following format
         x  z  value
    of the type that can be produced by the python routines

    and creates a file with the values xp yp

  Arguments:		

  Returns:

  Notes:

  History:
	01jul04	ksl	Modified so that multiple inclination could be handled
			in the same go	

 ************************************************************************/
#include	<math.h>
#include	<stdio.h>
#include	<strings.h>
#include        "fitsio.h"

#define		NPTS 1000
#define 	LINELENGTH 	132
#define         RADIAN		57.29578
#define 	IMAP 		128

struct basis
{
  double a[3][3];
}
base;

struct points
{
  double x[3];
}
pts[NPTS];

  short map[128][128];

main (argc, argv)
     int argc;
     char *argv[];
{

  FILE *fptr, *optr, *fopen ();
  char root[LINELENGTH], pffile[LINELENGTH], answer[LINELENGTH];
  char logfile[LINELENGTH];
  char line[LINELENGTH], infile[LINELENGTH];
  char ofile[LINELENGTH];
  char ffile[LINELENGTH];
  double rwd, rdisk, theta, thetamin, thetamax, dtheta;
  int n, nang;
  double lmn[3], x[3], y[3], z[3];
  double s,size;
  double length ();
  int ii,jj;


  if (argc == 2)
    {
      get_root (root, argv[1]);
    }
  else
    {
      strcpy (answer, "");
      printf ("Parameter file: ");
      fgets (answer, LINELENGTH, stdin);
      get_root (root, answer);
    }

  strcpy (pffile, root);
  strcat (pffile, ".pf");
  printf ("Parameter file name %s\n", pffile);
  strcpy (logfile, root);
  strcat (logfile, ".log");

  Log_init (logfile);

  opar (pffile);

  rddoub ("rstar(cm)", &rwd);
  rddoub ("disk.radmax(cm)", &rdisk);
  rddoub ("inclination_min(deg)", &thetamin);
  rddoub ("inclination_max(deg)", &thetamax);
  rdint ("num_inclinations", &nang);
  rdstr ("points.file", infile);
  rddoub ("size",&size);

  dtheta = (thetamax - thetamin) / nang;


  for (n = 0; n <= nang; n++)
    {
      theta = thetamin + dtheta * n;

/* Calculate the direction cosines to the observer */
      lmn[0] = sin (theta / RADIAN);
      lmn[1] = 0;
      lmn[2] = cos (theta / RADIAN);

      y[0] = 0.;
      y[1] = 1.;
      y[2] = 0.;

/* Create the transformation matrix so that any vector can be projected
toward the observer */

      create_basis (lmn, y, &base);

/* Open an output file */
      strcpy (ofile, "");
      sprintf (ofile, "%s.%02d.out", root, n);
      optr = fopen (ofile, "w");
      strcpy (ffile, "");
      sprintf (ffile, "%s.%02d.fit", root, n);

/* Create the disk */

      disk (rwd, rdisk, root, n);



/* Zero the array */
   for(ii=0;ii<IMAP;ii++){
	for(jj=0;jj<IMAP;jj++){
map[ii][jj]=0;
}
}
      
/* Open the input file.  Exit if it is not opened successfully */
      if ((fptr = fopen (infile, "r")) == NULL)
	{
	  printf ("Failed to open file %s\n", infile);
	  exit (0);
	}

/* Read and print the input file */
      while (fgets (line, LINELENGTH, fptr) != NULL)
	{
	  printf ("%s\n", line);
	  sscanf (line, "%lf %lf %lf", &x[0], &x[1], &x[2]);

	  if (x[2] > 0)
	    {
	      project_to (&base, x, y);
	      fprintf (optr, "%.2e %.2e\n", y[1], y[2]);
              ii=y[1]/(2.*size)*IMAP+IMAP/2;
              jj=y[2]/(2.*size)*IMAP+IMAP/2;
              if(ii>IMAP-1)ii=IMAP-1;
              if(ii<0)ii=0;
              if(jj>IMAP-1)jj=IMAP-1;
              if(jj<0)ii=0;
              map[jj][ii]+=1;
	    }
	  else
	    {
	      s = (-x[2]) / lmn[2];
	      vmove (x, lmn, s, y);
	      printf ("test %e %e %e\n", y[0], y[1], y[2]);
	      if (length (y) > rdisk)
		{
		  project_to (&base, x, y);
		  fprintf (optr, "%.2e %.2e\n", y[1], y[2]);
              ii=y[1]/(2.*size)*IMAP+IMAP/2;
              jj=y[2]/(2.*size)*IMAP+IMAP/2;
              if(ii>IMAP-1)ii=IMAP-1;
              if(ii<0)ii=0;
              if(jj>IMAP-1)jj=IMAP-1;
              if(jj<0)ii=0;
//              map[ii][jj]+=1;
              map[jj][ii]+=1;
		}
	    }


	}

      write_fits(ffile);
      fclose (optr);
    }
}

#define NPHI 100
#define NR   10

int
disk (r1, r2, root, m)
     double r1, r2;
     char root[];
     int m;
{
  double dphi, dr, phi, r;
  double x[3], xout[3];
  char ofile[LINELENGTH];

  FILE *dptr, *fopen ();
  int i, n;

/* Open an output file */
  strcpy (ofile, "");
  sprintf (ofile, "%s.%02d.disk", root, m);
printf("%s",ofile);
  dptr = fopen (ofile, "w");

  dphi = 360. / RADIAN / NPHI;
  dr = (r2 - r1) / NR;
  for (i = 0; i <= NR; i++)
    {
      r = r1 + dr * i;
      for (n = 0; n <= NPHI; n++)
	{
	  phi = dphi * n;
	  x[0] = r * cos (phi);
	  x[1] = r * sin (phi);
	  x[2] = 0;
	  project_to (&base, x, xout);
	  if (n == 0)
	    {
	      fprintf (dptr, "relocate %.2e %.2e\n", xout[1], xout[2]);
	    }
	  else
	    {
	      fprintf (dptr, "draw %.2e %.2e\n", xout[1], xout[2]);
	    }

	}

    }
}


write_fits(name)
char name[];
{
fitsfile *fptr;
int status, ii, jj;
long fpixel=1, naxis=2, nelements;
long naxes[2]={IMAP,IMAP};

status=0;
fits_create_file(&fptr,name,&status);

fits_create_img(fptr, SHORT_IMG, naxis, naxes, &status);

nelements=naxes[0]*naxes[1];

fits_write_img(fptr, TSHORT, fpixel, nelements, map[0], &status);

fits_close_file(fptr,&status);

fits_report_error(stderr, status);

return(0);
}
