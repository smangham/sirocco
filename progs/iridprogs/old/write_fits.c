
/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:

  Description:	This program 

  Arguments:		

  Returns:

  Notes:

  History:
	01aug	ksl	Replaced ndim with xdim so ndim could be used in
			python.h
	

 ************************************************************************/

#include        <math.h>
#include        <stdio.h>
#include        <strings.h>
#include 	<string.h>
#include        "fitsio.h"


int
write_fits (name, array, xdim)
     char name[];
     float array[][];
     int xdim;
{
  fitsfile *fptr;
  int status;
  long fpixel = 1, naxis = 2, nelements;
  long naxes[2];
  char filename[132];

  naxes[0] = naxes[1] = xdim;

  strcpy (filename, name);
  strcat (filename, ".fit");

  status = 0;
  fits_create_file (&fptr, filename, &status);

  fits_create_img (fptr, FLOAT_IMG, naxis, naxes, &status);

  nelements = naxes[0] * naxes[1];

  fits_write_img (fptr, TFLOAT, fpixel, nelements, array, &status);

  fits_close_file (fptr, &status);

  fits_report_error (stderr, status);

  return (0);
}
