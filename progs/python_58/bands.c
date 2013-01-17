
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************
                                       Space Telescope Science Institute

Synopsis:

   
Arguments:		
Returns:
 
 
Description:	

		
Notes:

History:
**************************************************************/


/* Next bits are for banding of the frequency intervals, which is an attempt to improve the
spacing of frequencies to assure there are enough ionizing photons to stabilize photoabsorption.
01dec ksl (python40)
*/

// Actual structures are in python.h.  Here for reference only.
//#define NBANDS 10
//struct bands
//{
//  double f1,f2;
//  double min_fraction;
//  double nat_fraction;                // The fraction of the accepted luminosity in this band
//  double used_fraction;
//  double f;                   //The "luminosity" within a band
//  double weight;
//  int nphot;
//}
//band[NBANDS];

//int nbands;           // Actual number of bands in use
/* End of banding parameters */




/***********************************************************
                                       Space Telescope Science Institute

Synopsis:
This is the routine that initializes the bands.  It is hardwired
at present 

   
Arguments:		
	mode	0	Use temperature to define a single band
		1	Use f1 and f2 to define a single band
		2	Use t,f1 and f2, and hardwired bands
			to define multiple bands
Returns:
 
 
Description:	

		
Notes:

History:
	02jul	ksl	Adapted from an similar routine used to set
			up bands in photon_gen.  This may eventually
			replace that routine
**************************************************************/


int
init_bands (t, f1, f2, mode, band)
     double t;			// A temperature which can be used to set abolute limits on the bands
     double f1, f2;		// frequency limits that can overide any other limits
     int mode;			// A switch used for determining how the bands are to be populated
     struct xbands *band;

{
  if (mode == 0)
    {
      band->nbands = 1;
//              band->f1[0] = C / 10000e-8;     /*10000 A */
      band->f1[0] = BOLTZMANN * t / H * 0.05;
      band->f2[0] = BOLTZMANN * t / H * 20.;
      band->min_fraction[0] = 1.0;
    }
  else if (mode == 1)
    {
      band->nbands = 1;
      band->f1[0] = f1;
      band->f2[0] = f2;
      band->min_fraction[0] = 1.0;
    }
  else if (mode == 2)
    {
      band->nbands = 4;
      band->f1[0] = f1;
      band->f2[0] = band->f1[1] = 13.599 / HEV;
      band->f2[1] = band->f1[2] = 24.588 / HEV;
      band->f2[2] = band->f1[3] = 54.418 / HEV;
      band->f2[3] = f2;
      band->min_fraction[0] = 0;
      band->min_fraction[1] = 0.1;
      band->min_fraction[2] = 0.1;
      band->min_fraction[3] = 0.1;
      if (f1 > band->f2[0])
	{
	  Error ("init_bands: f1 (%e) > 13.599/HEV)\n", f1);
	  exit (0);
	}
      if (f2 < band->f2[2])
	{
	  Error ("init_bands: f2 (%e) < 54.418/HEV)\n", f2);
	  exit (0);
	}

    }
  else
    {
      Error ("Init bands: Unknown mode %d\n", mode);
      mytrap ();
      //exit(0);
    }


  return (0);
}
