


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:


  Arguments:		

  Returns:

  Notes:

  History:
	

 ************************************************************************/
#include	<math.h>
#include	<stdlib.h>
#include	<stdio.h>
#include	<strings.h>
#include	<string.h>
#include         "models.h"
#include	"log.h"
#define    BIG 1e32




/* Get the data 
 *
 
History
	0708	ksl	Removed the lines that prevented reading of the first two lines
			of data, which I think was leftover from when I was using
			both kslfit and specfit
	0809	ksl	Simplified from the existing version in kslfit.  I have removed
			edges in the structue, since we really only need that for 
			regridding, and I want to be able to do it more generally.


 * */
int
get_data (filename)
     char filename[];
{
  FILE *ptr;
  char dummy[LINELEN];
  int n;

  if ((ptr = fopen (filename, "r")) == NULL)
    {
      Error ("get_data: Could not open filename %s\n", filename);
      exit (0);
    }
  strcpy (spec.name, filename);

  n = 0;
  while (n < NWAVES && (fgets (dummy, LINELEN, ptr)) != NULL)
    {
      if ((dummy[0] != '#'))
	{
	  sscanf (dummy, "%e %e %e", &spec.w[n], &spec.f[n], &spec.e[n]);
	  n++;
	}
    }
  spec.nw = n;

  fclose (ptr);

  /* This is a request that models which are read in be regridded to
   * this wavelength scale.
   */

  spec.regrid_spectra = 1;

  return (0);
}

/* A structure to hold information for regrid */

struct ReGrid
{
//OLD  float *w;
//OLD  float *ww;
  int n;
  int nn;

  float w_edge[NWAVES];
  float ww_edge[NWAVES];

} RRGrid;



/* Note that we are using spectrum structures here.  
 * spec_in is  */
int
regrid (mod_in, mod_out)
     struct Model *mod_in;
     struct Model *mod_out;
{
  float *w, *ww;
  float *f, *ff;
  float *w_edge, *ww_edge;
  float x, y, wbot, wtop;
  int n, nwaves, nnwaves;
  int i, imin, imax;

  /* All of the next variables have good and final
   * values, ecept for mod_out.f
   */

  w = mod_in->w;
  ww = mod_out->w;
  f = mod_in->f;
  ff = mod_out->f;
  nwaves = mod_in->nwaves;
  nnwaves = mod_out->nwaves;
  w_edge = RRGrid.w_edge;
  ww_edge = RRGrid.ww_edge;

  /* Check whether we need to recreate edges */

//OLD  if (w != RRGrid.w)
//OLD    {
//OLD      Log ("Recreating input spectrum edges\n");
      w_edge[0] = w[0] - 0.5 * (w[1] - w[0]);
      for (n = 1; n < nwaves; n++)
	{
	  w_edge[n] = 0.5 * (w[n + 1] + w[n]);
	}
      w_edge[n] = w[n - 1] + 0.5 * (w[n - 1] - w[n - 2]);
//OLD      RRGrid.w = w;
//OLD    }


//OLD  if (ww != RRGrid.ww)
//OLD    {
//OLD      Log ("Recreating output spectrum edges\n");

      ww_edge[0] = ww[0] - 0.5 * (ww[1] - ww[0]);
      for (n = 1; n < nnwaves; n++)
	{
	  ww_edge[n] = 0.5 * (ww[n + 1] + ww[n]);
	}
      ww_edge[n] = ww[n - 1] + 0.5 * (ww[n - 1] - ww[n - 2]);
//OLD      RRGrid.ww = ww;
//OLD    }


  /* At this point we know the edges are correct and we can follow the 
   * basic prcdure we have followed previously.
   */

  imin = 0;
  for (n = 0; n < nnwaves; n++)
    {
      ff[n] = 0;
      /* identify the first interval that contributs to the flux
       */

      while (w_edge[imin] <= ww_edge[n] && imin < nwaves)
	imin++;
      imin--;
      imax = imin;
      while (w_edge[imax] < ww_edge[n + 1] && imax < nwaves)
	imax++;
      imax--;

      // Now calculate the average flux in the wavelength interval

      x = y = 0;
      for (i = imin; i <= imax; i++)
	{
	  wbot = w_edge[i];
	  if (wbot < ww_edge[n])
	    wbot = ww_edge[n];
	  wtop = w_edge[i + 1];
	  if (wtop > ww_edge[n + 1])
	    wtop = ww_edge[n + 1];
	  x += f[i] * (wtop - wbot);
	  y += (wtop - wbot);
	}
      if (y > 0)
	ff[n] = x / y;
    }

  return (0);


}
