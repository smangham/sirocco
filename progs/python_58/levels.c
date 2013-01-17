
/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:
	levels (w, mode) calculates the fractional occupation numbers of 
	the various levels of atomic configurations as designated in
	the atomic data files

  Description:	

	mode	0	LTE with t_r
		1	LTE with t_e
		2	Non-LTE (reduced by weighted BB)

  Arguments:		

  Returns:

  Notes:

	At present (01oct), levels is called through ion_abundances in 
	python.

  History:
	01sep23	ksl	Began work
	01oct10	ksl	Modified so modes matched python ionization modes
			more exactly.
	01dec03	ksl	Modified to simplify so modes match those of
			nebular concentrations
	01dec12	ksl	Modified to react to changes which split "nlte"
			and "lte" levels.  Levels is explicitly for so-
			called "nlte" levels, which are tracked in the
			Wind structure
	04Apr   SS      If statement added to avoid this routine changing
                        macro atom level populations.
        04May   SS      The if statment added above is modified for the case
                        of all "simple" ions.
	06may	ksl	57+ -- Modified to make use of plasma structue

 ************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

int
levels (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;
{
  double t, weight;
  int n, m, mm;
  int mm_ground;
  double kt;
  double z;

  int ireturn;

  if (mode == 0)		// LTE with t_r 
    {
      t = xplasma->t_r;
      weight = 1;
    }
  else if (mode == 1)		// LTE with t_e
    {
      t = xplasma->t_e;
      weight = 1;
    }
  else if (mode == 2)		// non_LTE with t_r and weights
    {
      t = xplasma->t_r;
      weight = xplasma->w;
    }
  else
    {
      Error ("levels: Could not calculate levels for mode %d\n", mode);
      exit (0);
    }

/* Next calculation should be almost identical to that contained in the routine
partition.  An error in one place implies one in the other.  It seems like
that this routine should be subsumed into that one ???? */
  kt = BOLTZMANN * t;

  for (n = 0; n < nions; n++)
    {
      if (ion[n].nlte > 0)
	/* Extra if statement added to prevent changing levden of macro atom populations (SS, Apr04) */
	{
	  if (ion[n].macro_info == 0 || geo.macro_ioniz_mode == 0)

	    {			// Then calculate levels for this ion
	      z = xplasma->partition[n];	/* N.B. partition functions will most likely have
						   been calculated from "lte" levels, at least for now ?? */
	      mm = ion[n].first_nlte_level;
	      mm_ground = mm;	//store the ground state index - allow for gs energy neq 0 (SS)
	      xplasma->levden[mm] = config[mm].g / z;	// Assumes first level is ground state
	      for (m = 1; m < ion[n].nlte; m++)
		{
		  mm++;
		  xplasma->levden[mm] =
		    weight * config[mm].g *
		    exp ((-config[mm].ex + config[mm_ground].ex) / kt) / z;
		}
	    }
	}
    }
  return (0);

  return (ireturn);
}
