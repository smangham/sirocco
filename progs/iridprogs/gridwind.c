

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"

#include "python.h"

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	This file contains routines for handling the interaction
	between the GridPtr and WindPtr array
  
 Arguments:		

	
 Returns:
  
Description:	

Notes:


History:
 
**************************************************************/



/***********************************************************
                                       Space Telescope Science Institute

 Synopsis: create_maps (ichoice)
  
 Arguments:		
     int ichoice	!0 (or true) then the size of the plasma structure is 
			just given by the cells with non-zero volume
			0 (or false) then the size of the palma structure
			is the same as the wind structue

	
 Returns:
  
Description:	

Notes:


History:
	06may	ksl	Created as part of the attempt to reduce the
			overall size of the structures
 
**************************************************************/


int
create_maps (ichoice)
     int ichoice;
{
  int i, j;
  j = 0;
  if (ichoice)
    {				// map to reduce space
      for (i = 0; i < NDIM2; i++)
	{
	  wmain[i].nwind = i;
	  if (wmain[i].vol > 0)
	    {
	      wmain[i].nplasma = j;
	      plasmamain[j].nplasma = j;
	      plasmamain[j].nwind = i;
	      j++;
	    }
	  else
	    {
	      wmain[i].nplasma = NPLASMA;
	    }
	}
      if (j != NPLASMA)
	{
	  Error
	    ("create_maps: Problems with matching cells -- Expected %d Got %d\n",
	     NPLASMA, j);
	  exit (0);
	}

    }

  else				// one plama cell for each wind cell
    {
      for (i = 0; i < NDIM2; i++)
	{
	  wmain[i].nwind = i;
	  wmain[i].nplasma = i;
	  plasmamain[i].nplasma = i;
	  plasmamain[i].nwind = i;

	}

    }

  plasmamain[NPLASMA].nplasma = NPLASMA;
  plasmamain[NPLASMA].nwind = -1;
  return (0);
}



/***********************************************************
                                       Space Telescope Science Institute

 Synopsis: calloc_wind (nelem)
  
 Arguments:		

	
 Returns:
  
 Description:	

 Notes:


 History:
	06may	ksl	57a -- Coded
 
**************************************************************/


int
calloc_wind (nelem)
     int nelem;
{

  wmain = (WindPtr) calloc (sizeof (wind_dummy), nelem + 1);

  if (wmain == NULL)
    {
      Error
	("There is a problem in allocating memory for the wind structure\n");
      exit (0);
    }
  else
    {
      Log
	("Allocated %10d bytes for each of %5d elements of      w totaling %10.1f Mb\n",
	 sizeof (wind_dummy), nelem, 1.e-6 * nelem * sizeof (wind_dummy));
    }

  return (0);
}



/***********************************************************
                                       Space Telescope Science Institute

 Synopsis: calloc_plasma (nelem)
  
 Arguments:		

	
 Returns:
  
Description:	

Notes:


History:
 
**************************************************************/


int
calloc_plasma (nelem)
     int nelem;
{

// Allocate one extra element to store data where there is no volume

  plasmamain = (PlasmaPtr) calloc (sizeof (plasma_dummy), (nelem + 1));
  geo.nplasma = nelem;

  if (plasmamain == NULL)
    {
      Error
	("There is a problem in allocating memory for the plasma structure\n");
      exit (0);
    }
  else
    {
      Log
	("Allocated %10d bytes for each of %5d elements of plasma totaling %10.1f Mb \n",
	 sizeof (plasma_dummy), (nelem + 1),
	 1.e-6 * (nelem + 1) * sizeof (plasma_dummy));
    }

/* Now allocate space for storing photon frequencies -- 57h*/
  photstoremain =
    (PhotStorePtr) calloc (sizeof (photon_store_dummy), (nelem + 1));

  if (photstoremain == NULL)
    {
      Error
	("There is a problem in allocating memory for the photonstore structure\n");
      exit (0);
    }
  else
    {
      Log
	("Allocated %10d bytes for each of %5d elements of photonstore totaling %10.1f Mb \n",
	 sizeof (photon_store_dummy), (nelem + 1),
	 1.e-6 * (nelem + 1) * sizeof (photon_store_dummy));
    }

  return (0);
}


int
check_plasma (xplasma, message)
     PlasmaPtr xplasma;
     char message[];
{
  if (xplasma->nplasma == NPLASMA)
    {
      Error ("check_plasma -- %s \n", message);
      return (1);		// true
    }
  else
    return (0);
}



/***********************************************************
               Space Telescope Science Institute

 Synopsis: calloc_macro (nelem)
  
 Arguments:		

	
 Returns:
  
Description:	

Notes:

060805 -- ksl -- Begin notes on changing the macro structures

At present, we allocate macromain if there are any macrolevels,
but the size of macromain is determined by NLEVELS_MACRO and
NBBJUMPS or NBFJUMPS.  These are not set dynamically, and so
one should consider how to do this, since for the current values
these are set to 500 for NLEVELS_MACRO, 20 for NBBJUMPS, and 
NBFJUMPS, which means they are roughly 1000 per level.  The
first thing to do about this is to reduce them to values that
are more like those in our datafiles.

If further work is needed there are several possibilties, 
each of which has some disadvantages relative to the current
appraoch.  The advantage of the current approch is that each
plasma cell has a separate macro structure element and that
these are addressed rather logically, e.g
	macromain[n].jbar[macro_level][nbf]
It is pretty easy to see what is going on in arrays.  If one
wants to dynamically allocate say the number of macro_levels,
then one will create a new structure called for example 
macro2 to contain jbar, jbar_old etc, and then to map to the
appropriate elements of it  from the macro_array.  There are 
several ways to do this.

One would be to follow the approach used in other portions of
the atomic structures to add t add the place to start in the
macro2 structure as an integer.  So for example, if there
were 35 macrolevels read in, one could simply add a number
that told one where to start in the macro2 array, much as
is done in looking for firstion and lastion for an element.

In that case, the addressing would look somelike
	macro2[macromain[n].start_level]].jbar[nbf]

This is not crazy since we know how many levels have
been read in, and so the start_levels will increase
by the number of levels.  An advantage of this approach
is that you can dump the structues to disk, and read
them back with another program and eveything will be
properly indexed.

Another appraoch that would work would be to create
pointers instead of indexes.  This may be faster. In
this case the addresing would look something like.
	macro[n].lev_ptr[nlevel]->jbar[nbf]
where lev_ptr was defined as *Macro2Ptr within the
macro structure.

A final approach would be to dynamically resize the
way the macro array itselve is addressed within each
routine. This can also be done.  One simply creates a that
parallels the macro structue within each routine, e.g

struct macro xmacro[][nlevel]

and addrsses things as

        xmacro[n][nlev].jbar[nbf]

This also looks fairly straightforward, in c, and reads pretty
esaily.  


At present, I have not done anything about this

End notes on resizing -- 060805


History:
	06jun	ksl	Coded
	0608	ksl	57h -- Recoded to eliminate the creation
			of the array when there are no macro_atoms.
			Previously the array had been created but
			not written out.
 
**************************************************************/


int
calloc_macro (nelem)
     int nelem;
{


  if (nlevels_macro == 0)
    {
      Log ("Allocated no space for macro since nlevels_macro==0\n");
      return (0);
    }
// Allocate one extra element to store data where there is no volume

  macromain = (MacroPtr) calloc (sizeof (macro_dummy), (nelem + 1));
  geo.nplasma = nelem;

  if (macromain == NULL)
    {
      Error
	("There is a problem in allocating memory for the macro structure\n");
      exit (0);
    }
  else if (nlevels_macro > 0)
    {
      Log
	("Allocated %10d bytes for each of %5d elements of macro  totaling %10.1f Mb \n",
	 sizeof (macro_dummy), (nelem + 1),
	 1.e-6 * (nelem + 1) * sizeof (macro_dummy));
    }
  else
    {
      Log ("Allocated no space for macro since nlevels_macro==0\n");
    }

  return (0);
}
