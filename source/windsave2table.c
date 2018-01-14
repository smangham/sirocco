
/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	windsave2table writes key variables in a wind save file to an astropy table    
		as calculated by python.  This is the main routine.

Arguments:		

	py_wind  windsave_root



Returns:
 
Description:	
	

	
Notes:

	The main difficulty with this program is that one needs to be consistent
	regarding the size of the arrays that one stuffs the variables into.  
	As now written, if one wants to access a variable in wmain, one needs to
	include and offset, generally called nstart.


History:
	150428	ksl	Adapted from routines in py_wind.c
	160216	ksl	Resolved issues with multiple domains
    1706    ksl Refactored so that windsave2table_subs could be
                called from within python

**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
#include "python.h"



int
main (argc, argv)
     int argc;
     char *argv[];
{



  char root[LINELENGTH], input[LINELENGTH];
  char outputfile[LINELENGTH];
  char windsavefile[LINELENGTH];
  char parameter_file[LINELENGTH];
  int create_master_table (), create_ion_table ();
  int do_windsave2table ();


  // py_wind uses rdpar, but only in an interactive mode. As a result 
  // there is no associated .pf file

  strcpy (parameter_file, "NONE");

  /* Next command stops Debug statements printing out in py_wind */
  Log_set_verbosity (3);

  if (argc == 1)
  {
    printf ("Root for wind file :");
    fgets (input, LINELENGTH, stdin);
    get_root (root, input);
  }
  else
  {
    strcpy (input, argv[argc - 1]);
    get_root (root, input);
  }


  printf ("Reading data from file %s\n", root);

  /* Now create the names of all the files which will be written */

  strcpy (windsavefile, root);
  strcpy (outputfile, root);

  strcat (windsavefile, ".wind_save");
  strcat (outputfile, ".txt");


/* Read in the wind file */

  if (wind_read (windsavefile) < 0)
  {
    Error ("py_wind: Could not open %s", windsavefile);
    exit (0);
  }


  printf ("Read wind_file %s\n", windsavefile);

  get_atomic_data (geo.atomic_filename);

  printf ("Read Atomic data from %s\n", geo.atomic_filename);

  do_windsave2table (root);
}
