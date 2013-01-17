/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	Python is a program designed to simulate the transfer of radiation in a wind.  It uses the
	Sobolev approximation.  It models a wind as a biconical flow.     
	
	This is the "main" routine for Python.  It's basic purpose is to gather the input variables 
	and to control the flow of the program
 
Arguments:		

Returns:
 
Description:	
	Python is far too complicated to describe.  Basically it simulates the radiative transfer
	of photons through the wind of a cataclysmic variable or a star.  The kinematic formulation for
	the CV wind is due to Schlossman and Vitello while that of the star is due to Castor & Larmors. 
	
	 Radiation from an optically thick disk, the WD star, a boundary layer and the wind itself
	 may be included
	
	There are 4 basic portions to the program which are easy to see in the main program.
	
	1. A data gathering stage
	
	2. A calculation of the state of ionization of the wind.
	
	3. A calculation of a detailed spectrum in which the ionization is held fixed.
	
		
Notes:
	The program has been designed and tested both on Suns and MacIntoshes (with Symentec's
	C compiler).  Some of its peculiarities are due to memory limitations and relatively small 
	stack sizes on the Mac.  When compiling, check to see that the global varible MAC is
	set properly in python.h.

History:
 	97jan   ksl	Coding on python began.
 	97jul	ksl	Added capability model the spectrum as a function of orbital phase.
 	97aug	ksl	Added boundary layer as additional source of radiation
 	97aug	ksl	Added capability to extract spectra with a specific number of scatters
 	97nov	ksl	Changed error and diagnostic logging
 	98july	ksl	Added radiation from the wind itself
 	98nov	ksl	Added spherical winds and greatly modified I/O  
	99jan	ksl	Added ability to enter parameter file (w/out) root on command line
	99jul	ksl	Begin eliminating MAC dependencies
	00sep	ksl	Began work on adding a new coronal possibility
	01mar	ksl	Added a knigge wind
        01sept	ksl	Addied a thierry wind
	01dec	ksl	Incoporated improved atomic data reading codes, and
			made a variety of other major changes to the way in
			which many of the subroutines operate.
	01dec	ksl	Revised way in which old windfile option works.  It
			was not working at all originally, because of a
			change I made in April.  The new method is a bit dangerous
			because it does change a good bit of geo, since all of the
			input data is read in.  In principle, his can create an inconsistency
			between some of the info in geo and some in the wind array.
			But it keeps the wind the same, which is what was desired.  
	01dec	ksl	Added a general capability to python(40) so that during the
			the ionization cycle photons do not necessarily have identical
			weights.
	02jan	ksl	Moved specific inputs for KWD wind to knigge.c
	02feb	ksl	Added anisotropic wind scattering
	02feb	ksl	41.1 -- Allowed user to chose between anisotropic and
			isotropic scattering in the wind.
	02feb	ksl	41.2 -- Added aditional options for choosing what photons
			to count
	02apr	ksl	43.1 -- Added additional radiative transfer options.  Modified main
			so that a spectrum is printed out after each cycle
	02jul	ksl	43.7 -- Considerable changes associated with photoionization
			and recombination were made.  The initial attemps to allow
			for using a more accurate detailed balance approach were
			included in the ionizatio routines.
	03dec	ksl	Added back some timeing ability
        04Mar   SS      Minor changes to set switch for macro atom method.
                        geo.line_mode=6 switches geo.rt_mode=2 (macro method)
                        and then geo.ine_mode back to =3. (SS)
        04Apr   SS      Minor change to set a switch for geo.line_mode=7 
                        which is the same as 6 but with anisotropic scattering
        04May   SS      geo.line_mode=8 added: this mode is the same as mode 6
                        BUT will treat all ions as simple - it is for 
                        test purposes only.
	04Jun	ksl	Added new disk + wind model for ysos.  In process, moved
			the variables describing the wind mdots to subroutines
	04Aug	ksl	52 -- Modified input variables to allow a variety of input models
			all ascii-based.
	04Aug	ksl	52 -- Modified input variables to allow for various disk illumination
			models
	04Aug	ksl	52 -- Modified to allow for vertically extended disks. Note that some of
			the input variables were rearranged in the process.
        04Aug   SS      Minor modifications to cope with finite disk thickness calculations.
	04dec	ksl	54a -- Minor mod to calculation of fraction of luminosity that hits
			disk so that correct values are printed out in situation where there
			are multiple subcycles.
	04dec	ksl	54e -- Minor mod to record master atomic file in geo, and, for the 
			case of fixed concentration to store the file where the concentrations
			are contaned.
	05apr	ksl	56 -- Added new variable geo.binary_system to avoid binary system
			questions in cases where one was modelling a single object.  Note
			that I did not remove seemingly superfluous quesions regarding
			phase because one can still extract photons in various azimuthal
			directions, and if we ever move to 3d situations one might want
			to deal with non-axially symmetric systems.
	05jul	ksl	56d -- Updated to include new windcone definitions.
	05jul	ksl	56d -- Completed basic work to allow vertically extended disk
	06may	ksl	57g -- Began work on modifying structures to lower memory requirments
			At empty wind cells have mostly been eliminated, and macroatoms
			have been split off into a separate structure so this is not 
			used.
	06jul	ksl	57h -- Have contined work to speed program up as much as possible
			by giving more control over the parameters that govern the speed
			such as SMAX_FRAC.
	06oct	ksl	58 -- This version is a cleaned up version of py57ib, which differs
			from 57h primarily in bb.c and pdf.c.   In addition, this version
			makes it possible to use the SV velocity law in what is otherwise
			a KWD description of the wind.  

 	
 	Look in Readme.c for more text concerning the early history of the program.

**************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"

#define LINELENGTH 132

#include "python.h"

int
main (argc, argv)
     int argc;
     char *argv[];
{
  WindPtr w;
  PhotPtr p;

  int icycles, wcycles, pcycles;
  double freqmin, freqmax;
  double swavemin, swavemax, renorm;
  int n, nangles, photons_per_cycle, subcycles;
  int iwind;

/* Next three lines have variables that should be a structure, or possibly we
should allocate the space for the spectra to avoid all this nonsense.  02feb ksl */

  double angle[NSPEC - 3], phase[NSPEC - 3];
  int scat_select[NSPEC - 3], top_bot_select[NSPEC - 3];
  double rho_select[NSPEC - 3], z_select[NSPEC - 3], az_select[NSPEC - 3],
    r_select[NSPEC - 3];

  char yesno[20];
  int select_extract, select_spectype;
  double tmax;
  char root[LINELENGTH], input[LINELENGTH], wspecfile[LINELENGTH],
    specfile[LINELENGTH], diskfile[LINELENGTH];
  char windradfile[LINELENGTH], windsavefile[LINELENGTH];
  char photfile[LINELENGTH], diagfile[LINELENGTH],
    old_windsavefile[LINELENGTH];
  char atomic_filename[LINELENGTH];
  char dummy[LINELENGTH];
  double xbl;
  double n_ioniz, lum_ioniz;
  int j, nn;
  double zz, zzz, zze, ztot;
  int icheck;
  FILE *fopen (), *qptr;

  int get_models ();		// Note: Needed because get_models cannot be included in templates.h
  char model_list[LINELENGTH];
  double theat;
  int disk_illum;
  int istandard, keep_photoabs;


  printf
    ("This program simulates radiative transfer in a (biconical) CV or (spherical) stellar wind\n");

  /* Determine whether input data should be read from a file or from the terminal.  */
  if (argc == 2)
    {
      strcpy (dummy, argv[1]);
    }
  else
    {
      printf ("Input file (interactive=stdin):");
      fgets (dummy, LINELENGTH, stdin);
    }
  get_root (root, dummy);

  if (strncmp (root, "dummy", 5) == 0)
    {
      printf
	("Proceeding to create rdpar file in dummy.pf, but will not run prog\n");
    }
  else if (strncmp (root, "stdin", 5) == 0 || strncmp (root, "rdpar", 5) == 0
	   || root[0] == ' ' || strlen (root) == 0)
    {
      strcpy (root, "mod");
      printf
	("Proceeding in interactive mode\n Output files will have rootname mod\n");
    }
  else
    {
      strcpy (input, root);
      strcat (input, ".pf");

      opar (input);
      printf ("Reading data from file %s\n", input);

    }

  /* Now create the names of all the files which will be written.  Note that some files
     have the same root as the input file, while others have a generic name of python.
     This is intended so that files which you really want to keep have unique names, while
     those which are for short-term diagnostics are overwritten.  ksl 97aug. */

  strcpy (basename, root);	//56d -- ksl --Added so filenames could be created by routines as necessary

  strcpy (wspecfile, root);
  strcpy (specfile, root);
  strcpy (windradfile, "python");
  strcpy (windsavefile, root);
  strcpy (photfile, "python");
  strcpy (diagfile, root);
  strcpy (diskfile, root);

  strcat (wspecfile, ".spec_tot");
  strcat (specfile, ".spec");
  strcat (windradfile, ".wind_rad");
  strcat (windsavefile, ".wind_save");
  strcat (photfile, ".phot");
  strcat (diagfile, ".diag");
  strcat (diskfile, ".disk.diag");


/* Define the sizes of the wind arrays */
  geo.coord_type = 1;
  geo.ndim = ndim = 30;
  geo.mdim = mdim = 30;
  NDIM = ndim;
  MDIM = mdim;
  NDIM2 = ndim * mdim;
  dfudge = 1e5;
  DFUDGE = dfudge;
  geo.adiabatic = 0;		// Default is set so that adiabatic cooling is not included in the wind
/* End of definition of wind arrays */

  /* Start logging of errors and comments */

  Log_init (diagfile);
  Log ("!!Python Version %s \n", VERSION);	//54f -- ksl -- Now read from version.h

/* Set initial values for everything in geo struct which basically defines the overall geometry */
  init_geo ();

/* Get the atomic data */

  strcpy (atomic_filename, "atomic/standard39");
  rdstr ("Atomic_data", atomic_filename);
  strcpy (geo.atomic_filename, atomic_filename);	//
  get_atomic_data (atomic_filename);

/* Initialize variables which are used in the main routine */
  wcycles = pcycles = 1;
  photons_per_cycle = 20000;

/* Initialize basis vectors for a cartesian coordinate system */
  x_axis[0] = 1.0;
  x_axis[1] = x_axis[2] = 0.0;
  y_axis[1] = 1.0;
  y_axis[0] = y_axis[2] = 0.0;
  z_axis[2] = 1.0;
  z_axis[1] = z_axis[0] = 0.0;

/* Gather input data */

  /* Describe the basic calculation in terms of the number of iterations which will
     be used to calculate the wind parameters and the number of iterations and wavelength
     range which will be used for the final spectrom.  Also describe the observer's views
     of the system */

  rdint
    ("Wind_type(0=SV,1=Sphere,2=Previous,3=Proga,4=Corona,5=knigge,6=thierry,7=yso)",
     &geo.wind_type);

  rdint ("photons_per_cycle", &photons_per_cycle);
  NPHOT = photons_per_cycle;	// For now set NPHOT to be be photons/cycle --> subcycles==1

  photons_per_cycle = (photons_per_cycle / NPHOT) * NPHOT;
  if (photons_per_cycle < NPHOT)
    photons_per_cycle = NPHOT;
  subcycles = photons_per_cycle / NPHOT;
  Log ("Photons_per_cycle adjusted to %d\n", photons_per_cycle);

  rdint ("Ionization_cycles", &wcycles);

  rdint ("spectrum_cycles", &pcycles);

  if (wcycles == 0 && pcycles == 0)
    exit (0);			//There is really nothing to do!

/* Allocate the memory for the photon structure now that NPHOT is established */

  p = (PhotPtr) calloc (sizeof (p_dummy), NPHOT);

  if (p == NULL)
    {
      Error
	("There is a problem in allocating memory for the photon structure\n");
      exit (0);
    }


  if (geo.wind_type == 2)
    {				//then we will start with a previously calculated model

      strcpy (old_windsavefile, "root_only");
      rdstr ("Old_windfile", old_windsavefile);
      strcat (old_windsavefile, ".wind_save");

/* Note that wind_read allocates space for w for previously calculated models. 
See the discussion in wind_read or py_wind for the tortured explanation
for why windread was changed from wind_read(w,old_windsavefile) to 
wind_read(old_windsavefile) 02apr ksl */

      wind_read (old_windsavefile);	//program will exit if unable to read the file
      geo.wind_type = 2;	// after wind_read one will have a different wind_type otherwise
      w = wmain;

    }
  else
    {
      /* Define the coordinate system for the grid and allocate memory for the wind structure */
      rdint
	("Coord.system(0=spherical,1=cylindrical,2=spherical_polar,3=cyl_var)",
	 &geo.coord_type);

      rdint ("Wind.dim.in.x_or_r.direction", &geo.ndim);
      if (geo.coord_type)
	{
	  rdint ("Wind.dim.in.z_or_theta.direction", &geo.mdim);
	  if (geo.mdim < 4)
	    {
	      Error
		("python: geo.mdim must be at least 4 to allow for boundaries\n");
	      exit (0);
	    }
	}
      else
	geo.mdim = 1;

    }


  {

    rdint
      ("Wind_ionization(0=on.the.spot,1=LTE,2=fixed,3=recalc,4=recalc[ground_multiplicities])",
       &geo.ioniz_mode);

    if (geo.ioniz_mode == 2)
      {
	rdstr ("Fixed.concentrations.filename", &geo.fixed_con_file[0]);
      }
    if (geo.ioniz_mode > 4)
      {
	Error ("python not up to date on ionization modes > 3\n");
	exit (0);
      }
/*Normally, geo.partition_mode is set to -1, which means that partition functions are calculated to take
full advantage of the data file.  This means that in calculating the partition functions, the information
on levels and their multiplicities is taken into account.   */


    geo.partition_mode = -1;

/* Setting geo.ioniz_mode to 4, results in a test
case that perverts this, so that only the ground state multiplicities are used.
is calculated */
    if (geo.ioniz_mode == 77)	// Use recalc but only ground state multiplicities
      {
	geo.ioniz_mode = 3;
	geo.partition_mode = 0;
      }

    rdint
      ("Line_transfer(0=pure.abs,1=pure.scat,2=sing.scat,3=escape.prob,6=macro_atoms,7=macro_atoms+aniso.scattering)",
       &geo.line_mode);

/* ?? Next section seems rather a kluge.  Why don't we specifty the underlying variables explicitly 
It also seems likely that we have mixed usage of some things, e.g ge.rt_mode and geo.macro_simple */

    /* For now handle scattering as part of a hidden line transfermode ?? */
    if (geo.line_mode == 4)
      {
	geo.scatter_mode = 1;	// Turn on anisotropic scattering
	geo.line_mode = 3;	//  Drop back to escape probabilities
	geo.rt_mode = 1;	// Not macro atom (SS)
      }
    else if (geo.line_mode == 5)
      {
	geo.scatter_mode = 2;	// Thermal trapping model
	geo.line_mode = 3;	// Single scattering model is best for this mode
	geo.rt_mode = 1;	// Not macro atom (SS) 
      }
    else if (geo.line_mode == 6)
      {
	geo.scatter_mode = 0;	// isotropic
	geo.line_mode = 3;	// Single scattering
	geo.rt_mode = 2;	// Identify macro atom treatment (SS)
	geo.macro_simple = 0;	// We don't want the all simple case (SS)
      }
    else if (geo.line_mode == 7)
      {
	geo.scatter_mode = 2;	// thermal trapping
	geo.line_mode = 3;	// Single scattering
	geo.rt_mode = 2;	// Identify macro atom treatment (SS)
	geo.macro_simple = 0;	// We don't want the all simple case (SS)
      }
    else if (geo.line_mode == 8)
      {
	geo.scatter_mode = 0;	// isotropic
	geo.line_mode = 3;	// Single scattering
	geo.rt_mode = 2;
	geo.macro_simple = 1;	//This is for test runs with all simple ions (SS)
      }
    else
      {
	geo.scatter_mode = 0;	//isotrpic
	geo.rt_mode = 1;	//Not macro atom (SS)
      }


/*57h -- Next line added for speed when no macro atoms to prevent bf calculation of
macro_estimaters.  Is this OK, Stuart??   */

    if (nlevels_macro == 0)
      geo.macro_simple = 1;	// Make everything simple if no macro atoms -- 57h

    //SS - initalise the choice of handling for macro pops.
    geo.macro_ioniz_mode = 0;

    // Determine what radiation sources there are
    rdint ("Star_radiation(y=1)", &geo.star_radiation);
    rdint ("Disk_radiation(y=1)", &geo.disk_radiation);
    rdint ("Boundary_layer_radiation(y=1)", &geo.bl_radiation);
    rdint ("Wind_radiation(y=1)", &geo.wind_radiation);

    if (!geo.star_radiation && !geo.disk_radiation && !geo.bl_radiation
	&& !geo.bl_radiation)
      {
	Error ("python: No radiation sources so nothing to do but quit!\n");
	exit (0);
      }
    /* 
       With the macro atom approach we won't want to generate photon 
       bundles in the wind so switch it off here. (SS)
     */

    if (geo.rt_mode == 2)
      {
	Log
	  ("python: Using Macro Atom method so switching off wind radiation.\n");
	geo.wind_radiation = 0;
      }


    if (geo.star_radiation)
      {
	rdint ("Rad_type_for_star(0=bb,1=models)_to_make_wind",
	       &geo.star_ion_spectype);
	if (geo.star_ion_spectype == 0)
	  geo.star_ion_spectype = -1;	// bb
	else
	  {
	    rdstr ("Model_file", model_list);
	    get_models (model_list, 2, &geo.star_ion_spectype);
	  }
      }
    else
      {
	geo.star_ion_spectype = -3;	// No radiation
      }


    if (geo.disk_radiation)
      {
	rdint ("Rad_type_for_disk(0=bb,1=models)_to_make_wind",
	       &geo.disk_ion_spectype);
	if (geo.disk_ion_spectype == 0)
	  geo.disk_ion_spectype = -1;	// bb
	else
	  {
	    rdstr ("Model_file", model_list);
	    get_models (model_list, 2, &geo.disk_ion_spectype);
	  }
      }
    else
      geo.disk_ion_spectype = -3;	// No radiation from disk


    if (geo.bl_radiation)
      {
	rdint ("Rad_type_for_bl__(0=bb,1=models)_to_make_wind",
	       &geo.bl_ion_spectype);
	if (geo.bl_ion_spectype == 0)
	  geo.bl_ion_spectype = -1;	// bb
	else
	  {
	    rdstr ("Model_file", model_list);
	    get_models (model_list, 2, &geo.bl_ion_spectype);
	  }
      }
    else
      {
	geo.bl_ion_spectype = -3;	// No radiation from bl
      }

    /* Describe the basic binary star system */

    geo.binary_system = 1;
    rdint ("Overall_system_type(0=single_object,1=binary)",
	   &geo.binary_system);
    rddoub ("mstar(msol)", &geo.mstar);
    rddoub ("rstar(cm)", &geo.rstar);
    geo.rstar_sq = geo.rstar * geo.rstar;
    if (geo.star_radiation)
      rddoub ("tstar", &geo.tstar);

    if (geo.binary_system)
      {
	rddoub ("msec(msol)", &geo.m_sec);
	rddoub ("period(hr)", &geo.period);
      }

/* Describe the disk */

    geo.disk_type = 1;
    rdint
      ("disk.type(0=no.disk,1=standard.flat.disk,2=vertically.extended.disk)",
       &geo.disk_type);
    if (geo.disk_type)		/* Then a disk exists and it needs to be described */
      {
	if (geo.disk_radiation)
	  {
	    rddoub ("disk.mdot(msol/yr)", &geo.disk_mdot);
	    disk_illum = 0;
	    rdint
	      ("Disk.illumination.treatment(0=no.rerad,1=high.albedo,2=thermalized.rerad,3=analytic)",
	       &disk_illum);
	  }
	else
	  {
	    geo.disk_mdot = 0;
	    disk_illum = 0;
	  }
	/* 04aug ksl Until everything is initialized we need to stick to a simple disk, 
	 * while teff is being set up..  This is because some of the
	 * models, e.g. knigge have wind structures that depend on teff.
	 */
	geo.disk_illum = 0;
	if (disk_illum == 3)
	  {
	    geo.disk_illum = 3;
	  }

	rddoub ("disk.radmax(cm)", &geo.diskrad);
	Log ("geo.diskrad  %e\n", geo.diskrad);

/* If diskrad <= geo.rstar set geo.disk_type = 0 to make any disk transparent anyway. */
	if (geo.diskrad < geo.rstar)
	  {
	    Log
	      ("Disk radius is less than star radius, so assuming no disk)\n");
	    geo.disk_type = 0;
	  }
      }

    geo.disk_z0 = geo.disk_z1 = 0.0;
    if (geo.disk_type == 2)
      {				/* Get the additional variables need to describe a vertically extended disk */
	rddoub ("disk.z0(fractional.height.at.diskrad)", &geo.disk_z0);
	rddoub ("disk.z1(powerlaw.index)", &geo.disk_z1);
      }
    if (geo.disk_type == 0)
      {				/* There is no disk so set variables accordingly */
	geo.disk_radiation = 0;
	geo.diskrad = 0;
      }


    /* Describe the boundary layer */
    if (geo.bl_radiation)
      {
	xbl = geo.lum_bl =
	  0.5 * G * geo.mstar * geo.disk_mdot * MSOL * MSOL / (geo.rstar *
							       YR);
	rddoub ("lum_bl(ergs/s)", &geo.lum_bl);
	Log ("OK, the bl lum will be about %.2e the disk lum\n",
	     geo.lum_bl / xbl);
	rddoub ("t_bl", &geo.t_bl);
      }
    else
      {
	geo.lum_bl = 0;
	geo.t_bl = 0;
      }

    /* Describe the wind */
    rddoub ("wind.radmax(cm)", &geo.rmax);
    rddoub ("wind.t.init", &geo.twind);

    /* Convert all inputs to cgs units */

    geo.mstar *= MSOL;
    geo.m_sec *= MSOL;
    geo.disk_mdot *= MSOL / YR;
    geo.period *= 3600;

    geo.diskrad_sq = geo.diskrad * geo.diskrad;

/* Now get parameters that are specific to a given wind model

   Note: When one adds a new model, the only things that should be read in and modified
are parameters in geo.  This is in order to preserve the ability to continue a calculation
with the same basic wind geometry, without reading in all of the input parameters.  
*/

    if (geo.wind_type == 1)
      {
	get_stellar_wind_params ();
      }
    else if (geo.wind_type == 0)
      {
	get_sv_wind_params ();
      }
    else if (geo.wind_type == 3)
      {
	get_proga_wind_params ();
      }
    else if (geo.wind_type == 4)
      {
	get_corona_params ();
      }
    else if (geo.wind_type == 5)
      {
	get_knigge_wind_params ();
      }
    else if (geo.wind_type == 6)
      {
	get_thierry_params ();
      }
    else if (geo.wind_type == 7)
      {
	get_yso_wind_params ();
      }
    else if (geo.wind_type != 2)
      {
	Error ("python: Unknown wind type %d\n", geo.wind_type);
	exit (0);
      }



/* Now define the wind cones generically.  Note that thetamin and
  thetamax are measured from disk plane, and so if thetamin = 90 one ends
  up with a cone that is constant in rho.  
  56d -- ksl -- updated windcone definitions  */

//OLD    windcone[0].r_zero = geo.wind_rho_min; // Original wincone definition
//OLD   windcone[0].drdz = tan (geo.wind_thetamin);     // Original definition

    if (geo.wind_thetamin > 0.0)
      {
	windcone[0].dzdr = 1. / tan (geo.wind_thetamin);	// new definition
	windcone[0].z = (-geo.wind_rho_min / tan (geo.wind_thetamin));	// new definition
      }
    else
      {
	windcone[0].dzdr = INFINITY;	// new definition
	windcone[0].z = -INFINITY;;	// new definition
      }


//OLD   windcone[1].r_zero = geo.wind_rho_max;  // old definition
//OLD    windcone[1].drdz = tan (geo.wind_thetamax);    // old definition
// 060908 -- 57h -- fixed a small error in next section in definition of outer windcone
    if (geo.wind_thetamax > 0.0)
      {
	windcone[1].dzdr = 1. / tan (geo.wind_thetamax);	// new definition
	windcone[1].z = (-geo.wind_rho_max / tan (geo.wind_thetamax));	// new definition
      }
    else
      {
	windcone[1].dzdr = INFINITY;	// old definition
	windcone[1].z = -INFINITY;;	// new definition
      }


    geo.rmax_sq = geo.rmax * geo.rmax;

    /* Calculate additonal parameters associated with the binary star system */

    if (geo.binary_system)
      binary_basics ();

    if (geo.wind_type != 2)
      define_wind ();		// Do not reinit if you want to use old windfile

  }
  w = wmain;

  /* Check that the parameters which have been supplied for the star, disk and boundary layer will
     allow generation of photons where that is appropriate */

  if (geo.tstar <= 0.0)
    geo.star_radiation = 0;
  if (geo.disk_mdot <= 0.0)
    geo.disk_radiation = 0;
  if (geo.t_bl <= 0.0 || geo.lum_bl <= 0.0)
    geo.bl_radiation = 0;

  /* Next block added by SS August 04 in case diskrad = 0 is used to mean no disk at any point. */

  if (geo.diskrad <= 0.0)
    {
      geo.disk_type = 0;
      geo.disk_radiation = 0;
    }

  if (geo.star_radiation)
    Log ("There is a star which radiates\n");
  else
    Log ("The star in the system does not radiate\n");

  if (!geo.disk_type)
    Log ("There is no disk in the system \n");
  else if (!geo.disk_radiation)
    Log ("The disk exists, but only absorbs photons\n");
  else
    Log ("There is a disk which radiates and absorbs\n");

  if (geo.bl_radiation)
    Log ("There is a boundary layer which radiates\n");
  else
    Log ("There is no boundary layer\n");



/* Describe the spectra which will be extracted and the way it will be extracted */

  nangles = 4;
  angle[0] = 10;
  angle[1] = 30.;
  angle[2] = 60.;
  angle[3] = 80.;
  for (n = 4; n < NSPEC - 3; n++)
    angle[n] = 45;
  for (n = 0; n < NSPEC - 3; n++)
    {
      phase[n] = 0.5;
      scat_select[n] = 1000;
      top_bot_select[n] = 0;
    }
  swavemin = 1450;
  swavemax = 1650;

  select_extract = 1;
  select_spectype = 1;

  if (pcycles > 0)
    {
      if (geo.star_radiation)
	{
	  rdint
	    ("Rad_type_for_star(0=bb,1=models,2=uniform)_in_final_spectrum",
	     &geo.star_spectype);
	  if (geo.star_spectype == 0)
	    geo.star_spectype = -1;	// bb
	  else if (geo.star_spectype == 2)
	    geo.star_spectype = -2;	// uniform
	  else
	    {
	      rdstr ("Model_file", model_list);
	      get_models (model_list, 2, &geo.star_spectype);
	    }
	}
      else
	{
	  geo.star_spectype = -3;	// No radiation from star
	}

      if (geo.disk_radiation)
	{
	  rdint
	    ("Rad_type_for_disk(0=bb,1=models,2=uniform)_in_final_spectrum",
	     &geo.disk_spectype);
	  if (geo.disk_spectype == 0)
	    geo.disk_spectype = -1;	// bb
	  else if (geo.disk_spectype == 2)
	    geo.disk_spectype = -2;	// uniform
	  else
	    {
	      rdstr ("Model_file", model_list);
	      get_models (model_list, 2, &geo.disk_spectype);
	    }
	}
      else
	{
	  geo.disk_spectype = -3;	// No radiation
	}

      if (geo.bl_radiation)
	{
	  rdint
	    ("Rad_type_for_bl__(0=bb,1=models,2=uniform)_in_final_spectrum",
	     &geo.bl_spectype);
	  if (geo.bl_spectype == 0)
	    geo.bl_spectype = -1;	// bb
	  else if (geo.bl_spectype == 2)
	    geo.bl_spectype = -2;	// uniform
	  else
	    {
	      rdstr ("Model_file", model_list);
	      get_models (model_list, 2, &geo.bl_spectype);
	    }
	}
      else
	{
	  geo.bl_spectype = -3;	// No radiation    
	}



      rddoub ("spectrum_wavemin", &swavemin);
      rddoub ("spectrum_wavemax", &swavemax);
      if (swavemin > swavemax)
	{
	  swavemax = swavemin;
	  swavemin = swavemax;
	}

      /* SS June 04: convert these to frequencies and store for use
         in computing macro atom and k-packet emissivities. */

      em_rnge.fmin = C / (swavemax * 1.e-8);
      em_rnge.fmax = C / (swavemin * 1.e-8);

      geo.matom_radiation = 0;	//initialise for ionization cycles - don't use pre-computed emissivities for macro-atom levels/ k-packets.



      rdint ("no_observers", &nangles);
      if (nangles < 1 || nangles > NSPEC - MSPEC)
	{
	  Error ("no_observers %d should not be > %d or <0\n", nangles,
		 NSPEC - MSPEC);
	  exit (0);
	}

      for (n = 0; n < nangles; n++)
	rddoub ("angle(0=pole)", &angle[n]);

      /* 05apr-ksl-56--For diagnositic reasons I have left questions regarding phase
       * even for systems which are not binaries.  Phase 0 in this case corresponds to
       * an extraction direction which is in the xz plane
       */

      for (n = 0; n < nangles; n++)
	rddoub ("phase(0=inferior_conjunction)", &phase[n]);

      rdint ("live.or.die(0).or.extract(anything_else)", &select_extract);
      if (select_extract != 0)
	{
	  select_extract = 1;
	  Log ("OK, extracting from specific angles\n");
	}
      else
	Log ("OK, using live or die option\n");

/* Select spectra with certain numbers of scatterings.  See extract 1997 aug 28 ksl */

      strcpy (yesno, "n");
      rdstr ("Select_specific_no_of_scatters_in_spectra(y/n)", yesno);
      if (yesno[0] == 'y')
	{
	  Log
	    ("OK n>MAXSCAT->all; 0<=n<MAXSCAT -> n scatters; n<0 -> >= |n| scatters\n");
	  for (n = 0; n < nangles; n++)
	    {
	      rdint ("Select_scatters", &scat_select[n]);
	    }
	}
      strcpy (yesno, "n");
      rdstr ("Select_photons_by_position(y/n)", yesno);
      if (yesno[0] == 'y')
	{
	  Log
	    ("OK 0->all; -1 -> below; 1 -> above the disk, 2 -> specific location in wind\n");
	  for (n = 0; n < nangles; n++)
	    {
	      rdint ("Select_location", &top_bot_select[n]);
	      if (top_bot_select[n] == 2)
		{
		  Log
		    ("Warning: Make sure that position will be in wind, or no joy will be obtained\n");
		  rddoub ("rho(cm)", &rho_select[n]);
		  rddoub ("z(cm)", &z_select[n]);
		  rddoub ("azimuth(deg)", &az_select[n]);
		  rddoub ("r(cm)", &r_select[n]);

		}
	    }
	}
    }

  /* Select the units of the output spectra.  This is always needed */

  rdint ("spec.type(flambda(1),fnu(2),basic(other)", &select_spectype);
  if (select_spectype == 1)
    {
      Log ("OK, generating flambda at 100pc\n");
    }
  else if (select_spectype == 2)
    {
      Log ("OK, generating fnu at 100 pc\n");
    }
  else
    Log ("OK, basic Monte Carlo spectrum\n");

/* Determine whether to produce additonal diagnositics */

  rdint ("Extra.diagnostics(0=no)", &diag_on_off);


/* 57h -- New section of inputs to provide more control over how the program is
run -- 07jul -- ksl
*/

  istandard = 1;
  SMAX_FRAC = 0.5;
  DENSITY_PHOT_MIN = 1.0;
  keep_photoabs = 1;
  rdint ("Use.standard.care.factors(1=yes)", &istandard);

  if (!istandard)
    {
      rddoub ("Fractional.distance.photon.may.travel", &SMAX_FRAC);
      rddoub ("Lowest.ion.density.contributing.to.photoabsorption",
	      &DENSITY_PHOT_MIN);
      rdint ("Keep.photoabs.during.final.spectrum(1=yes)", &keep_photoabs);
    }
/* Wrap up and save all the inputs */

  if (strncmp (root, "mod", 3) == 0)
    cpar ("mod.pf");
  else if (strncmp (root, "dummy", 5) == 0)
    {
      cpar ("dummy.pf");
      exit (0);
    }
  else
    cpar ("python.pf");

/* OK all inputs have been obtained at this point and the inuts have been copied to "mod.pf" or "python.pf" */

  if (diag_on_off)
    {
      /* Open a diagnostic file or files.  These are all fixed files */
      open_diagfile ();
    }

/* initialize the random number generator */
//      srand( (n=(unsigned int) clock()));  
  srand (1084515760);

/* Determine the frequency range which will be used to establish the ionization balance of the wind */

  freqmin = C / 400000e-8;	/*400000 A */

  tmax = TSTAR;
  if (geo.twind > tmax)
    tmax = geo.twind;
  if (geo.tstar > tmax)
    tmax = geo.tstar;
  if (geo.t_bl > tmax && geo.lum_bl > 0.0)
    tmax = geo.t_bl;
  if ((0.488 * tdisk (geo.mstar, geo.disk_mdot, geo.rstar)) > tmax)
    tmax = 0.488 * tdisk (geo.mstar, geo.disk_mdot, geo.rstar);
  freqmax = BOLTZMANN * tmax / H * 10.;
  if (freqmax < 2.0 * 54.418 / HEV)
    {
      Log ("Increasing maximum frequency to twice the Helium edge\n");
      freqmax = 2.0 * 54.418 / HEV;
    }
  else
    Log ("Maximum frequency %8.2e determined by T %8.2e\n", freqmax, tmax);



/* Next section is a kluge to initialize a second disk structure */
  disk_init (geo.rstar, geo.diskrad, geo.mstar, geo.disk_mdot, freqmin,
	     freqmax, 0, &geo.f_disk);
  for (n = 0; n < NRINGS; n++)
    {
      qdisk.r[n] = disk.r[n];
      qdisk.t[n] = disk.t[n];
      qdisk.g[n] = disk.g[n];
      qdisk.v[n] = disk.v[n];
      qdisk.heat[n] = disk.heat[n];
      qdisk.nphot[n] = disk.nphot[n];
      qdisk.nhit[n] = disk.nhit[n];
    }

/* 04aug -- ksl -- now that everything is initialized, we set geo.disk_illum
 */
  geo.disk_illum = disk_illum;


  timer ();
  Log ("Timer initiated\n");


/* THE CALCULATION OF THE IONIZATION OF THE WIND */

  geo.ioniz_or_extract = 1;	//SS July 04 - want to compute MC estimators during ionization cycles

  for (icycles = 0; icycles < wcycles; icycles++)
    {				/* This allows you to build up photons in bunches */


      Log ("!!Python: Begining cycle %d of %d for defining wind\n",
	   icycles, wcycles);

      /* Initialize all of the arrays, etc, that need initialization for each cycle
       */

      spectrum_init (freqmin, freqmax, nangles, angle, phase, scat_select,
		     top_bot_select, select_extract, rho_select, z_select,
		     az_select, r_select);

      /* Zero the arrays that store the heating of the disk */
      for (n = 0; n < NRINGS; n++)
	{
	  qdisk.heat[n] = qdisk.nphot[n] = 0;
	}

//06may -- Changed call as part of restructuring of all of the structures -- ksl
      wind_rad_init ();		/*Zero the parameters pertaining to the radiation field */

#if DEBUG
      ispy_init ("python", icycles);
#endif

      n_ioniz = 0.0;
      lum_ioniz = 0.0;
      ztot = 0.0;

      /* Now execute each subcycle */

      for (j = 0; j < subcycles; j++)
	{
	  Log
	    ("Subcycle %d of %d in Cycle %d of %d for defining wind structure\n",
	     j, subcycles, icycles, wcycles);

	  if (!geo.wind_radiation || (icycles == 0 && geo.wind_type != 2))
	    iwind = -1;		/* Do not generate photons from wind */
	  else if (j == 0)
	    iwind = 1;		/* Create wind photons and force a reinitialization of wind parms */
	  else
	    iwind = 0;		/* Create wind photons but do not force reinitialization */

	  /*photons_per_cycle is the number of photon bundles which will equal the luminosity; 
	   * 0 => for ionization calculation 
	   */

	  define_phot (p, w, freqmin, freqmax, photons_per_cycle, 0,
		       iwind, 1);


	  // Moved photon checks to separate routine ... but not sure they are needed at all
	  photon_checks (p, freqmin, freqmax);

	  zz = 0.0;
	  for (nn = 0; nn < NPHOT; nn++)
	    {
	      zz += p[nn].w;
	    }

	  Log
	    ("!!python: Total photon luminosity before transphot %8.2e\n",
	     zz);

	  ztot += zz;		// Total luminosity in all subcycles

	  /* kbf_need determines how many & which bf processes one needs to considere.  It was introduced
	   * as a way to speed up the program.  It has to be recalcuated evey time one changes
	   * freqmin and freqmax
	   */

	  kbf_need (freqmin, freqmax);

	  trans_phot (w, p, 0);

	  /*Determine how much energy was absorbed in the wind */
	  zze = zzz = 0.0;
	  for (nn = 0; nn < NPHOT; nn++)
	    {
	      zzz += p[nn].w;
	      if (p[nn].istat == P_ESCAPE)
		zze += p[nn].w;
	    }
	  Log
	    ("!!python: Total photon luminosity after transphot %8.2e (diff %8.2e). Radiated luminosity %8.2e \n",
	     zzz, zzz - zz, zze);

#if DEBUG
	  wind_rad_summary (w, windradfile, "a");
#endif
	  spectrum_create (p, freqmin, freqmax, nangles, select_extract);


	}

      /* End of the subcycle loop */


#if DEBUG
      ispy_close ();
#endif


/* 04mar-- ksl -- add section to dump heating of disk.  
 * 04aug-- ksl -- fixed several small problems with this
 * 		Factor of 2 in area calculation reflects
 * 		fact that disk has two sides
 * 04dec -- ksl -- created variable ztot so fractional heating
 * 		is correct if multiple subcyles
 */
      qptr = fopen (diskfile, "w");
      fprintf (qptr,
	       "# r       zdisk   t_disk   heat    nhit fraction_emitted  t_heat\n");
      for (n = 0; n < NRINGS; n++)
	{
	  theat =
	    qdisk.heat[n] / (2. * PI *
			     (qdisk.r[n + 1] * qdisk.r[n + 1] -
			      qdisk.r[n] * qdisk.r[n]));
	  theat = pow (theat / STEFAN_BOLTZMANN, 0.25);
	  fprintf (qptr, "%8.3e %8.3e %8.3e %8.3e %5d %8.3e %8.3e\n",
		   qdisk.r[n], zdisk (qdisk.r[n]), qdisk.t[n],
		   qdisk.heat[n], qdisk.nhit[n],
		   qdisk.heat[n] * NRINGS / ztot, theat);
	}
      fclose (qptr);

/* Completed writing file describing disk heating */

      Log
	("!!python: Number of ionizing photons %g lum of ionizing photons %g\n",
	 n_ioniz, lum_ioniz);

      wind_update (w);
      if (diag_on_off)
	{
	  strcpy (dummy, "");
	  sprintf (dummy, "python%02d.wind_save", icycles);
	  wind_save (w, dummy);
	  Log ("Saved wind structure in %s\n", dummy);
	}


    }

  if (wcycles > 0)
    {
      Log ("Finished creating spectra\n");
      spectrum_summary (wspecfile, "w", 0, 5, 0, 1.);
      phot_gen_sum (photfile, "w");	/* Save info about the way photons are created and absorbed
					   by the disk */
    }

  if (wcycles > 0 || geo.wind_type != 2)
    {
      wind_save (w, windsavefile);
      Log ("Saved wind structure in %s\n", windsavefile);
    }


  Log (" Completed wind creation.  The elapsed TIME was %f\n", timer ());

/* THE CALCULATION OF A DETAILED SPECTRUM IN A SPECIFIC REGION OF WAVELENGTH SPACE */

  freqmax = C / (swavemin * 1.e-8);
  freqmin = C / (swavemax * 1.e-8);


  /* Perform the initilizations required to handle macro-atoms during the detailed
     calculation of the spectrum.  

     Next lines turns off macro atom estimators and other portions of the code that are
     unnecessary during spectrum cycles.  */

  geo.ioniz_or_extract = 0;

/* 57h -- 07jul -- Next steps to speed up extraction stage */
  if (!keep_photoabs)
    {
      DENSITY_PHOT_MIN = -1.0;	// Do not calculated photoabsorption in detailed spectrum 
    }

  /*Switch on k-packet/macro atom emissivities  SS June 04 */

  if (geo.rt_mode == 2)
    {
      geo.matom_radiation = 1;
    }

  /* Finished initializations required for macro-atom approach */

  /* Calculate and store which bf processess need to be considered in each cell
   * Note that this is not macro-specific but is just to speed the program up.
   */

  kbf_need (freqmin, freqmax);


  spectrum_init (freqmin, freqmax, nangles, angle, phase, scat_select,
		 top_bot_select, select_extract, rho_select, z_select,
		 az_select, r_select);


  /* Begin cycle to create the detailed spectrum */

  for (icycles = 0; icycles < pcycles; icycles++)
    {				/* This allows you to build up photons in bunches */


#if DEBUG
      ispy_init ("python", icycles + 1000);
#endif

      Log ("!!Cycle %d of %d to calculate a detailed spectrum\n", icycles,
	   pcycles);

      if (!geo.wind_radiation)
	iwind = -1;		/* Do not generate photons from wind */
      else if (icycles == 0)
	iwind = 1;		/* Create wind photons and force a reinitialization of wind parms */
      else
	iwind = 0;		/* Create wind photons but do not force reinitialization */

      /*NPHOT*pcycles is the number of photon bundles which will equal the luminosity, 1 implies for spectrum calculation */
      define_phot (p, w, freqmin, freqmax, NPHOT * pcycles, 1, iwind, 0);

      for (icheck = 0; icheck < NPHOT; icheck++)
	{
	  if (sane_check (p[icheck].freq))
	    {
	      Error
		("python after define phot: unnatural frequency for photon %d\n",
		 icheck);
	    }
	}


      trans_phot (w, p, select_extract);

#if DEBUG
      wind_rad_summary (w, windradfile, "a");
#endif

      spectrum_create (p, freqmin, freqmax, nangles, select_extract);

/* Write out the detailed spectrum each cycle so that one can see the statistics build up! */
      renorm = ((double) (pcycles)) / (icycles + 1.0);
      spectrum_summary (specfile, "w", 0, nspectra - 1, select_spectype,
			renorm);
    }


  /* End cycle to calculate detailed spectrum */

  phot_gen_sum (photfile, "a");

/* 57h - 07jul -- ksl -- Write out the freebound information */

  fb_save ("recomb.save");

/* Finally done */

  Log ("Completed entire program.  The elapsed TIME was %f\n", timer ());
  return EXIT_SUCCESS;
}


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	init_geo initializes the geo structure to something that is semi-reasonable
Arguments:		

Returns:
 
Description:	
	Initial values for all of the variables that are not part of the individual
	wind descriptions that are actully read(!) into the program should be 
	created here.  The derived values are not needed.

		
Notes:


History:
 	98dec	ksl	Coded and debugged.  Much of code was copied from old main routine for
			python
	04dec	ksl	This is probably still not completely up to date, but have
			added some initializations 
**************************************************************/

int
init_geo ()
{
  geo.wind_type = 0;		// Schlossman and Vitello

  geo.star_ion_spectype = geo.star_spectype
    = geo.disk_ion_spectype = geo.disk_spectype
    = geo.bl_ion_spectype = geo.bl_spectype = 0;

  geo.log_linear = 0;		/* Set intervals to be logarithmic */

  geo.rmax = 1e11;
  geo.rmax_sq = geo.rmax * geo.rmax;
  geo.rstar = 7e8;
  geo.rstar_sq = geo.rstar * geo.rstar;
  geo.mstar = 0.8;
  geo.m_sec = 0.4;
  geo.period = 3.2;
  geo.tstar = 40000;
  geo.twind = 40000;
  geo.wind_mdot = 1.e-9;

  geo.ioniz_mode = 3;		/* default is on the spot and find the best t */
  geo.line_mode = 3;		/* default is escape probabilites */

  geo.star_radiation = 0;	/* 1 implies star will radiate */
  geo.disk_radiation = 0;	/* 1 implies disk will radiate */
  geo.bl_radiation = 0;		/*1 implies boundary layer will radiate */
  geo.wind_radiation = 0;	/* 1 implies wind will radiate */

  geo.disk_type = 0;		/*1 implies existence of a disk for purposes of absorbtion */
  geo.diskrad = 2.4e10;
  geo.disk_mdot = 1.e-8;

  geo.t_bl = 100000.;
  geo.mstar = 0.8;


  strcpy (geo.atomic_filename, "none");	// 54e
  strcpy (geo.fixed_con_file, "none");	// 54e
  return (0);
}

/*
Perform some simple checks on the photon distribution just produced.

History:
	01	ksl	Removed from main routine
	02jul	ksl	Loosened frequency limits to reflect the
			fact that in some cases, e.g. those in
			which the photon distribution has been split
			into small energy segments, Doppler shifts
			move photons out of that region.


*/
int
photon_checks (p, freqmin, freqmax)
     PhotPtr p;
     double freqmin, freqmax;
{
  int nnn, nn;
  double lum_ioniz;
  int n_ioniz;


  nnn = 0;


  /* Next two lines are to allow for fact that photons generated in
   * a frequency range may be Doppler shifted out of that range, especially
   * if they are disk photons generated right up against one of the frequency
   * limits
   * 04aug--ksl-increased limit from 0.02 to 0.03, e.g from 6000 km/s to 9000 km/s
   * shift.
   */

  freqmax *= (1.03);
  freqmin *= (0.97);

  for (nn = 0; nn < NPHOT; nn++)
    {
      if (H * p[nn].freq > ion[0].ip)
	{
	  lum_ioniz += p[nn].w;
	  n_ioniz += p[nn].w / (H * p[nn].freq);
	}
      if (sane_check (p[nn].freq) != 0)
	{
	  Error ("photon_check: freq  %d %e %e %e sane_check failure\n", nn,
		 p[nn].freq, freqmin, freqmax);
	  p[nn].freq = freqmax;
	  nnn++;
	  if (nnn > 100)
	    {
	      Error
		("photon_check: Exiting because too many bad photons generated\n");
	      exit (0);
	    }
	}
      if (p[nn].freq < freqmin || freqmax < p[nn].freq)
	{
	  Error ("photon_check: freq  %d %e %e %e\n", nn, p[nn].freq,
		 freqmin, freqmax);
	  p[nn].freq = freqmax;
	  nnn++;
	  if (nnn > 100)
	    {
	      Error
		("photon_check: Exiting because too many bad photons generated\n");
	      exit (0);
	    }
	}
    }

  return (0);
}
