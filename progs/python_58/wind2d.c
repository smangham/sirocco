#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	define_wind  initializes the structure which characterizes the wind.  

Arguments:		
	WindPtr w;			The structure which defines the wind in Python
 
Returns:
 
Description:

	This the basic routine which initiallizes the wind structure, calculating the
	salient properties of the wind in each grid cell. 
		
Notes:  This routine has been modified a lot as multiple coordinate systems were implemented
        Now it is mostly a driver for routines found elsewhere.
	

History:
 	97jan      ksl	Coding on python began.
 	97sep13	ksl	Added the possibility of fixing the concentrations of the ions to specific
			values
 	97nov19 ksl	Fixed dmdot_da to account for the fact that mass is lost from both sides
                        of the disk.
 	98mar20	ksl	Fixed calculation of wind_midx and wind_midz.  New code is actually
 			in wind_complete
 	98may24	ksl	Modified calculations of wind volumes to force the edges to have zero volume
 			because the edge points are actually only used for interpolation and photons
 			do not interact in them.
 	98dec	ksl	Eliminated wind element .rzero, Added capability to define linear or logarithmic
 			intervals on the fly  Isolated all dependences on Shlossman and Vitello 
 			prescription of the wind to subroutine calls, which are in the file sv.c
 	98dec22	ksl	Added routines to allow definition of a one d stellar wind.  Weights inside
 			star artificially set to 0.5
	99dec29 ksl	Added routines for proga's wind 
	00oct06	ksl	Added routines for a corona
	01mar13	ksl	Added routines for knigge's wind
	01dec26	ksl	Deleted wind_summary routine since I haven't looked at this output for several 
			years (largely replaced by py_wind)
	02feb	ksl	Made changes to allow unequal dimensions in two directions
        04May   SS      Factor of two increase in the cell volumes to allow for space above and below plane.
	04aug	ksl	Began removing portions of this routine which were coordinate sytem independent
			to individial routines in wind.c.  Also redefined what is meant by .inwind.
        04Aug   SS      Modified the calulation of whether a cell is in the wind because it was 
                        missing cells when a finite disk is included (the cell corners could lie in the disk
                        but there could still be genuine wind material within the cell).
	05jul	ksl	56d -- Added checks of volume calculations for two d models as part of general
			effort to properly incoporate a vertically extneded disk.
	06may	ksl	Modified order in which items were calculated to prepare for splitting
			of the Wind structure.  Moved the creation of thw arrays here as well
			Allowed for wheich variables need to go into plasma

**************************************************************/


int
define_wind ()
{

  int i, j, n;
  double nh, rrstar;
  double x[3];
  double mdotbase, mdotwind, rr;
  int ierr;
  int n_vol, n_inwind, n_part;

  int nwind;
  int nplasma;

  WindPtr w;
  PlasmaPtr xplasma;

  /* In order to interpolate the velocity (and other) vectors out to geo.rmax, we need
     to define the wind at least one grid cell outside the region in which we want photons
     to propagate.  This is the reason we divide by NDIM-2 here, rather than NDIM-1 */

  NDIM = ndim = geo.ndim;
  MDIM = mdim = geo.mdim;
  NDIM2 = NDIM * MDIM;

  calloc_wind (NDIM2);
  w = wmain;



  if (geo.coord_type == SPHERICAL)
    {
      spherical_make_grid (w);
    }
  else if (geo.coord_type == CYLIND)
    {
      cylind_make_grid (w);
    }
  else if (geo.coord_type == RTHETA)
    {
      rtheta_make_grid (w);
    }
  else if (geo.coord_type == CYLVAR)
    {
      cylvar_make_grid (w);
    }
  else
    {
      Error ("define_wind: Don't know how to make coordinate type %d\n",
	     geo.coord_type);
    }

  for (n = 0; n < NDIM2; n++)
    {
      /* 04aug -- ksl -52 -- The next couple of lines are part of the changes
       * made in the program to allow more that one coordinate system in python 
       */
      model_velocity (w[n].x, w[n].v);
      model_vgrad (w[n].x, w[n].v_grad);
    }


  wind_complete (w);

/* wind_complete has to precede the volume calculations because the 1d vectors are used
in the volume calculations.  wind_complete itself is just a driver routine.  It has to
be called out as a separate routine, because the 1d vectors are not saved and have to be
recreated when a windfile is read into the program
06may ksl */


  /* Now define the valid volumes of each cell and also determine whether the cells are in all
   * or partially in the wind.
   *
   * Note - 05apr --55d -- ksl.  Previously there was a separate calculation here of whether a cell
   * was in the wind and its volume.  Given the increasingly complicated geometry when thick
   * disks were first allowed, Stuart had gone back to determining whether a cell was all or
   * partly in the wind by checking a bunch of positions in the cell.  But this is almost idetnaical 
   * to calculating the volume of the cell, and therefore I have moved this functionality
   * in the the respective volumes calculation.  At least then we will do the calculation the
   * same way both times.  . 
   */


  if (geo.coord_type == SPHERICAL)
    {
      spherical_volumes (w);
    }
  else if (geo.coord_type == CYLIND)
    {
      cylind_volumes (w);
    }
  else if (geo.coord_type == RTHETA)
    {
      rtheta_volumes (w);
    }
  else if (geo.coord_type == CYLVAR)
    {
      cylvar_volumes (w);
    }
  else
    {
      Error
	("wind2d.c: Don't know how to make volumes for coordinate type %d\n",
	 geo.coord_type);
    }


  n_vol = n_inwind = n_part = 0;
  for (n = 0; n < NDIM2; n++)
    {
      if (w[n].vol > 0.0)
	n_vol++;
      if (w[n].inwind == 0)
	n_inwind++;
      if (w[n].inwind == 1)
	n_part++;
    }

  Log
    ("wind2d: %3d cells of which %d are in inwind, %d partially in_wind, & %d with pos. vol\n",
     NDIM2, n_inwind, n_part, n_vol);

/* 56d --Now check the volume calculations for 2d wind models */
  if (geo.coord_type != SPHERICAL)
    {
      for (n = 0; n < NDIM2; n++)
	{
	  wind_n_to_ij (n, &i, &j);
	  if (i < (NDIM - 1) && j < (MDIM - 1))
	    {
	      n_inwind = 0;
	      if (where_in_wind (w[n].x) == 0)
		n_inwind++;
	      if (where_in_wind (w[n + 1].x) == 0)
		n_inwind++;
	      if (where_in_wind (w[n + MDIM].x) == 0)
		n_inwind++;
	      if (where_in_wind (w[n + MDIM + 1].x) == 0)
		n_inwind++;
	      if (w[n].vol == 0 && n_inwind > 0)
		{
		  Error
		    ("wind2d: Cell %d has %d corners in wind, but zero volume\n",
		     n, n_inwind);
		}
	      if (w[n].inwind == 1 && n_inwind == 4)
		{
		  Error
		    ("wind2d: Cell %d has 4 corners in wind, but is only partially in wind\n",
		     n);
		}
	    }
	}

    }


/* Now create the second structure, the one that is sized only to contain cells in the wind */

  if (CHOICE)
    {
      NPLASMA = n_vol;


    }
  else				/* Force NPLASMA to equal NDIM2 (for diagnostic reasons) */
    {
      NPLASMA = NDIM2;

    }
  calloc_plasma (NPLASMA);
  xplasma = plasmamain;
  create_maps (CHOICE);

  calloc_macro (NPLASMA);

/* 06may -- At this point we have calculated the volumes of all of the cells and it should
be optional which variables beyond here are moved to structures othere than Wind */


/* Now calculate parameters that need to be calculated at the center of the grid cell */


  for (n = 0; n < NPLASMA; n++)
    {

      nwind = plasmamain[n].nwind;
      stuff_v (w[nwind].xcen, x);

      plasmamain[n].rho = model_rho (x);


      nh = plasmamain[n].rho * rho2nh;
      plasmamain[n].t_r = geo.twind;
/* Initialize variables having to do with converence in initial stages */
      plasmamain[n].gain = 0.5;
      plasmamain[n].dt_e_old = 0.0;
      plasmamain[n].dt_e = 0.0;
      plasmamain[n].t_e = plasmamain[n].t_e_old = 0.9 * plasmamain[n].t_r;	//Lucy guess
      rrstar =
	1. - (geo.rstar * geo.rstar) / (x[0] * x[0] + x[1] * x[1] +
					x[2] * x[2]);
      if (rrstar > 0)
	{
	  plasmamain[n].w = 0.5 * (1 - sqrt (rrstar));
	}
      else
	plasmamain[n].w = 0.5;	//Modification to allow for possibility that grid point is inside star
      /* Determine the initial ionizations, either LTE or  fixed_concentrations */
      if (geo.ioniz_mode != 2)
	{			/* Carry out an LTE determination of the ionization */
	  ierr = ion_abundances (&plasmamain[n], 1);
	}
      else
	{			/* Set the concentrations to specified values */
	  ierr = ion_abundances (&plasmamain[n], 2);
	}
      if (ierr != 0)
	{
	  Error
	    ("wind_define after ion_abundances: cell %d rho %8.2e t_r %8.2 t_e %8.2e w %8.2e\n",
	     n, plasmamain[n].rho, plasmamain[n].t_r, plasmamain[n].t_e,
	     plasmamain[n].w);
	}
    }


/* Calculate the the divergence of the wind at the center of each grid cell */
  wind_div_v (w);

/* Now calculate the adiabatic cooling.  Note: adiabatic cooling is not used in
 * the program at present.  There are issues concerning how to incorporate it
 * into the macro atom approach, as well as questions concerning the appropriate
 * calculation.  If changes are made to this, then they must also be made in
 * the corresponding portion of wind_updates.  04nov -- ksl
 */

/*06may -- ksl -- This is awkward because liminsotities are now part of palam structure) */
  for (i = 0; i < NPLASMA; i++)
    {
      if (geo.adiabatic)
	{
	  nwind = plasmamain[i].nwind;
	  plasmamain[i].lum_adiabatic =
	    adiabatic_cooling (&w[nwind], plasmamain[i].t_e);
	}
      else
	plasmamain[i].lum_adiabatic = 0.0;
    }
  /* Calculate one over dvds */
  dvds_ave ();
  wind_check (w, -1);		// Check the wind for reasonability

/*
     Check that m_dot is correct.  This calculation is very approximate.  It only calculates mdot
     flowing up through the grid, i.e. any material that leaks out the side will not be counted.  A
     much better approach is needed, e.g. to calculate the mass flow in a spherical shell at the edge
     of the disk.  In addition, it uses a very inexact way to determine which grid points are in the
     wind.  Today, I simply modified the routine so the flux was calculated mid-way in grid space
     through the wind. ??? ksl 01mar15  

     04aug -- I have not written the routine to do the check that with cooridnate systems
     other than spherical, we are getting the correct mdot.  It should be possible to do this
     in a coordinate system independent manner simply by inteegrating using interpolated values
     and so this would be a good thing to fix.  But it is not used elsewhere so I have skipped
     it for now.

     05apr -- Looked at this again.  It is clearly a straightforwrd thing to do use the rho function
     which is below. One would simply integrate over various surface integrals to make it happen.  
*/

  if (geo.coord_type == CYLIND)
    {
      mdotwind = 0;
      mdotbase = 0;
      for (i = 0; i < NDIM - 1; i++)
	{
	  n = i * MDIM + (MDIM / 2);
// rr is router * router - rinner *rinner for this shell
	  rr =
	    w[n + MDIM].x[0] * w[n + MDIM].x[0] + w[n +
						    MDIM].x[1] * w[n +
								   MDIM].
	    x[1] - (w[n].x[0] * w[n].x[0] + w[n].x[1] * w[n].x[1]);
	  if (w[i * MDIM].inwind == 0)
	    {
	      nplasma = w[i * MDIM].nplasma;
	      mdotbase +=
		plasmamain[nplasma].rho * PI * rr * w[i * MDIM].v[2];
	    }
	  if (w[n].inwind == 0)
	    {
	      nplasma = w[n].nplasma;
	      mdotwind += plasmamain[nplasma].rho * PI * rr * w[n].v[2];
	    }
	}

      //Factor of 2 to allow for wind on both sides of disk
      Log ("m dot wind: Desired %g   Base %g Calculated %g\n",
	   geo.wind_mdot, 2. * mdotbase, 2. * mdotwind);
      mdot_wind (w, 1.e6, geo.wind_rmax / 2.);
    }
  else
    {
      Error
	("wind2d.c: Don't know how to check wind mdot for this coordtype %d\n",
	 geo.coord_type);
    }

  return (0);
}


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
 	where_in_grid locates the 1-d grid position of the photon. 

 Arguments:		
	double x[];
 Returns:
 	where_in_grid normally  returns the cell number associated with
 		a postion.  If the photon is in the grid this will be a positive
 		integer < NDIM*MDIM.
 	photon is inside the grid        -1
	photon is outside the grid       -2
 Description:	
	
		
 Notes:
	Where_in grid does not tell you whether the photon is in the wind or not. 

	What one means by inside or outside the grid may well be different
	for different coordinate systems.

	ksl--Not entirely clear why I don't make wig_x a vector.  Programming
	is not elegant at all.

 History:
 	97jan	ksl	Coding on python began.
	98jul	ksl	Removed WindPtr as one of arguments since not used
	98dec	ksl	Modified so that the argument is simply a position
	99jan	ksl	Modified to check that we are not trying to short
			circuit if we are asking for where_in_grid of same
			position as previously
	04aug	ksl	52a -- Now almos a shell, as moved to allow
			multiple coordinate systems
	05apr	ksl	55d -- Added spherical as a possiblity
 
**************************************************************/

int wig_n;
double wig_x, wig_y, wig_z;
int
where_in_grid (x)
     double x[];
{
  int n;
  double fx, fz;

  if (wig_x != x[0] || wig_y != x[1] || wig_z != x[2])	// Calculate if new position
    {

      if (geo.coord_type == CYLIND)
	{
	  n = cylind_where_in_grid (x);
	}
      else if (geo.coord_type == RTHETA)
	{
	  n = rtheta_where_in_grid (x);
	}
      else if (geo.coord_type == SPHERICAL)
	{
	  n = spherical_where_in_grid (x);
	}
      else if (geo.coord_type == CYLVAR)
	{
	  n = cylvar_where_in_grid (x, 0, &fx, &fz);
	}
      else
	{
	  Error ("where_in_grid: Unknown coord_type %d\n", geo.coord_type);
	  exit (0);
	}

/* Store old positions to short-circuit calculation if asked for same position more
than once */

      wig_x = x[0];
      wig_y = x[1];
      wig_z = x[2];
      wig_n = n;
    }

  return (wig_n);
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	vwind_xyz(w,p,v) finds the velocity vector v for the wind in cartesian 
	coordinates at the position of the photon p.
 
 Arguments:		
	WindPtr w;
	PhotPtr p;
	double v[];

Returns:
	0 	on successful completion
		the value of where in grid if the photon outside the wind grid.
 
Description:
	This routine carries out a bilinear interpolation of the wind velocity.  The velocities
	at the edges of the grid cells must have been defined elsewhere (e.g in wind_define).
	
	If the position is outside the wind grid, then the velocities will be returned as 0.	
	
Notes:
	For all of the 2d coordinate systems, the positions in the WindPtr array are
	caclulated in the xz plane.  As a result, the velocities in that array are
	calcuated there as well.  Hence, to obtain the velocity in cartesion coordinates
	one needs to do a simple rotation of the velocity vector.  

	For spherical (1d) coordinates the situation is complex, the philosopy that
	has been adopted is to let the program run even if the wind model is intrinsically
	2d.   The xyz positions of the grid are defined are defined so that
	are in the xz plane at an angle of 45 degrees to the x axis.  This was done
	so that if one used, say an sv wind, in spherical coordinates, one would select
	the wind velocities at a plausible location.   But it implies that the velocities
	that have been accepted are not purely radial as one might expect.  The choice
	that has been made is to assume that v[1] is preserved as a non zero value in
	making the projection.  An alternative might have been to preserve the angular
	rotation rate.
	
	Regardless, the fundamental challenge is to make sure the program works correctly 
	for the spherically symmetric case, e.g a simple stellar wind, which it did not
	before 56b.  

	If we ever implement a real 3d coordinate system, one will need to look at
	this routine again.

History:
 	97jan	ksl	Coding on python began.
	02feb	ksl	Modified to use new general purpose fraction routine.  
	04aug	ksl	52a -- Made a shell for choosing what coordinate
			system to calculate this in
	05jun	ksl	56b -- Despite the note above nothing had been done
			to fix vwind_xyz to handle spherical coordinates
			properly.  It was OK for any two-d system.  This 
			is now fixed, but see the notes above.
	06may	ksl	57+ -- Changed call to elimatnate passing the
			Wind array.  Use wmain instead.
 
**************************************************************/
int ierr_vwind = 0;

int
vwind_xyz (p, v)
     PhotPtr p;
     double v[];
{
  int i;
  double rho, r;
  double vv[3];
  double ctheta, stheta;
  double x, frac[4];
  int nn, nnn[4], nelem;


  // gives the correct result in all coord systems.


  /* 56d -- 05jul -- I had considerable problem with the
   * next routine for cylvar coords.  Coord_fraction must
   * produce a plausible result in all cases
   */

  coord_fraction (0, p->x, nnn, frac, &nelem);



  for (i = 0; i < 3; i++)
    {

      x = 0;
      for (nn = 0; nn < nelem; nn++)
	x += wmain[nnn[nn]].v[i] * frac[nn];

      vv[i] = x;
    }

  rho = sqrt (p->x[0] * p->x[0] + p->x[1] * p->x[1]);

  if (geo.coord_type == SPHERICAL)
    {				// put the velocity into cylindrical coordinates on the xz axis
      x = sqrt (vv[0] * vv[0] + vv[2] * vv[2]);	//v_r
      r = length (p->x);
      vv[0] = rho / r * x;
      vv[2] = p->x[2] / r * x;
    }
  else if (p->x[2] < 0)		// For 2d coord syatems, velocity is reversed if the photon is in the lower hemisphere.
    vv[2] *= -1;


  if (rho == 0)
    {				// Then we will not be able to project from a cylindrical ot a cartesian system
      Error
	("vwind_xyz: Cannot determine an xyz velocity on z axis. Returnin 0,0,v[2]\n");
      v[0] = v[1] = 0;
      v[2] = vv[2];
      return (0);
    }

  /* Now we have v in cylindrical coordinates, but we would like it in cartesian coordinates.  
     Note that could use project_from_cyl_xyz(p->x,vv,v) ?? */

  ctheta = p->x[0] / rho;
  stheta = p->x[1] / rho;
  v[0] = vv[0] * ctheta - vv[1] * stheta;
  v[1] = vv[0] * stheta + vv[1] * ctheta;
  v[2] = vv[2];

  if (sane_check (v[0]) || sane_check (v[1]) || sane_check (v[2]))
    {
      Error ("vwind_xyz: %f %f %f\n", v[0], v[1], v[2]);
    }


  return (0);
}


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	wind_div_v calculates the divergence of the velocity at the center of all the grid cells.  
Arguments:		

Returns:
 
Description:	
	This is one of the initialization routines for the wind.   
	Because the routine calls vwind_xyz, one must have previously 
	populated w[].v.  Also as a result of this fact, the routine 
	is not tightly tied to the SV prescription.
Notes:
	It is used to calculated PdV cooling of a grid cell 

	The divergence, like othe scalar quantities, does not
	need to be "rotated" differently in different coordinate
	systems

History:
	98mar20	ksl	First coded and debugged 
 	99dec1	ksl	Explicitly typed as integer and made the return 0
			Previously it was untyped, but returned div/EPSILON
			which made no sense to me.
	04mar	ksl	This routine traditionally has generated
			lots of errors.  One reason for could have beeen
			that the step size was too large.  So I made this
			smaller.  But the other reason is that the div v
			is often calculated outside the wind cone, as in
			kwd winds.  In those cases, I put an extrapolation
			of the wind velocity inside the wind cone, but
			this is not necessarily physical.  If the div v
			is interpolated this could be wrong.  I have 
			made some minor modifications though to suppress
			multiple errors if the div is being calculated
			outside the wind.
	05apr	ksl	55d -- Switched to using center positions for
			this calculation, rather than wind_midx, and 
			wind_midz.  Basically this was because I wanted
			to make things as independent of these two 
			arrays as possible in preparation for allowing
			more arbitrarily spaced grids.  Furthermore
			the previously formulatin was incorrect for
			rtheta components since wind_midz was actually
			theta.  It was also helpful for spherical models, 
			since it avoids the necessity of going from 1d to
			2d specification.

 
**************************************************************/

int wind_div_err = (-3);
int
wind_div_v (w)
     WindPtr w;
{
  int icell;
  double x_zero[3], v_zero[3], v[3];
  struct photon ppp;
  double div, delta;
  double length ();
  double xxx[3];
  for (icell = 0; icell < NDIM2; icell++)
    {
      /* Find the center of the cell */

      stuff_v (w->xcen, x_zero);

      delta = 0.01 * x_zero[2];	//new 04mar ksl

      /* Calculate v at midpoint */
      stuff_v (x_zero, ppp.x);
      vwind_xyz (&ppp, v_zero);
      /* Calculate dv_x/dx at this position */
      stuff_v (x_zero, ppp.x);
      ppp.x[0] += delta;
      vwind_xyz (&ppp, v);
      div = xxx[0] = (v[0] - v_zero[0]) / delta;
      /* Calculate dv_y/dy */
      stuff_v (x_zero, ppp.x);
      ppp.x[1] += delta;
      vwind_xyz (&ppp, v);
      div += xxx[1] = (v[1] - v_zero[1]) / delta;
      /* Calculate dv_z/dz */
      stuff_v (x_zero, ppp.x);
      ppp.x[2] += delta;
      vwind_xyz (&ppp, v);
      div += xxx[2] = (v[2] - v_zero[2]) / delta;
      w[icell].div_v = div;
      if (div < 0 && (wind_div_err < 0 || w->inwind == 0))
	{
	  Error
	    ("wind_div_v: div v %e is negative in cell %d. Major problem if inwind (%d) == 0\n",
	     div, icell, w->inwind);
	  wind_div_err++;
	}
    }

  return (0);
}


/* 
find the density of the wind at x

History:
	02feb	ksl Changed the calculation to use the general purpose routine fraction.
			This made a slight difference, an improvement I believe in how
			the program performs at the boundaries.  In the new version
			if you ask for rho at a position that is below say wind_midz
			one obtains the density at the center of the cell. In the old
			version it tried to extrapolate.  As far as I can determine
			the difference was relatively modest.
	04aug	ksl	52a -- modified to carry out a coordinate system independend
			determination of rho.  
	06may	ksl	57+ -- Modified for plasma structure but this is probably not 
			what one wants in the end.  There is no reason for the call
			to include w at all, since it is not used..

*/

double
rho (w, x)
     WindPtr w;
     double x[];
{
  int n, where_in_grid ();
  double dd;
  double frac[4];
  int nn, nnn[4], nelem;
  int nplasma;


  if (where_in_wind (x) != 0)	//note that where_in_wind is independent of grid.
    return (0.0);

  n = coord_fraction (1, x, nnn, frac, &nelem);

  if (n < 0)
    {
      dd = 0;
    }
  else
    {

      dd = 0;
      for (nn = 0; nn < nelem; nn++)
	nplasma = w[nnn[nn]].nplasma;
      dd += plasmamain[nplasma].rho * frac[nn];

    }


  return (dd);
}

/* Next is a routine to check the wind mass loss rate.  The mdot is checked in a plane running
from 0 to r, and in a sphere of radius r 

04aug -- ksl -- this looks coordinate system independnet but need to check
*/

#define NSTEPS 100
int
mdot_wind (w, z, rmax)
     WindPtr w;
     double z;			// The height (usually small) above the disk at which mdot will be calculated
     double rmax;		// The radius at which mdot will be calculated
{
  struct photon p;		//needed because vwind_xyz has a call which involves a photon
  double r, dr, rmin;
  double theta, dtheta;
  double den, rho ();
  double mdot, mplane, msphere;
  double x[3], v[3], q[3], dot ();
// Calculate the mass loss rate immediately above the disk
  rmin = geo.rstar;
  dr = (rmax - rmin) / NSTEPS;
  mdot = 0;
  p.x[1] = x[1] = 0.0;
  p.x[2] = x[2] = z;
  for (r = dr / 2; r < rmax; r += dr)
    {
      p.x[0] = x[0] = r;
      den = rho (w, x);
      vwind_xyz (&p, v);
      mdot += 2 * PI * r * dr * den * v[2];
//mdot+=2*PI*r*dr;
    }
  mplane = 2 * mdot;		// Factor of two because wind is in both hemispheres
// Calculate the mass loss rate in a sphere
  dtheta = PI / (2 * NSTEPS);
  mdot = 0;
  p.x[1] = x[1] = 0.0;
  q[1] = 0;
  for (theta = dtheta / 2; theta < (PI / 2.); theta += dtheta)
    {
      p.x[0] = x[0] = r * sin (theta);
      p.x[2] = x[2] = r * cos (theta);
      q[0] = sin (theta);
      q[2] = cos (theta);
      den = rho (w, x);
      vwind_xyz (&p, v);
      mdot += 2 * PI * rmax * rmax * sin (theta) * dtheta * den * dot (v, q);
    }
  msphere = 2 * mdot;		// Factor of two because wind is in both hemispheres
  Log ("Wind Mdot Desired %e Plane %e Sphere(%e) %e\n",
       geo.wind_mdot, mplane, rmax, msphere);
  return (0);
}


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	get_random_location is simply will produce a at a random place in
	a cell n.  
Arguments:		

Returns:
 
Description:	
	This is really just a driver routine for coordinate system specific
	routines
Notes:

History:
	04aug	ksl	52a -- created as part of project to allow multiple
			coordinate systems in python
	05apr	ksl	55d -- Added spherical option

 
**************************************************************/

int
get_random_location (n, x)
     int n;			// Cell in which to create postion
     double x[];		// Returned position
{
  int cylind_get_random (), rtheta_get_random ();

  if (geo.coord_type == CYLIND)
    {
      cylind_get_random_location (n, x);
    }
  else if (geo.coord_type == RTHETA)
    {
      rtheta_get_random_location (n, x);
    }
  else if (geo.coord_type == SPHERICAL)
    {
      spherical_get_random_location (n, x);
    }
  else if (geo.coord_type == CYLVAR)
    {
      cylvar_get_random_location (n, x);
    }
  else
    {
      Error ("get_random_location: Don't know this coord_type %d\n",
	     geo.coord_type);
      exit (0);
    }

  return (0);
}
