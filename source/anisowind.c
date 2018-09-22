
/***********************************************************/
/** @file  anisowind.c
 * @author ksl
 * @date   March, 2018
 *
 * @brief  Routines to implement anisotropic scattering in
 * the wind
 *
 *
 * Python supports several ways to determine a new photon
 * direction when a photon scatters in thie wind.  These
 * include 
 *
 * * SCATTER_MODE_ISOTROPIC -isotropic scattering (randvec),
 * * SCATTER_MODE_ANISOTROPIC - anisotropic scattering
 * (randwind) and
 * * SCATTER_MODE_THERMAL - thermally-broadened anisotropic scattering
 * (randwind_thermal_trapping).
 *
 ***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "atomic.h"
#include "python.h"


struct Cdf cdf_randwind_store[100];
CdfPtr cdf_randwind;
struct photon phot_randwind;



/*
randwind (xyz, lmn, north) is a routine to calculate a new direction
for the photon for specific presciptions regarding the anisotropy
of scattering.  To use a different prescription the only thing that
should have to be modified is vrandwind.

Here

 xyz is the position of the photon in cartesion coordinates
 north is the
 the direction of maxium velocity gradient in
	cylindrical coordinates

and
 lmn is the new direction of the photon in cartesian coordinates
	which is returned is returned.

02feb	ksl	Modified jumps in randwind so evenly spaced.
02feb	ksl	Mad major modifications to the way in which
		the direction was projected back to the
		cartesian frame, and also to the generation
		of the photon direction in the azimuthal
		direction.  The routine was seriously
		screwed up.
02may	ksl	Made modifications to generalize anisowind
		and to incorporate changes associated
		with creating a new (constant source function)
		based scattering function.  The routines are
		now much more integrated in the sense that
		only one routine vrandwind needs to be
		changed to incorporate a new scattering
		model (as long as the model is axisymmetric)
		about some direction..presumable the direction
		of tau_min.  N.B. There are some dependences
		on the order in which the routines are called,
		i.e. you must call randwind before vrandwind,
		since randwind initializes some variables that
		one needs in vrandwind.
02may	ksl	Changed calls to avoid having to pass tau_min
		to the routine.  Hopefully this will help
		eliminate questions of what has been calculated
		when, or at least isolate it to the routines
		in this file.
02may	ksl	Added a check to prevent tau from going too low
		when calculating pdf_rand.  ?? It might be more
		computationally efficient to switch to isotropic
		in this case ??
15aug	ksl	Minor changes for multiple
17jun	nsh - changed references from PDFs to CDFs - because
		this is what they actually are!
*/

#define  TAU_TOP   10.0         /* tau above which we assume the angular distribution function
                                   no longer changes significantly, e.g above which the surface behaves
                                   like a solid surface */
#define  TAU_BOT    0.01


double tau_randwind = -1000.;



/**********************************************************/
/**
 * @brief      calculate a new direction
 * for the photon which is being scattered in the wind
 * in the SCATTER_MODE_ANISOTROPIC scattering mode.
 *
 * @param [in] PhotPtr  p   The photon (at the position it is scattering, and the resonance causing the scatter indicated)
 * @param [out] double  lmn[3]   The new photon direction
 * @param [in] double  north[3]   the direction of maximum velocity gradient in the wind cell where the photon is located
 * @return     Normally returns 0, -1 if nres for the photon is not positive (which would indicate this was not a resonance scatter)
 *
 * @details
 * This routine is called whenever one wants to generate a new
 * direction in  the SCATTER_MODE_ANISOTROPIC scattering mode.
 *
 *
 * Given the photon position and the resonance,
 * the routine * first calcuates the optical
 * depth in the direction of the maximum dv_ds for that
 * cell.
 *
 * The routine then calls, make_cdf_randwind to create a
 * comulative distribution function for scattering of as
 * a function of the angle from the direction of minimum
 * optical depth.
 *
 * Using this cdf, it randomly selects a direction for
 * the photon.
 *
 * These calculations are all carried out as if the
 * photon were in the x-z plane, because that is where
 * the velocity gradient is defined.  As a result
 * one then has to trasform the new photon direction
 * back to the location of the photon
 *
 * ### Notes ###
 * To use a different prescription the only thing that
 * should have to be modified is vrandwind.
 *
 * It might have been more logical to have transformed
 * the direction of the maximum velcocity gradient
 * to the position of the photon.
 *
 *
 **********************************************************/

int
randwind (p, lmn, north)
     PhotPtr p;
     double lmn[3], north[3];
{

  double xyz[3];
  double n;                     /* the individual direction cosines in the rotated frame */
  double q;
  double phi;
  double xlmn[3], dummy[3];
  double north_xyz[3];
  struct basis nbasis;
  int k;
  double tau;


  if (p->nres < 0)
  {
    Error ("randwind: Cannot calculate line scattering without nres %d of line\n", p->nres);
    return (-1);
  }

/* Get the position of the photon in the appropriate domain, and calculate
 * the optical depth in the direction of the maximum velocity gradient,
 * which will be the direction of minimum optical depth */

  k = where_in_grid (wmain[p->grid].ndom, p->x);

  tau = sobolev (&wmain[k], p->x, -1., lin_ptr[p->nres], wmain[k].dvds_max);

  stuff_v (p->x, xyz);

  /* 180415 - ksl - There is a lot of belt and suspenders going on here.  In fact,
   * make_cdf_randwind whould be quite fast as it just selects
   * sets a pointer for cdf_rand_wind based on tau
   */

  if (tau > TAU_TOP)
    tau = TAU_TOP;
  if (tau < TAU_BOT)
    tau = TAU_BOT;

  if (fabs (tau - tau_randwind) > 0.01)
  {                             // Create or select the appropriate cdf
    make_cdf_randwind (tau);
    stuff_phot (p, &phot_randwind);
  }

  /* Having selected the appropriate cdf find the polar (theta direction
   * for the scatter
   */

  xlmn[0] = n = cdf_get_rand (cdf_randwind);
  if (sane_check (n))
  {
    Error ("randwind:sane_check of cdf_get_rand returned %f\n", n);
    xlmn[0] = n = cdf_get_rand (cdf_randwind);
  }

  q = sqrt (1. - n * n);

  /* Now get the azimuthal direction for the scattered photon */

  phi = 2. * PI * random_number (0.0, 1.0);

  xlmn[1] = q * cos (phi);
  xlmn[2] = q * sin (phi);


/* So at this point we have the direction cosines of the photon in a
   coordinate system in which "north" is along the FIRST axis. "north"
   as presented is a vector relative to the basis vectors for the cylindrical
   coordinate system.  So
   we need to express "north" in cartesian coordinates.
*/


  project_from_cyl_xyz (xyz, north, north_xyz);

/* Given north in cartesion coordinates, we can now follow the
same procedure as in randvcos, namely, create a basis in which
north is along the FIRST axis.  This is the basis in which we
have calculated the photon direction.  The other two directions
don't really matter, but to create a basis we need to ensure
that the two axes given in create_basis are not parallel to
one another.  */

  cross (north_xyz, z_axis, dummy);
  if (length (dummy) > 0.0)
    create_basis (north_xyz, z_axis, &nbasis);
  else
    create_basis (north_xyz, x_axis, &nbasis);

/* Now project from the coordinate system described by nbasis to
the cartesian frame */

  project_from (&nbasis, xlmn, lmn);

  if (sane_check (lmn[0]) || sane_check (lmn[1]) || sane_check (lmn[2]))
  {
    Error ("Randwind:sane_check NAN problem lmn %f %f %f\n", lmn[0], lmn[1], lmn[2]);
  }

  return (0);
}


#define VRANDWIND_FLOOR 1e-5
#define VRANDWIND_TAU 11.51


/**********************************************************/
/**
 * @brief      the function that is used to generate the cdf for anisotropic scattering.
 *
 * @param [out] double  x   The direction cosine with respect to the direction
 * of tau_min
 *
 * @return     A number proportional to the probability that a photon will
 * scatter in a direction with a particular direction cosine
 *
 * @details
 * The routine calculates
 * a function that is proportional to the probability
 * density.  dP/dcos(theta) for the photon.
 *
 * ### Notes ###
 * This routine is used by cdf_gen_from_func (via make_cdf_randwind, when that routine is
 * called from randwind)
 *
 **********************************************************/

double
vrandwind (x)
     double x;
{
  double abscos, z;

  abscos = fabs (x);
  if (abscos == 0.0)            // Then at 90 degrees
  {
    return (VRANDWIND_FLOOR);
  }


/* The next lines generate a model for scattering in which
   the source function is assumed to be constant through the
   scattering region.  It is also related to escape probs.  */

  z = (1. - exp (-tau_randwind / abscos)) * abscos;

/* These generate the probability density for the Eddington approximation
    a = 0.5;
    b = 1.5;
    z = x * (a * (1. + b * x));


   The next line would generate a completely isotropic distribution
   z=0.1;

      Completely isotropic
*/
  if (sane_check (z))
  {
    Error ("vrandwind:sane_check tau_randwind %f x %f\n", tau_randwind, x);
    return (0.0);
  }


  return (z);


}


#define REWEIGHTWIND_TAU_MAX 100.
int reweightwind_init = 1;      //TRUE to start
double reweightwind_zmax;


/**********************************************************/
/**
 * @brief      calculates the weight of a photon in the wind that
 * is forced to scatter along a specific line of sight.
 *
 * @param [in, out] PhotPtr  p   The photon being scattered
 * @return     The number corresponding to the reweighting is returned.  The
 * revised weight is also stored in p
 *
 * @details
 * p must contain the position and NEW direction for the photon
 * being scattered.  In general, p must also be up to date w. r. t.
 * the number of the resonance that generated the photon.
 *
 *
 * ### Notes ###
 * This routine is called when extracting photons along a specific
 * line of sight.
 *
 * wmain.dv_ds, and wmain.lmn are maximum value of dvds, and
 * the direction in a cell (that is
 * in the positive xz plane)
 *
 *
 **********************************************************/

double
reweightwind (p)
     PhotPtr p;
{
  double ctheta, x, z;
  double north_xyz[3];
  int project_from_cyl_xyz ();
  double xyz[3], lmn[3], north[3], tau;
  double delta[3];
  int k;

  vsub (p->x, phot_randwind.x, delta);

/* Idea here is that if photon has moved from position where it was you must recalc.
 * Since typically we extract photons in a range of inclinations, this routine
 * is called multiple times for the same photon when it scatters.  The if statement
 * here is intended to capture this, and only to recalculate the cdf if the photon
 * has moved from where it was previously.
 *
 */

/** @bug 180414 - ksl - It is not obvious that using DFUDGE as the way to choose whether
 * to recalculate the cdf for scatter is correct.  One would
 * obtain the wrong answer if the photon moved a small distance and hit a different
 * resonance.
 */
  k = where_in_grid (wmain[p->grid].ndom, p->x);

  if ((x = length (delta)) > DFUDGE)
  {
    tau = sobolev (&wmain[k], p->x, -1., lin_ptr[p->nres], wmain[k].dvds_max);
    make_cdf_randwind (tau);    // Needed for the normalization
    stuff_phot (p, &phot_randwind);
  }
  else
    tau = 0.0;

  stuff_v (p->x, xyz);
  stuff_v (p->lmn, lmn);
  stuff_v (wmain[k].lmn, north);

/* We need to know the cos of the angle between the direction
of maximum velocity gradient and the new direction of
the photon.  north is a vector expressed in terms of
a cylindrical coordinate system (at the point xyz) so one
way to proceed is to project north to xyz coords which is
what we do here */

  project_from_cyl_xyz (xyz, north, north_xyz);

  ctheta = fabs (dot (lmn, north_xyz));


//
/**
 * @bug This is a note related to an XXX comoment which read as follows:
 * Factor of 2 needed because the interval is from -1 to 1
 * XXX It's definitely needed to make a uniform distribution work
 * but is this reason really right
 *
 * 180415 - ksl - My investigation of this indicates that the cdf was generated
 * from -1 to 1, and futhermore that cdf_gen_from_function which is used to
 * create the cdf always sets the norm to 1.  So it is not clear to me why
 * this is the reweghting.  It is also nto clear to me why the reweighting to
 * be less than a factor of 2, which is what if stateement requries.  Basically,
 * I do not understand what is going on here at all.
 */

  z = 2. / cdf_randwind->norm;
  x = vrandwind (ctheta) * z;

  if (sane_check (x) || x > 2.0)
  {
    Error ("Reweightwind:sane_check x %f tau %f ctheta %f z %e cdf_randwind->norm %f\n", x, tau, ctheta, z, cdf_randwind->norm);
    x = 2.0;
  }

  p->w *= x;
  return (x);

}




/**********************************************************/
/**
 * @brief      Generate the cumulative distribution functions
 * needed for selecting a direction for a scattered photon, and
 * then select the cdf needed ro a particulat tau
 *
 *
 * @param [in] double  tau   The optical depth in the direction
 * of the maximum velocity gradient associated with
 * the scattering event
 *
 * @return     Always returns 0
 *
 * @details
 * The first time this routine is entered it generates a
 * set of cdfs correspoding to various values of tau which
 * are logrithmically spaced between 0.01 and 10.  These
 * cddfs are stored in an array of cdf structures
 *
 * It then selects the closest cdf and sets cdf_rand_wind
 * to point to this particular cdf
 *
 * On subsequent calls, the routine merely chooses which cdf
 * is most appropriate for a given tau
 *
 * ### Notes ###
 *
 * The probability densities are calculated in vrandwind
 *
 * The basic idea of this routine is to calculate a series
 * of cdfs the first time the routine is called that are
 * spaced in a sensible fashion so that one does not have
 * to regenerate the cdfs everytime one wants to calulate
 * the probability a photon will be scattered in certain
 * direction.
 *
 * On subsequent calls the routine simply chooses which
 * of the cdfs to use by setting cdf_randwind to one of
 * the precalculate cdfs.
 *
 *
 **********************************************************/

int init_make_cdf_randwind = 1;
int make_cdf_randwind_njumps;
double make_cdf_randwind_jumps[180];
#define LOGTAUMIN -2.
#define LOGTAUMAX 1.
double cdf_randwind_dlogtau;


int
make_cdf_randwind (tau)
     double tau;
{
  int jj;
  int echeck;
//OLD  int cdf_gen_from_func ();
//OLD  double vrandwind ();
  double xtau;
//OLD  double log10 ();

/* Initalize jumps the first time routine is called */

  if (init_make_cdf_randwind)
  {
    make_cdf_randwind_njumps = 0;
    for (jj = -88; jj <= 88; jj += 2)
    {
      make_cdf_randwind_jumps[make_cdf_randwind_njumps] = sin (jj / 57.29578);
      make_cdf_randwind_njumps++;
    }
    cdf_randwind_dlogtau = (LOGTAUMAX - LOGTAUMIN) / 99.;
    for (jj = 0; jj < 100; jj++)
    {
      xtau = pow (10., LOGTAUMIN + cdf_randwind_dlogtau * jj);
      tau_randwind = xtau;      /* This is passed to vrandwind by an external variable */
      if ((echeck =
           cdf_gen_from_func (&cdf_randwind_store[jj], &vrandwind, -1.0, 1.0, make_cdf_randwind_njumps, make_cdf_randwind_jumps)) != 0)
      {
        Error ("make_cdf_randwind: return from cdf_gen_from_func %d\n", echeck);
      }
    }
    init_make_cdf_randwind = 0;
  }

  if (sane_check (tau))
  {
    Error ("make_cdf_randwind:sane_check Need proper tau (%e) to make cdf_randwind\n", tau);
    tau = 10.;                  // Forces something close to isotropic
  }

  jj = (log10 (tau) - LOGTAUMIN) / cdf_randwind_dlogtau + 0.5;

  if (jj < 0)
    jj = 0;
  if (jj > 99)
    jj = 99;

  cdf_randwind = &cdf_randwind_store[jj];
  tau_randwind = tau;           // This is passed to vrandwind as an external variable

  return (0);
}


//OLD /***************************************************************
//OLD
//OLD                       University of Southampton
//OLD
//OLD Synopsis:
//OLD   randwind_thermal_trapping is the routine which chooses
//OLD   a new anisotropic direction in geo.scatter_mode = SCATTER_MODE_THEMAL
//OLD
//OLD Arguments:
//OLD
//OLD
//OLD Returns:
//OLD   0 for success. Also modifies the photon ptr p
//OLD   to reflect new direction (p->lmn), and nnscat, which
//OLD   should be copied to the phoiton structure after calling
//OLD   this routine.
//OLD
//OLD Description:
//OLD   This routine uses a rejection method to choose a direction
//OLD   so that the probability distribution of directions generated
//OLD   reflects the probability of escape along each direction in accordance
//OLD   with the sobolev optical depth.
//OLD
//OLD Notes:
//OLD
//OLD History:
//OLD   1406  Moved code here from photo_gen_matom and scatter to avoid
//OLD         duplication
//OLD
//OLD
//OLD ****************************************************************/



/**********************************************************/
/**
 * @brief      is the routine which chooses
 *   a new anisotropic direction in geo.scatter_mode = SCATTER_MODE_THERMAL
 *
 * @param [in,out] PhotPtr  p   The photon being scattered
 * @param [out] int *  nnscat   The number of times the phton
 * scattered internally before escaping the local scattering region
 *
 * @return     0 for success. Also modifies the photon ptr p
 *   to reflect new direction (p->lmn), and nnscat, which
 *   should be copied to the phoiton structure after calling
 *   this routine.
 *
 * @details
 * This routine uses a rejection method to choose a direction
 * so that the probability distribution of directions generated
 * reflects the probability of escape along each direction in accordance
 * with the sobolev optical depth.
 *
 * ### Notes ###
 * The resonance that caused the scatter must be stored in the
 * photon bundle.
 *
 * The name of the routine is something of a misnomer. Pure
 * sobolev optical depths are used.  The temperature in the
 * cell does not come into the calculation.
 *
 * Unlike randwind, this routine does not explicitly need to to
 * make a transformation from the xz plane to the location of the
 * photon.  This is because dvds is calculated directly using dwind_ds
 *
 **********************************************************/

int
randwind_thermal_trapping (p, nnscat)
     PhotPtr p;
     int *nnscat;
{
  double tau_norm, p_norm;
  double tau, dvds, z, ztest;
  double z_prime[3];
  WindPtr one;

  /* find the wind pointer for the photon */
  one = &wmain[p->grid];

  /* we want to normalise our rejection method by the escape
     probability along the vector of maximum velocity gradient.
     First find the sobolev optical depth along that vector */
  tau_norm = sobolev (one, p->x, -1.0, lin_ptr[p->nres], one->dvds_max);

  /* then turn into a probability. Note that we take account of
     this in trans_phot before calling extract */
  p_norm = p_escape_from_tau (tau_norm);

  /* Throw error if p_norm is 0 */
  if (p_norm <= 0)
    Error ("randwind_thermal_trapping: p_norm is %8.4e in cell %i", p_norm, one->nplasma);

  ztest = 1.0;
  z = 0.0;

  /* JM 1406 -- we increment nnscat here, and it is recorded in the photon
     structure. This is done because we actuall have to multiply the photon weight
     by 1/mean escape probability- which is nnscat. this is done in trans_phot.c
     before extract is called.
   */
  *nnscat = *nnscat - 1;

  /* rejection method loop, which chooses direction and also calculated nnscat */
  while (ztest > z)
  {
    *nnscat = *nnscat + 1;
    randvec (z_prime, 1.0);     /* Get a new direction for the photon (isotropic */
    stuff_v (z_prime, p->lmn);  // copy to photon pointer


    /* generate random number, normalised by p_norm with a 1.2 for 20%
       safety net (as dvds_max is worked out with a sample of directions) */
    ztest = random_number (0.0, 1.0) * p_norm;

    dvds = dvwind_ds (p);
    tau = sobolev (one, p->x, -1.0, lin_ptr[p->nres], dvds);

    z = p_escape_from_tau (tau);        /* probability to see if it escapes in that direction */
  }

  return (0);
}
