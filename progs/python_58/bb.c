
/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	
double planck(t,freqmin,freqmax)	returns the frequency of a "pseudo-randomly generated photon. 
double t,freqmin,freqmax;		temperature of BB, and frequency limits
  
 					The frequency is pseudo random in the following limited sense.  
					The photons are returned are weighted so that the energy distribution 
 					of a planck function is approximately reproduced. Within an energy 
 					bin the photon frequency is uniform.  

double
emittance_bb (freqmin, freqmax, t)
     double freqmin, freqmax, t;
					Calculate the emittance of a bb between freqmin and freqmax

double integ_planck_d(alphamin,alphamax)
double alphamin,alphamax;

					returns the integral of the dimensionless blackbody function
					between alphamin and alphamax.  

					Because we cannot afford to keep reintegrating this integral, 
					the following procedure is used.  The first time the function is
					called the integ_planck[] is filled.  integ_plank[n]  contains the 
					integral of the dimensionless planck function from  ALPHAMIN to 
					ALPHAMAX.  Therefore if one wants to obtain the integral of the 
					dimensionless bb function, one simply interpolates this array


init_integ_planck_d() 			performs the integration of the dimensionless bb function
					in order to fill the arry integ_planck.  This routine needs only to be
					 called onetime. The actuall integration is done with the numerical 
					 recipes routine "qromb" .

Arguments:		

Returns:
 
Description:	
		The way this works currently is a little obscure.  The primary routines are emittance_bb or planck.  
		They call the other routines

		The first call to either of these routines (surely emittance_bb) results in a call to integ_planck_d, 
		which in turn calls integ_planck_init.  This initialization routine turn populates the array integ_planck
		, which contains the integral of the bb function in an array.

		bb_emittance continues to access the array integ_plank through integ_planck_d every 
		future time is is called.

		planck does the same thing albeit more indirectly. It sets up a pdf each time new frequency 
		limits are placed on it.  planck therefore really uses the
		pdf.

		
Notes:
		It seems possible that some of these routines could be combined in a clever way
		
		Bug/Problem:  04mar  ksl
			The program saves time by only integrating dimensionless bb function one 
			time.  Then you use pdf_limit to handle frequency limiting.  But we have introduced
			photon banding to force the program to generate photons in certain portions of the
			spectrum, where they are not expected on average.  This is a problem, because the program is
			written so that all pdfs have the same dimension.  Expanding this number without modifying
			the pdf routines will cause the entire program to slow down.  On the other hand, increasing
			the number of fixed points in the pdf, with bb_set, will limit the quality of the BB spectra
			produced in the region where most of the energy flux is.  An additional problem, is that
			no matter what one does, when you limit the range to the very ends of the pdf, then one
			gets a very poor represntation of the cdf.  This problem appeared in some of the disks models
			that I tried, when I was producing IR spectra.

			Redesigning all of the pdf generating software is undesirable, since it would require 
			changes of many routines.   

			One possibility that might work would be to create multiple dimensionless pdfs, one for
			the nominal case 0-30, one for the high end, say 20-100, and one for the low end, say
			0-1.  The problem really is that on the high end particularly the probablilty distribution
			is hugely changing over the wavelength band.  

			Another possibility would be to use pdf_gen_from_array, since this can grab an array
			of any reasonable size (currently 4000).  There is a good bit of work to be done, and
			so ideally one would still not generate a new pdf every time.    

			The problem that we are finding with planck_d also must existe in integ_planck_d at some 
			level, since this is also called exteranally, and uses ALPHAMAX to limit the integration.
			
			At present, when we get the error I was trying to solve one simply returns the minimum
			frequency, on the assertion that the probablility density is dropping preciptiouls with
			frequency.  Maybe the correct thing to do with that particular problem is simply to
			make an approximation of the frequency width and use that.  Something similar might
			be a better approch for the low frequency limit as well.

			It turns out that in the restricted case where slpha=nu/kT is constrained to be large, then
			then the randomly selected photon frequency is given by
				alpha = alpha_min + ln (1-rand) === alpha_min + ln (rand)
			where rand is a uniforly generated random number.

			Similarly in the low frequency limit 
			        alpha= alphamax * (1-rand)**1/3 ==== alpha_max * rand**1/3

			This is what is now implemented.

			0608--57h -- There was a bug in the routine that Stuart had implemented 
			pdf_get_rand_limit_bbir (&pdf_bb); to handle the
			case where things were in the low temperature limit, which caused the routine to give
			photons outside the allowed frequency range. He was trying to correct the way
			the program handles low enery photons, but the problem also occurs at high enegies.
			The answer to the problem really is to break the geneation into three regions, 
			a low, a reasonble and a high frequency regime. 

			I think the right way to do this is as follows.  Decide where the boundaries between
			the three rebimes are.  Create a pdf which the first inteval encompasses the
			entire low energy regime, and the last interval the entire high energy regime.
			Generate a random number.  If this does not end up in either of the edge
			bins then take a standard linear interpolation.  If it does however, do a second
			calculation for the edge.

			06aug -- At least a partial fix to this problem is now implemented in the routines
			below.  There is additional work to be done for the high frequency limit.  

History:
 	97	ksl	The various portions of this file were written in 1997 as part of the python effort
 	98mar	ksl	Combined the routines here from bb_gen.c and planck.c to reduce the total number of
 			files.
 	98mar23	ksl	planck(t,freqmin,freqmax) modified to use new pdf routines
 	98apr25	ksl	integ_planck_d modified to fix an error when the ranges of alphamin and alphamax
 			were out of bounds
	98jul27	ksl	planck moddified to force denser selection of gridpoints for small and high freq relative
			to kt
	01dec	ksl	Increased number of points at high end as part of exercise to force selection
			of some high frequency photons (python_40)
	04mar	ksl	Modified to provide a better treatment when the frequency range is
			above or below the frequency limits of the cdf that was calculated 
			initially.
	05jul	ksl	56d -- Relegated writing the pdf file to a case where you were debugging
	06aug	kls	57h -- Recommented and reorganized the file, so that could complete a partial
			fix to a problem Sturat had identified.  The partial fix was generating 
			photons outsde the min and maximum frequency

**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "atomic.h"
#include "python.h"

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis: planck returns the frequency of a "pseudo-randomly generated photon.  The
   frequency is pseudo random in the following limited sense.  The photons are returned
   are weighted so that the energy distribution of a planck function is approximately
   reproduced. Within an energy bin the photon frequency is uniform.  
	

Arguments:		

Returns:
 
Description:	

Notes:

`	06aug - ksl -- I am fairly sure that one can eliminate the sections
	in this routine that have to do with high and low temperuature limits
	in favor of makeing a few additional modifications to the new routine
	bb_pdf_get_rand_limit below.  

	06sep - ksl -- I really think what should be done in this program is to
	rewrite so that one generates a random number at the start.  If this 
	random numbers leave you in a portion of the range that is the very high 
	or low frequency limit, you shoud use it to generate a frequency from
	that interval with one of the approsimations that we have.  If it leaves
	you in the main part of the interval you should use the pdf that was 
	initailly generated. This is not what is currently implemented.

	This would require determining where to set the high and low frequency limits.
	Preusmably what one owould do is find alphamin and alphamax where the low and
	high grequence limits occcured, e.g where the probability was less than say a 0.1
        or greater than 0.9.  One would consider these the boundaries.  Anytime the random
	number came up in these limtis one would use the low or the high fequeny nunbers.

	Note that alpha small is the low freqency limit, alpha large corresponds o
	high freqencies relative to the T.
		
History:
	06aug	ksl	57h --Recommented, as moved to fix a problem
			that was returning photons that were occassionally outside
			of the frequency limits.  The problem had been introduced
			to handle the YSO case, where Stuart had found the uniform
			distribution used to interpolate between frequency points
			annoyoing
	06aug	ksl	Eliminated call to Stuart's attempt to better handle
			the low frequency limit, in favor of a new routine
			that I hope fixes the problem of the frequency
			limits.  
	06sep	ksl	Hopefully correct the remaining problems with this,
			but there is a significant issue with all of our pdfs after 
			we have created pdf arrays, we simply create a uniform 
			distribution within an interval.
	06sep	ksl	57i -- Have simplfied this to elimiante for now any lo or hi
			freque4ncy special cases
**************************************************************/

#define ALPHAMIN 0.4		// Region below which we will use a low frequency approximation
#define ALPHAMAX 30.		// Region above which we will use a high frequency approximation
#define ALPHABIG 100.		//  Region over which can maximmally integrate the Planck function
#define NMAX 		1000

int ninit_planck = 0;
struct Pdf pdf_bb;

double old_t = 0;
double old_freqmin = 0;
double old_freqmax = 0;
double alphamin, alphamax;
double cdf_bb_lo, cdf_bb_hi, cdf_bb_tot;	// The precise boundaries in the the bb cdf 
double cdf_bb_ylo, cdf_bb_yhi;	// The places in the CDF defined by freqmin & freqmax
double lo_freq_alphamin, lo_freq_alphamax, hi_freq_alphamin, hi_freq_alphamax;	//  the limits to use for the low and high frequency values

// bb_set is thae array that pdf_gen_from_func uses to esablish the 
// specific points in the cdf of the dimensionless bb function.
double bb_set[] = {
  0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
  10.3, 11.3, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
  19., 20., 21., 22., 23., 24., 25., 26., 27., 28., 29.
};


int error_bb_hi = 0;
int error_bb_lo = 0;

double
planck (t, freqmin, freqmax)
     double t, freqmin, freqmax;
{
  FILE *fopen ();
  double freq, alpha, y;
  double planck_d (), pdf_get_rand_limit ();
  double get_rand_pow ();
  int pdf_gen_from_func (), pdf_to_file (), echeck;
  int pdf_limit ();


  /*First time through create the array containing the proper boundaries for the integral of the BB function,
     Note calling pdf_gen_from func also difines ylo and yhi */
  if (ninit_planck == 0)
    {				/* First time through p_alpha must be initialized */
      if ((echeck =
	   pdf_gen_from_func (&pdf_bb, &planck_d, ALPHAMIN, ALPHAMAX, 29,
			      bb_set)) != 0)
	{
	  Error ("Planck: on return from pdf_gen_from_func %d\n", echeck);
	}
      /* We need the integral of the bb function outside of the regions of interest as well */


      cdf_bb_tot = qromb (planck_d, 0, ALPHABIG, 1e-8);
      cdf_bb_lo = qromb (planck_d, 0, ALPHAMIN, 1e-8) / cdf_bb_tot;	//position in the full cdf of low frequcny boundary
      cdf_bb_hi = 1. - qromb (planck_d, ALPHAMAX, ALPHABIG, 1e-8) / cdf_bb_tot;	//postion in fhe full hi frequcny boundary

//      pdf_to_file (&pdf_bb, "pdf.out");
      ninit_planck++;

    }


/* If temperatures or frequencies have changed limit the region of the pdf in which we are interested */

  if (t != old_t || freqmin != old_freqmin || freqmax != old_freqmax)
    {

      alphamin = H * freqmin / (BOLTZMANN * t);
      alphamax = H * freqmax / (BOLTZMANN * t);

      old_t = t;
      old_freqmin = freqmin;
      old_freqmax = freqmax;

//      freq = BOLTZMANN * t / H * ALPHAMIN;
//      Log ("Boundary %f\n", C * 1e8 / freq);

      cdf_bb_ylo = cdf_bb_yhi = 1.0;
      if (alphamin < ALPHABIG)
	{
	  cdf_bb_ylo = qromb (planck_d, 0, alphamin, 1e-8) / cdf_bb_tot;	//position in the full cdf of current low frequcny boundary
	  if (cdf_bb_ylo > 1.0)
	    cdf_bb_ylo = 1.0;
	}
      if (alphamax < ALPHABIG)
	{
	  cdf_bb_yhi = qromb (planck_d, 0, alphamax, 1e-8) / cdf_bb_tot;	//position in the full cdf of currnt hi frequcny boundary
	  if (cdf_bb_yhi > 1.0)
	    cdf_bb_yhi = 1.0;
	}

/* I still need to determine what all of the limits are set */

/* These variables are not always used */
      lo_freq_alphamin = alphamin;	// Never used if 
      lo_freq_alphamax = alphamax;
      if (lo_freq_alphamax > ALPHAMIN)
	lo_freq_alphamax = ALPHAMIN;

      hi_freq_alphamax = alphamax;
      hi_freq_alphamin = alphamin;
      if (hi_freq_alphamin < ALPHAMAX)
	hi_freq_alphamin = ALPHAMAX;


      if (alphamin < ALPHAMAX && alphamax > ALPHAMIN)
	{
	  pdf_limit (&pdf_bb, alphamin, alphamax);
	}




    }


  y = rand () / (MAXRAND);

  y = cdf_bb_ylo * (1. - y) + cdf_bb_yhi * y;	// y is now in an allowd place in the cdf

/* There are 3 cases to worry about
	The case were everything is in the low frequency limit
	The case where everything is in the normal limit
	The case where some photons are in the low regime and some are
	in the normal regime
*/

  if (y <= cdf_bb_lo || alphamax < ALPHAMIN)
    {
      alpha = get_rand_pow (lo_freq_alphamin, lo_freq_alphamax, 2.);
    }
  else if (y >= cdf_bb_hi || alphamin > ALPHAMAX)
    {
      alpha = get_rand_exp (hi_freq_alphamin, hi_freq_alphamax);
    }
  else
    {
      alpha = pdf_get_rand_limit (&pdf_bb);
    }

  freq = BOLTZMANN * t / H * alpha;
  if (freq < freqmin || freqmax < freq)
    {
      Error ("planck: freq %g out of range %g %g\n", freq, freqmin, freqmax);
    }
  return (freq);
}



/***********************************************************
	Space Telescope Science Institute

 Synopsis: get_rand_pow obtains a random number between x1 and x2 for a power law densitity distribution with
	index alpha
	

Arguments:		

Returns:
 
Description:	
	this array. 
Notes:
		
History:
	06sep	ksl	Coded       
**************************************************************/
double
get_rand_pow (x1, x2, alpha)
     double x1, x2, alpha;
{
  double r;
  double a;

  r = rand () / MAXRAND;

  x1 = pow (x1, alpha + 1.);
  x2 = pow (x2, alpha + 1.);

  a = (1. - r) * x1 + r * x2;

  a = pow (a, 1. / (alpha + 1.));
  return (a);
}


/***********************************************************
	Space Telescope Science Institute

 Synopsis: get_rand_exp wobtains a random number between x1 and x2 for a exp law densitity distribution with
	index alpha. 
	

Arguments:		

Returns:
 
Description:	
	this array. 
Notes:
		
History:
	06ocd	ksl	Coded       
**************************************************************/
double
get_rand_exp (alpha_min, alpha_max)
     double alpha_min, alpha_max;
{
  double r;
  double x1, x2;
  double a;

  r = rand () / MAXRAND;

  x1 = exp (-alpha_min);
  x2 = exp (-alpha_max);

  a = (1. - r) * x1 + r * x2;

  a = log (a);
  return (-a);
}

/***********************************************************
	Space Telescope Science Institute

 Synopsis: integ_plank_d(alphamin,alphamax) returns the integral of the dimensionless blackbody function
   between alphamin and alphamax.  
	

Arguments:		

Returns:
 
Description:	
	To save computing time, the routine actually accesses an array that contains the
	the integral of the bb function from 0 to x.

   	The first time the function is
	called the integ_planck[] is filled.  integ_plank[n]  contains the integral of the
	dimensionless planck function from  ALPHAMIN to ALPHAMAX.  Therefore if one
	wants to obtain the integral of the dimensionless bb function, one simply interpolates
	this array. 
Notes:
		
History:
	06aug	ksl	Recommented
**************************************************************/


double integ_planck[NMAX + 1];
int i_integ_planck_d = 0;
double
integ_planck_d (alphamin, alphamax)
     double alphamin, alphamax;
{
  double x, z1, z2;
  int n;
  int init_integ_planck_d ();
  if (i_integ_planck_d == 0)
    {				/*First time through integ_planck must be defined */
      init_integ_planck_d ();
      i_integ_planck_d++;
    }

  x = (alphamin - ALPHAMIN) / (ALPHAMAX - ALPHAMIN) * NMAX;
  if (x <= 0.0)
    z1 = 0.0;
  else if (x >= (NMAX))
    {
      return (0.0);		/* Because the minimum frequency is too high */
    }
  else
    {
      n = x;
      x -= n;
      z1 = integ_planck[n] * (1. - x) + integ_planck[n + 1] * x;
    }

  x = (alphamax - ALPHAMIN) / (ALPHAMAX - ALPHAMIN) * NMAX;
  if (x < 0.0)
    {
      return (0.0);		/* Because the maximum frequency is too low */
    }
  else if (x >= (NMAX))
    {
      z2 = integ_planck[NMAX];
    }
  else
    {
      n = x;
      x -= n;
      z2 = integ_planck[n] * (1. - x) + integ_planck[n + 1] * x;
    }

  return (z2 - z1);
}


/***********************************************************
       Space Telescope Science Institute

 Synopsis: init_integ_planck_d calulates integrals of the dimensionaless bb function and
 	stores it in the array integ_planck.  Each value of integ_plank contains
	the integration from 0 to x
	

Arguments:		

Returns:
 
Description:	
   This routine needs to be called once, since ALPHAMIN and ALPHAMAX are
   hardcoded. 

   The actual integration is done with the numerical recipes routine "qromb" 

Notes:
		
History:
	97+	ksl	Originally coded
	06aug	ksl	Revised comments

**************************************************************/

int
init_integ_planck_d ()
{
  double x;
  double planck_d (), qromb ();
  int n;
  integ_planck[0] = 0;
  for (n = 1; n <= NMAX; n++)
    {
      x = ALPHAMIN + n * (ALPHAMAX - ALPHAMIN) / NMAX;
// 1e-7 is the fractional accuracy in my modified version of qromb -- ksl
      integ_planck[n] = qromb (planck_d, 0.0, x, 1e-7);
    }

  return (0);
}


/***********************************************************
	Space Telescope Science Institute

 Synopsis:
"planck_d" is the dimensionless bb function.  The total emittance
   is related to planck_d as follows:
	

Arguments:		

Returns:
 
Description:	
   F_nu= 2*PI* (kT/h**3)*planck_d(h*freq/ k T)
Notes:
		
History:
	06aug	ksl	Recommented
**************************************************************/

#define EPSILON	1.e-6

double
planck_d (alpha)
     double alpha;
{
  double x;
  if (alpha < EPSILON || alpha > ALPHABIG)
    return (0);
  x = (alpha * alpha * alpha) / (exp (alpha) - 1);
  return (x);
}

// Calculate the emittance of a bb between freqmin and freqmax
// Should integrate to sigma 
double
emittance_bb (freqmin, freqmax, t)
     double freqmin, freqmax, t;
{
  double alphamin, alphamax, q1;
  double integ_planck_d ();
  q1 =
    2. * PI * (BOLTZMANN *
	       BOLTZMANN * BOLTZMANN * BOLTZMANN) / (H * H * H * C * C);
  alphamin = H * freqmin / (BOLTZMANN * t);
  alphamax = H * freqmax / (BOLTZMANN * t);
  return (q1 * t * t * t * t * integ_planck_d (alphamin, alphamax));
}

#undef NMAX
#undef ALPHAMIN
#undef ALPHAMAX
