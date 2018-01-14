/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
 
 	These are a series of routines which are designed to generate a random sampling
 	of a probability distribution using cumulative distribution functions.
 	
 	The main routines are
 	
		cdf_gen_from_func	(&cdf,&func,xmin,xmax,njumps,jump)		
			Generate a cumulative distrubution function (cdf) from a function
 	          
		cdf_gen_from_array(&cdf,x,y,n_xy,xmin,xmax,njumps,jump)		
			Generate a cumulative distribution function (cdf) from an
 			array of values

		cdf_get_rand(&cdf)				
			Generate a single sample from a cdf defined by either of the first
 			two routines 	
 	
 	Sometimes it is expensive computationally to recalculate the distribution function
 	every time you change the limits.   An example of this is the distribution function
 	for a bb where the distribution function is quite smooth (and can be represented in a 
 	dimensionless manner) but calculating the function and the distribution function
 	can be expensive.  The following two routines allow one to sample a part of a
 	pdf.  Note that pdf_gen-from_func or pdf_gen_from_array should have been
 	called previously!
 
		cdf_limit(&pdf,xmin,xmax)  sets limit1 and limit2 so that one can generate
 			distributions from a limited range within a previously created cdf.

		cdf_get_rand_limit(&pdf) will samples the cdf but limits the range of
			the return to lie between xmin and xmax from pdf_limit.

	
	There are some internal routines

		gen_array_from_func (func, xmin, xmax, pdfsteps) generates a cdf from a function
			on values linearly spaced between xmin and xmax.  pdfsteps gies the
			number of points.

 	There is also a simple diagnostic routine
		cdf_check(&cdf)

	which returns a nonzero number if certain checks of the cdf fail
	and a routine  

		cdf_to_file(&cdf,filename)				
 		  	to write a cdf structure to a file


								
Arguments for cdf_gen_from_func
	CdfPtr cdf;		A ptr to a cdf structure (see below)
	double (*func)()	The probablility density function to be integrated 
				to create the pdf
	double xmin,xmax;	The range over which the function is to be 
				integrated
	double jump[];		Specifix values of x, perhaps a Lyman edge, to be
				included in the pdf_structure.  These must be in
				ascending order but do not have to lie between xmin
				and xmax
	int njumps;		The size of jumps[]

Arguments for cdf_gen_from_array are same as above except
 	
	double x[],y[]		The one dimensional arrays containg the places x where the probability 
				density function density has been evaluated.  The array need not 
				be uniformly spaced. y is the value of the pdf at the point x 
	int n_xy		The size of x[] and y[]

Arguments for pdf_get_rand
	CdfPtr cdf		cdf (See below)	is a structure which contains the cumulative
				distribution function.

Arguments for cdf_limit
	double xmin,xmax	Here xmin and xmax will cause cdf[].limit1 and cdf[].limit2 to be set
				in such a way that when you call pdf_get_rand_limit the distribution
				will be sampled only between xmin and xmax.  Watch out though
				because if xmin and xmax lie outside of the original value of xmin
				and xmax then the distribution will not extend this far.							
Arguments for pdf_to_file
	PdfPtr cdf
	char filename[]		If you can't figure this out, you are in deep trouble!
				
Returns:
	All of the routines return 0 on successful completion, except cdf_get_rand which returns
	a random value between xmin and xmax.  Multiple calls to cdf_get_rand should regenerate
	the probability density function.
	
	If there is an error the routines usually send an error message to the screen and logfile
	and cause the program to stop!  Occasionally the program will send a nonzero return to
	indicate an error.
 
Description:

	The purpose of these routines is first to generate a cumulative distribution function
	and then allow one to sample that distribution function.  The procedure is either
	to call pdf_gen_from_func if you have a function whose probability distribution you
	want (as might be the case if you were generating bremsstrahlung photons and you knew
	the formula)  or cdf_gen_from_array (as might be the case if you were trying to
	generate Monte Carlo spectra from precalculated Kurucz models).
	
	The result of calling either of these routines will be to populate a structure of
	the form Cdf
		cdf->x[]  will contain NPDF values of x between xmin and xmax.  The values will be
			spaced so that it is almost but not quite true that if you generate a uniform
			number n between 0 and NPDF you could just read off the value pdf[n].x and
			that by doing this multiple times you could regenerate the distribution function
		cdf->y[]  will contain values between 0 and 1 and is the integral of the cumulative
			distribution function from xmin to pdf[].x
			
	Once a cdf has been generated, one samples the cdf by calls to cdf_get_rand which 
	creates a random number between 0 and 1, finds the elements in cdf which surround
	the random number and interpolates to return a value x. 
	
	It is possible to force specific values of x to appear in the pdf.  This is desirable
	if there are edges where the probability density changes rapidly, as for example
	at the Lyman edge in a stellar spectrum.  By using jumps properly you can assure
	that edges will be sharp in the regenerated distribution (since otherwse cdf_get_rand
	will linearly interpolate between two points on the opposite sides of the edge).

Notes:
	Must be compiled with log.c since it makes use of the "python" logging routines for
	errors.

	It would probably be desirable to modify this so that the size of the parallel 
	arrays being generated could be determined from inside the program.  In that case
	NPDF would have to become the maximum size. - NSH stis was tried in Jul 2017 - turns
	out to be very hard!


	Error??: For Kurucz models, The probability distribution beyond the He edge  
	is often very peaked at the long wavelength end and rapidly changing.  In this 
	instance, some bizarre results can be obtained because NPDF is not enough to 
	properly populate the pdf.  One could increase NPDF, but this would slow 
	everything down.  This is an another argument that the whole routine pdf_gen_from_array
	neeeds to be rewritten to avoid using NPDF altogether.  One way to do this
	would be simply to use the values of x and y in creating an inital pdf.  There
	is however a question of resampling. - NSH Jul 2017 - this has been addressed
	by removing the sampling aspect. Turns out with the faster machines we dont
	need to be so careful about making it easy to find the right point in a CDF.
	

History:
	97aug	ksl 	Initial versions of these and similar routines coded as part of python project.
	98mar	ksl	Recoded all of the pdf routines so that jumps were incorporated and so
			that they all use the pdf structure. 
	98mar25	ksl	Fixed a problem with gen_from_func which occured when there were
			discontinues changes in the pdf which resulted in multiple values of
			the pdf being associated with the same probability.  Also, fixed I had
			not done proper accounting when there were specific values at which
			I want the pdf recorded. 
	98jul27	ksl	Fixed a problem with renormalization inf pdf_gen_from_func for the case when
			one wanted the pdf at specific values.  Discovered, and have fixed a problem
			in pdf_gen_from_array and pdf_gen_from_func associated with requireing the
			pdf at specific values.  Note: while I have tested this for pdf_gen_from func
			I have not really tested it for pdf_gen_from_array and this needs to be done.!!
			Added a routine pdf_check which looks for errors in a pdf_array.
	98jul27 ksl	Added check to assure that if one does put in points at specific values
                        that the probability density was not zero over this interval which would
                        produce two y values that are exactly the same.	
	00dec30	ksl	Fixed small error to change order of truth calculation in a while statement that
			was causing the program to fail on Proga's alpha
	06sep	ksl	57i -- Previously, one assumed the probability density was effectively
			uniform between the places where the cumulative distribution function
			was calculated.  Now we attempt to make a correction for this.
	06oct	ksl	58 -- After a lot of testing, I decided to baseline the version
			of pdf that makes a first order correction for non-uniformity of
			the probability distribution within intervals where the cdf is
			calculated.  The version here corresponds to 57ib.  I still don't understand
			some of the timing issues associated with this, as there are circumstances
			where this more excact calculation seems to take less time to exeucute
			in python than the old version.  The old version should still be in the
			directory (and is called pdf_uniform.c)
	10oct	ksl	Changed the way in wich pdf_gen_from function works so it will change the
			the number of points that are used in the array pdf_array in situations
			where the binning is too course.  In the process, eliminated pdf_init.  
			It was only called by one routine.  
	17Jul   nsh Modified to make everything refer to CDFs as they should be, also changed
			cdf_gen_from_array to work on a supplied array - no jumps.

 
**************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "atomic.h"
#include "python.h"
#include "models.h"
#include <gsl/gsl_sort.h>
#include <gsl/gsl_interp.h>


/*  The structure is defined in python.h.  Here for reference only */
//#define NPDF 200

//typedef struct Pdf {
//      double x[NPDF+1];
//      double y[NPDF+1];
//      double d(NPDF + 1];      
//      double limit1,limit2;
//      double norm;           //The scaling factor which would renormalize the pdf
//      int npdf    
//} *PdfPtr,pdf_dummy;


#define PDFSTEPS 10000          // This is the initial value of PDFSTEPS
int pdf_steps_current;          // This is the value of pdfsteps at this point in time
int init_pdf = 0;
double *pdf_array;

/* Generate a pdf structure from a function.  

Description:

This routine stores values of the function in  cdf_array

	02feb	ksl	Modified to handle a problem in which the
			pdf was not monotonic because of the 
			interaction between jumps and desired
			values of the pdf.  This is probably more
			like pdf_gen_from_array.  At this point
			it is not entirely clear to me that 
			pdf_gen_from function should not simple
			call pdf_gen_from array...but that is
			not currently the case 
	10oct	ksl	Modified (and reorganized) in an attempt to calculate a 
			finer grid of points in pdf_array, in situations
			where array was so course that we did not
			adequately sample the CDF.

*/
int
cdf_gen_from_func (cdf, func, xmin, xmax, njumps, jump)
     CdfPtr cdf;
     double (*func) (double);
     double xmin, xmax;
     double jump[];
     int njumps;
{
  double xstep;
  double y;
  int j, m, mm, n;
  int njump_min, njump_max;
  int icheck, pdfsteps;
  int cdf_check (), calc_cdf_array ();
  double gen_array_from_func (), delta;

  njump_min = njump_max = 0;
  /* Check the input data before proceeding */
  if (xmax <= xmin)
  {
    Error ("pdf_gen_from_func: xmin %g <= xmax %g\n", xmin, xmax);
    exit (0);
  }

  if (njumps > 0)
  {
    for (j = 1; j < njumps; j++)
    {
      if (jump[j] <= jump[j - 1])
      {
        Error ("pdf_gen_from_func: jump[%d]=%g <=jump[%d]=%g out of order\n", j, jump[j], j - 1, jump[j - 1]);
        exit (0);
      }
    }

    njump_min = 0;
    while (njump_min < njumps && jump[njump_min] <= xmin)
      njump_min++;
    njump_max = 0;
    while (njump_max < njumps && jump[njump_max] < xmax)
      njump_max++;
    /* So at this point njump_min will point to the first jump which is betewen xmin and
       xmax or it will equal to njumps in which case there were no jumps which were greater
       than xmin.

       Similarly njump_max will be the point to the first jump above xmax or if there are
       no jumps above xmax, then it will be njumps. */

    njumps = njump_max - njump_min;
  }

  /* OK, all the input data seems OK */

  /* Construct what is effectively is the definite integral from xmin to x. Note
     however that currently pdf_array[0] corresponds awkwardly to the integral
     from xmin to xstep.  

     There code increases the size of the array if the steps appear so large that they
     will affect the digitization.  

     Problems might occur if there are very large jumps in the function. 

     10oct - ksl - Incorprated into the code in attempting to create a power law.
   */

  delta = 1.;
  n = 0;
  pdfsteps = PDFSTEPS;

  while (n < 3)
  {
    delta = gen_array_from_func (func, xmin, xmax, pdfsteps);
    if (delta < 0.1 / FUNC_CDF)
      break;
    pdfsteps *= 10;
    n = n + 1;
  }


  xstep = (xmax - xmin) / pdfsteps;

/* So at this point pdf_array contains an unnormalized version of the CDF for the function */


  cdf->x[0] = xmin;
  cdf->y[0] = 0;

  n = 0;                        //This is the position in pdf_array
  mm = 1;                       //This is the index to a desired value of y, with no jumps
  j = njump_min;                //This refers to the jumps
  for (m = 1; m < FUNC_CDF; m++)
  {
    y = (float) mm / (FUNC_CDF - njumps);       // Desired value of y ignoring jumps

    while (pdf_array[n] < y && n < pdfsteps)    // Work one's way through pdf_array
    {
      if (j < njump_max && jump[j] <= xmin + (n + 1) * xstep)
      {
        cdf->x[m] = xmin + (n + 1) * xstep;     //Not exactly jump but close
        cdf->y[m] = pdf_array[n];
        j++;                    //increment the jump number
        m++;                    //increment the pdf structure number
      }
      n++;
    }

    /* So at this point pdf_array[n-1] < x and pdf_array[n]>x */
    cdf->x[m] = xmin + (n + 1) * xstep;
    cdf->y[m] = pdf_array[n];
    mm++;                       // increment the number associated with the desired y ignoring jumps
    /* So pdf->y will contain numbers from 0 to 1 */
  }

  cdf->x[FUNC_CDF] = xmax;
  cdf->y[FUNC_CDF] = 1.0;
  cdf->norm = 1.;               /* pdf_gen_from array produces a properly nomalized cdf and so the
                                   normalization is 1.  110629 ksl */
  cdf->ncdf = FUNC_CDF;

/* Calculate the gradients */
  if (calc_cdf_gradient (cdf))
  {
    Error ("cdf_gen_from_func: Errro returned from calc_cdf_gradient\n");
  }                             // 57ib 


  /* Check the pdf */
  if ((icheck = cdf_check (cdf)) != 0)
  {
    Error ("cdf_gen_from_function: error %d on cdf_check\n", icheck);
  }
  return (icheck);

}


/* 

gen_array_from_func


This is a routine which is called by cdf_gen_from_func which simply calculates the cumulative
distribution of the function in equally spaced steps between xmin and xmax.  The CDF is properly
normalized.


pdfsteps is the number of steps that are calculated.  The routine returns the largest change
in the cdf between any two points in the grid.

Notes:

10oct - ksl -This is rather brute force.  An alternative would have been to have increased the density of points
only in the region where it was needed.  This could be done by looking at the regions where the
sampling was poor in the initial calculation and addeding new points in in a new calculations.  But
this would require one to pass both the values of the CDF and the positions where the values were
taken.  This would be  not be hard but would require more extensive modifications to pdf_gen_from func than 
currently.  

11jun	ksl	The name of this function is rather misleading reflecting its history.  At one point
		it simple generated the value of the function at the various points in the grid.  That
		is no longer the case, and has led to confusion.  

History

10oct	ksl(python_69)	Coded to enable one to adaptively increase the density of points
11jun	ksl(69d)	This routine creates a properly normalized CDF.

*/

double
gen_array_from_func (func, xmin, xmax, pdfsteps)
     double (*func) (double);
     double xmin, xmax;
     int pdfsteps;
{

  double x, z, xstep;
  double sum;
  int m, n;
  int idiag;
  double delta;

  xstep = (xmax - xmin) / pdfsteps;
  idiag = 0;

  /* First, allocate an array for internal use, if that
     has not already been done */


  if (init_pdf == 0 || pdf_steps_current < pdfsteps)
  {

    if (pdf_array != NULL)
      free (pdf_array);

    if ((pdf_array = calloc (sizeof (x), pdfsteps)) == NULL)
    {
      Error ("pdf: Could not allocate space for pdf_array\n");
      exit (0);
    }
    init_pdf = 1;
    pdf_steps_current = pdfsteps;
  }


  for (n = 0; n < pdfsteps; n++)
  {
    x = xmin + (n + 0.5) * xstep;
    if ((z = (*func) (x)) < 0 || z > VERY_BIG || sane_check (z))
    {
      Error ("pdf_gen_from_func: probability density %g < 0 at %g\n", z, x);
    }
    if (n == 0)
      pdf_array[0] = z;
    else
      pdf_array[n] = pdf_array[n - 1] + z;
    /* Check to see if the integral seems too large */
    if (pdf_array[n] > 1.0e100 && idiag == 0)
    {
      idiag = 1;
      Log ("pdf_gen_from_func: ZOWIE  n %d z %g pdf_array[n] %g x %g\n", n, z, pdf_array[n], x);
      for (m = 0; m < n; m += 100)
      {
        z = (*func) (x);
        Log ("pdf_gen_from_func: zowie n %d x %g z %g pdf_array %g\n", m, x = xmin + (0.5 * m) * xstep, z, pdf_array[m]);
      }
    }
  }
  /* Thus, pdf_array is proportional to  the definite integral from xmin to x 
     (where x=xmin+(n+1)*xstep) */

  sum = pdf_array[pdfsteps - 1];

  if (sane_check (sum))
  {
    Error ("pdf_gen_from_func:sane_check Sum %f is NaN\n", sum);
  }

  /* Renormalize the array so that this really is a cumulative distribution function */
  delta = 0;
  for (n = 0; n < pdfsteps; n++)
  {
    pdf_array[n] /= sum;
    if (n > 0)
    {
      x = pdf_array[n] - pdf_array[n - 1];
      if (x > delta)
      {
        delta = x;
      }
    }

  }


  return (delta);
}

// End of function


/*  
cdf_gen_from_array

The one dimensional arrays containg the places x where the probability
density y has been evaluated.  The array need not be uniformly spaced.  
The assumption made is that y is the probability at the point x and that 
one can linearly interpolate between points 

Notes:


	06sep -- There are differences between this and cdf_gen_from_func
	that look historical -- ksl


History:
	98jul	ck	Added check to see if pdf, really cdf, is all zeros
	98jul27 ksl	Changed the logic in the way the pdf array is created.  
			The old method was not correct.  
	02jul	ksl	Rewrote much of routine to try to handle situations in
			which the array implied a distribution which was very skewed 
			toward either end. Specifically, I eliminated the old approach
			that did a pseudo integration along the lines of the one
			in pdf_gen_from_func. The new approach should also speed the
			routine up significantly, and eliminate some of the other
			errors which have plagued earlier versions, e.g. non-monotonic
			pdfs and values just out of bounds.  (See the comments in
			python_43/pdf.c)
	04mar	ksl	A familiar problem of the pdf being out of order appeared 
			in pdf_gen_from_array when I increased the number of 
			elements in the pdf array.  I found an error in the way q 
			for jumps was calculated and added checks (and a shuffle)
			for situations where the cdf was out of order.
	090122	ksl	Increase size of PDF_ARRAY to reflect the fact that we are
			tending towards models with more points in them and updated 
			the error message associated with allowing too few points
			so that problems with this would be easier to update in
			future.
	1405	JM -- Increased PDF array for use with disk14 models
	1707	NSH -- Removed code for jumps - we now just supply a suitable unscaled pdf, 
			The code expects data between xmin and xmax - it scan cope with 
			leading ro training zeros, but does require at least two non zero
			points, otherwise it cannt make a CDF.
*/

#define PDF_ARRAY  28000

double pdf_x[PDF_ARRAY], pdf_y[PDF_ARRAY], pdf_z[PDF_ARRAY];
int pdf_n;

int
cdf_gen_from_array (cdf, x, y, n_xy, xmin, xmax)
     CdfPtr cdf;
     double x[], y[];
     int n_xy;
     double xmin, xmax;
{
  int allzero;
  int nmin, nmax, cdf_n;
  int n;
  double sum;
  int echeck;

/* Perform various checks on the inputs */

  if (n_xy > NCDF)
  {
    Error ("cdf_gen_from_array: supplied data %i is larger than default aray size %i - increase NCDF\n", n_xy, NCDF);
    exit (0);
  }


  if (xmax < xmin)
  {
    Error ("cdf_gen_from_array: xmin %g <= xmax %g\n", xmin, xmax);
    exit (0);
  }

  allzero = 0;

  for (n = 1; n < n_xy; n++)
  {
    if (x[n] <= x[n - 1])
    {
      Error ("cdf_gen_from_array: input x not in ascending order at element %5d/%5d  %11.6e %11.6e\n", n, n_xy, x[n - 1], x[n]);
      exit (0);
    }
    if (y[n] < 0)
    {
      Error ("cdf_gen_from_array: probability density %g < 0 at element %d\n", y[n], n);
      exit (0);
    }
    else if (y[n] > 0)
    {
      allzero = 1;
    }
  }







  /* ksl - exiting if xmin and xmax do not conform to the limits of the array should not be a reqquirement
   */

//OLD  if (x[0] != xmin || x[n_xy - 1] != xmax)
//OLD    {
//OLD      Error
//OLD   ("cdf_gen_from_array: input array does not run exactly from xmin to xmax %e != %e or %e != %e\n",
//OLD    x[0], xmin, x[n_xy - 1], xmax);
//OLD      for (n = 0; n < n_xy; n++)
//OLD   printf ("%i %i %e %e %e %e\n", n_xy, n, xmin, xmax, x[n], y[n]);
//OLD      exit (0);

//OLD    }

  /* ksl - Our requrements are as follows.  We want the cdf to have boundaries given by xmin and xmax, unless
   * the array does not encompass this, in which case we will only include what the array contains.  We also
   * want to suppress regions at the end of the cdf where the pdf was 0
   */


  nmin = 0;
  while (x[nmin] < xmin && nmin < n_xy)
  {
    nmin++;
  }

  if (nmin == n_xy)
  {
    Error ("cdf_gen_from_array: xmin greater tahn all array values\n");
    exit (0);
  }

  // if xmin is equal to ove of the values in the x array then nmin will point to that value, but otherwise it 
  // will be slightly larger and we will need to fix this
  // if nmin=0 then it is also possible that the xmin was below the array
  nmax = nmin;
  while (x[nmax] < xmax && nmax < n_xy)
  {
    nmax++;
  }

  /* In general, nmax should be one past the last element */

  if (nmax == nmin)
  {
    Error ("cdf_gen_from_array: nmin and mnax are identical which is not desirable\n");
    exit (0);
  }


  /* Now we want to deal with situations where xmin and xmax are between elements. For now
   * we do not worry about interpolating*/

  if (nmin > 0 && x[nmin] > xmin)
  {
    nmin--;
    x[nmin] = xmin;
  }

  if (x[nmax - 1] > xmax)
  {
    x[nmax - 1] = xmax;
  }

  /* So at this point we want the array running from nmin to nmax-1 */





  /* The next two checks look to see if there is a part of the CDF that is all zeros as the start or end of the distribution
   *
   * ksl - It was not obvious why we car if part of the array is 0, and it could just as well happen in the middle of the array as at the ends 
   * for now I am commented all of this out

   Start first */

//OLD  n = 0;
//OLD  nmin = -1;
//OLD  while (y[n] == 0.0)
//OLD    {
//OLD      nmin = n;
//OLD      n++;
//OLD    }

//OLD  nmin++;



  //Now at the end

//OLD  n = n_xy - 1;
//OLD  nmax = n_xy;
//OLD  while (y[n] == 0.0)
//OLD    {
//OLD      nmax = n;
//OLD      n--;
//OLD    }

//OLD  nmax--;

// ksl I don;t believe the first two cases can happen at this point
  if (nmax == nmin)
  {
    Error ("cdf_gen_from_array - only one point in supplied PDF\n");
    exit (0);
  }

  if (xmax < x[0] || xmin > x[n_xy - 1] || allzero == 0)
  {                             // These are special (probably nonsensical) cases
    cdf->x[0] = xmin;
    cdf->y[0] = 0.;
    cdf->x[1] = xmax;
    cdf->y[1] = 1.;
    cdf->norm = 1.0;
    cdf->ncdf = 1;
    Error ("cdf_gen_from_array: all y's were zero or xmin xmax out of range of array x-- returning uniform distribution %d\n", allzero);

  }
  else
  {
/* So at this point, have unscaled probability density in pdf_x, pdf_y for the points
 * specified by the input array but we want the cumulative distribution
 * We also have assured that there is one value of pdf_x that corresponds to all of the jumps
 */

    // The following lines perform an integration via the trapezoid rule - each point contains
    // the integral up to that poont, so it starts at 0 and ends at the total

    cdf_n = (nmax - nmin);

    cdf->x[0] = x[nmin];
    cdf->y[0] = 0.0;
    for (n = 1; n < cdf_n + 1; n++)
    {
      cdf->x[n] = x[nmin + n];
      cdf->y[n] = cdf->y[n - 1] + 0.5 * (y[nmin + n - 1] + y[nmin + n]) * (x[nmin + n] - x[nmin + n - 1]);
    }

    cdf->y[cdf_n] = cdf->y[n - 1] + y[nmin + cdf_n] * (x[nmax] - x[nmax - 1]);

    // sum = cdf->y[cdf_n - 1];       //the total integrated pdf
    sum = cdf->y[cdf_n];        //the total integrated pdf

    // for (n = 1; n < cdf_n; n++)
    for (n = 1; n <= cdf_n; n++)
    {
      cdf->y[n] /= sum;         //this is now a cdf - we go from 0 to 1.
    }

    cdf->ncdf = cdf_n;
    cdf->x[cdf->ncdf] = x[nmax];
    cdf->y[cdf->ncdf] = 1.0;
    cdf->norm = sum;

  }


/* Calculate the gradients */
  if (calc_cdf_gradient (cdf))
  {
    Error ("cdf_gen_from_array: Error returned from calc_cdf_gradient\n");
  }                             // 57ib 
  if ((echeck = cdf_check (cdf)) != 0)
  {
    Error ("cdf_gen_from_array: error %d on cdf_check\n", echeck);
    for (n = 0; n < n_xy; n++)
      printf ("pdf_n=%i %e %e\n", pdf_n, x[n], y[n]);
    for (n = 0; n < cdf->ncdf + 1; n++)
      printf ("cdf_n=%i %e %e\n", n, cdf->x[n], cdf->y[n]);


    exit (0);
  }

  cdf_to_file (cdf, "foo.diag");
  return (echeck);

}





/* 
cdf_get_rand


Our searching routine assumes we can predict roughly where in the
structure the right values will be 

History:

	06sep	ksl	57i -- Modified to account for the fact
			that the probability density is not 
			uniform between intervals.
    17jun   ksl Modified for variable sizes of distribution function	
	17jul	nsh Modified to use gsl routine to find the right point in the CDF.
*/

double
cdf_get_rand (cdf)
     CdfPtr cdf;
{
  double x, r;
  int i, j;
  double q;
  double a, b, c, s[2];
  int xquadratic ();
/* Find the interval within which x lies */
  r = rand () / MAXRAND;        /* r must be slightly less than 1 */
//  i = r * cdf->ncdf;          /* so i initially lies between 0 and the size of the pdf array -1 */
//  while (cdf->y[i + 1] < r && i < cdf->ncdf - 1)
//    i++;
//  while (cdf->y[i] > r && i > 0)
//    i--;

  i = gsl_interp_bsearch (cdf->y, r, 0, cdf->ncdf);


/* Now calculate a place within that interval */
  q = rand () / MAXRAND;
  a = 0.5 * (cdf->d[i + 1] - cdf->d[i]);
  b = cdf->d[i];
  c = (-0.5) * (cdf->d[i + 1] + cdf->d[i]) * q;
  if ((j = xquadratic (a, b, c, s)) < 0)
  {
    Error ("cdf_get_rand: %d\n", j);
  }
  else
  {
    q = s[j];
    if (q < 0 || q > 1)
    {
      Error ("cdf_get_rand: q out of range %d  %f\n", j, q);
    }
  }

  x = cdf->x[i] * (1. - q) + cdf->x[i + 1] * q;
  if (!(cdf->x[0] <= x && x <= cdf->x[cdf->ncdf]))
  {
    Error ("cdf_get_rand: %g %d %g %g\n", r, i, q, x);
  }
  return (x);
}



/* 
We have defined the cdf in such a way that one can generate a random number 0 and
1, find the place in the y array that corresponds to this and map back onto the 
x values (which is what is returned. We now want to limit the range of the returned
x values, but we need to know in terms of the original array where thiese values
are and the generate a random number between some minimum (>0) and maximum (<1) whichi
will allow one to find the x that corresponds. The point of this routine is to
find that minimum and maximum value. Having done this one can call pdf_get_rand_limit


Error Some of this is not quite right.  It's intended to let us generate limits in
the CDF but we need x too

	06sep	ksl	57i -- Added the xlimits to the cdf
    17jun ksl Modified for variable sizes of distribution function
*/

int
cdf_limit (cdf, xmin, xmax)
     CdfPtr cdf;
     double xmin, xmax;
{
  int i;
  double q;
  if (cdf->y[cdf->ncdf] != 1.0)
  {
    Error ("cdf_limit: cdf not defined!)");
    exit (0);
  }
  if (xmin >= cdf->x[cdf->ncdf])
  {
    Error ("cdf_limit: xmin %g > cdf->x[cdf->ncdf] %g\n", xmin, cdf->x[cdf->ncdf]);
//      exit (0);
  }
  if (xmax <= cdf->x[0])
  {
    Error ("cdf_limit: xmax %g < cdf->x[0] %g\n", xmax, cdf->x[0]);
    exit (0);
  }

/* Set the limits for the minimum */

  if (xmin <= cdf->x[0])
  {
    cdf->limit1 = 0;
    cdf->x1 = cdf->x[0];
  }
  else
  {
    cdf->x1 = xmin;
    i = 0;
    while (xmin > cdf->x[i])
    {
      i++;
    }
    q = (xmin - cdf->x[i - 1]) / (cdf->x[i] - cdf->x[i - 1]);
    cdf->limit1 = cdf->y[i - 1] + q * (cdf->y[i] - cdf->y[i - 1]);
  }

/* Now set the limits for the maximum */

  if (xmax >= cdf->x[cdf->ncdf])
  {
    cdf->limit2 = 1.0;
    cdf->x2 = cdf->x[cdf->ncdf];
  }
  else
  {
    cdf->x2 = xmax;
    i = cdf->ncdf;
    while (xmax <= cdf->x[i])
    {
      i--;
    }
    q = (xmax - cdf->x[i]) / (cdf->x[i + 1] - cdf->x[i]);
    cdf->limit2 = cdf->y[i] + q * (cdf->y[i + 1] - cdf->y[i]);
  }

  return (0);
}

/*
cdf_get_rand_limit (cdf)

History
	06sep	ksl	57h -- Modified to account for the fact
			that the probability density is not 
			uniform between intervals.
    17jul   ksl Modified for variable lengths of pdfs

*/
double
cdf_get_rand_limit (cdf)
     CdfPtr cdf;
{
  double x, r;
  int i, j;
  double q;
  double a, b, c, s[2];
  int xquadratic ();
  r = rand () / MAXRAND;        /* r must be slightly less than 1 */
  r = r * cdf->limit2 + (1. - r) * cdf->limit1;
  i = r * cdf->ncdf;
  while (cdf->y[i + 1] < r && i < cdf->ncdf - 1)
    i++;
  while (cdf->y[i] > r && i > 0)
    i--;
  while (TRUE)
  {
    q = rand () / MAXRAND;
    a = 0.5 * (cdf->d[i + 1] - cdf->d[i]);
    b = cdf->d[i];
    c = (-0.5) * (cdf->d[i + 1] + cdf->d[i]) * q;
    if ((j = xquadratic (a, b, c, s)) < 0)
    {
      Error ("pdf_get_rand: %d\n", j);
    }
    else
    {
      q = s[j];
      if (q < 0 || q > 1)
      {
        Error ("cdf_get_rand_limit: q out of range %d  %f\n", j, q);
      }
    }

    x = cdf->x[i] * (1. - q) + cdf->x[i + 1] * q;
    if (cdf->x1 < x && x < cdf->x2)
      break;
  }

  return (x);
}

/* 

Write the full structure of the cumulattive distribution function to
a file

110628	ksl	Added lines to write out all of the information contained
		in the pdf structure in an attempt to debug an error in
		anisowind
*/

int
cdf_to_file (cdf, filename)
     CdfPtr cdf;
     char filename[];
{
  FILE *fopen (), *fptr;
  int n;
  fptr = fopen (filename, "w");
  fprintf (fptr, "# limits (portion.to.sample)   %10.4g %10.4g\n", cdf->limit1, cdf->limit2);
  fprintf (fptr, "# x1 x2  Range(to.be.returned) %10.4g %10.4g\n", cdf->x1, cdf->x2);
  fprintf (fptr, "# norm   Scale.factor          %10.4g \n", cdf->norm);
  fprintf (fptr, "#n x y  1-y d\n");
  for (n = 0; n <= cdf->ncdf; n++)
    fprintf (fptr, "%3d %14.8e	%14.8e %14.8e  %14.8e\n", n, cdf->x[n], cdf->y[n], 1. - cdf->y[n], cdf->d[n]);
  fclose (fptr);
  return (0);
}

/* This is a simple routine to check some basics of the
   cumulative distribution function.

   History:
   98oct        ksl     Error returns improved so they can be interpreted.
   Also the requirement for a monotonic distribution
   function has been relaxed, in the sense that two
   successive values in x and y can be the same.  This was
   allowed because there was a case in which I had
   required values of the pdf at specific points but
   in the actual event, the partial distribution function
   between these two points was everywhere zero and hence
   the cumulative distribution function had two y values
   the same.  This can cause an error in pdf_get_rand
   in principle.
 */

int
cdf_check (cdf)
     CdfPtr cdf;
{
  int n;
  double x, y;
  int hcheck, icheck, jcheck, kcheck, fcheck;
  hcheck = icheck = jcheck = kcheck = fcheck = 0;
  x = cdf->x[0];
  y = cdf->y[0];
  if (y != 0.0)
  {
    Error ("cdf_check: cumulative distribution function should start at 0 not %e\n", y);
    hcheck = 1;
    printf ("%e %e %e %e\n", cdf->x[0], cdf->y[0], cdf->x[1], cdf->y[1]);
  }
  if (cdf->y[cdf->ncdf] != 1.0)
  {
    Error ("cdf_check: cumulative distribution function should end at 1 not %e\n", cdf->y[cdf->ncdf - 1]);
    icheck = 1;
  }

  for (n = 1; n < cdf->ncdf + 1; n++)
  {
    // Note the equal sign here 
    if (x <= cdf->x[n])
      x = cdf->x[n];
    else
    {
      jcheck = 1;
      Error ("cdf_check: x problem n %d x %f pdf->x[n] %f\n", n, x, cdf->x[n]);
    }
    // Note the equal sign here 
    if (y <= cdf->y[n])
      y = cdf->y[n];
    else
    {
      kcheck = 1;
      Error ("cdf_check: y problem n %d y %f pdf->y[n] %f\n", n, y, cdf->y[n]);
    }
  }
  if (jcheck == 1)
  {
    Error ("cdf_check: x values in pdf should be monotonic\n");
  }
  if (kcheck == 1)
  {
    Error ("cdf_check: cumulative distribution function should be monotonic\n");
  }

  if (hcheck != 0 || icheck != 0 || jcheck != 0 || kcheck != 0)
  {
    cdf_to_file (cdf, "cdf.diag");
    fcheck = hcheck * 1000 + icheck * 100 + jcheck * 10 + kcheck;
    Error ("cdf_check %d %d %d %d %d\n", fcheck, hcheck, icheck, jcheck, kcheck);
  }

/* Next section is a dangerous attempt to repair the cdf  */

  if (hcheck != 0)
    cdf->y[0] = 0.0;
  if (icheck != 0)
    cdf->y[cdf->ncdf] = 1.0;
  if (jcheck != 0)
  {
    for (n = 0; n < cdf->ncdf; n++)
    {
      if (cdf->x[n] >= cdf->x[n + 1])
        cdf->x[n + 1] = cdf->x[n] + 1.e-20;
    }
  }
  if (kcheck != 0)
  {
    for (n = 0; n < cdf->ncdf; n++)
    {
      if (cdf->y[n] >= cdf->y[n] + 1)
        cdf->y[n + 1] = cdf->y[n] + 1.e-20;
    }
  }

  return (fcheck);
}


/***********************************************************
                Space Telescope Science Institute
                                                                                             
Synopsis:
	calc_cdf_gradient(cdf)
                                                                                             
Arguments:
                                                                                             
Returns:
                                                                                             
Description:

	Calculate the gradient of a cdf at each point
	allow one to calculate a first order correction
	to a uniform distibution between the points
	where the gradient has been calculated.
                                            
	XXX NSH - the way that the ends are dealt with is bot
		great - we should really try to come up with an
		extrapolation rather than just fill in the same gradients
		for the second and penultimate cells into the 
		forst and last cells.
											 											
Notes:

	d is calculated by differencing the cdf.  Thus
	this routine should work both for cdf's calculated
	from functions and from arrays
                                                                                             
History:
	06sep	ksl	57i -- Initiated to try an improve
			the calculation of the distribuion
			within a pdf interval.
	06nov	ksl	58b: Fixed problem occuring when there
			were two points in cdf with same x
	17jul	nsh modified name to better relfect what it does
                                                                                             
**************************************************************/

int
calc_cdf_gradient (cdf)
     CdfPtr cdf;
{
  int n, istat;
  double dx1, dx2, dy1, dy2;

  istat = 0;
  for (n = 1; n < cdf->ncdf; n++)
  {
    dy1 = cdf->y[n] - cdf->y[n - 1];
    dx1 = cdf->x[n] - cdf->x[n - 1];
    dy2 = cdf->y[n + 1] - cdf->y[n];
    dx2 = cdf->x[n + 1] - cdf->x[n];
    if (dx1 != 0.0 && dx2 != 0.0)
    {
      cdf->d[n] = 0.5 * (dy1 / dx1 + dy2 / dx2);
    }
    else if (dx1 != 0.0)
    {
      cdf->d[n] = dy1 / dx1;
    }
    else if (dx2 != 0.0)
    {
      cdf->d[n] = dy2 / dx2;
    }
    else
    {
      cdf->d[n] = 0.0;
      Error ("calc_cdf_gradient: dx1 and dx2 both 0 at  %3d %11.6e\n", n, cdf->x[n]);
      istat = 1;
    }

  }
  /* Fill in the ends */
  cdf->d[0] = cdf->d[1];


  //XXX NSH - this routine does not treat the ends of the array very well - it just asumes the same gradient - we could do better...

//NSH - a couple of other ways of doing this - linear interpolation from the previous points, or just the gradient of the current interval - none works great.    

//        cdf->d[0]=cdf->d[1]-((cdf->d[2]-cdf->d[1])/(cdf->x[2]-cdf->x[1]))*(cdf->x[1]-cdf->x[0]);
//      cdf->d[0]=(cdf->y[1] - cdf->y[0])/(cdf->x[1] - cdf->x[0]);
  cdf->d[cdf->ncdf] = cdf->d[cdf->ncdf - 1];
  return (istat);
}


int
cdf_array_fixup (x, y, n_xy)
     double *x, *y;
     int n_xy;
{
  int n, m;
  size_t *order;
  double *xx, *yy;

  order = calloc (sizeof (size_t), n_xy);
  xx = calloc (sizeof (double), n_xy);
  yy = calloc (sizeof (double), n_xy);

  gsl_sort_index (order, x, 1, n_xy);

  for (n = 0; n < n_xy; n++)
  {
    xx[n] = x[order[n]];
    yy[n] = y[order[n]];
  }

  /* So now I know the order, but there could be two x values that are the same */

  m = 0;
  x[0] = xx[0];
  y[0] = yy[0];
  for (n = 0; n < n_xy; n++)
  {
    if (m == 0 && xx[n] > 0)
    {
      x[0] = x[n];
      y[0] = y[n];
      m++;
    }
    else if (xx[n] > xx[n - 1])
    {
      x[m] = xx[n];
      y[m] = yy[n];
      m++;
    }
  }
  free (order);
  free (xx);
  free (yy);



  return (m);
}
