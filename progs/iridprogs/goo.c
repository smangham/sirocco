

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "my_linalg.h"

/*****************************************************************************/

/*
alpha_sp_e - to govern the calculation of the spontaneous recombination rate (energy
weighted).

Very similar to alpha_sp

Modified May 05 2004 by SS to work for continua of simple elements as well as macro 
elements.
*/
#define ALPHA_SP_CONSTANT 5.79618e-36

double
alpha_sp_e (cont_ptr, one)
     struct topbase_phot *cont_ptr;
     WindPtr one;
{
  double alpha_sp_e_value;
  double fthresh, flast;
  double qromb ();
  double alpha_sp_e_integrand ();

  temp_ext = one->t_e;		//external for use in alph_sp_integrand
  cont_ext_ptr = cont_ptr;	//"
  fthresh = cont_ptr->freq[0];	//first frequency in list
  flast = cont_ptr->freq[cont_ptr->np - 1];	//last frequency in list
  alpha_sp_e_value = qromb (alpha_sp_e_integrand, fthresh, flast, 1e-4);

  /* The lines above evaluate the integral in alpha_sp. Now we just want to multiply 
     through by the appropriate constant. */
  if (cont_ptr->macro_info == 1 && geo.macro_simple == 0)
    {
      alpha_sp_e_value = alpha_sp_e_value * config[cont_ptr->nlev].g
	/ config[cont_ptr->uplev].g * pow (one->t_e, -1.5);
    }
  else				//case for simple element
    {
      alpha_sp_e_value = alpha_sp_e_value * config[cont_ptr->nlev].g / ion[cont_ptr->nion + 1].g * pow (one->t_e, -1.5);	//g for next ion up used
    }

  alpha_sp_e_value = alpha_sp_e_value * ALPHA_SP_CONSTANT;

  return (alpha_sp_e_value);
}



/******************************************************************************/

/* alpha_sp_e_integrand. This returns the integrand for alpha_sp_e at a chosen
   frequency*/

double
alpha_sp_e_integrand (freq)
     double freq;		//frequency 
{
  double fthresh;
  double x;
  double integrand;
  double tt;

  fthresh = cont_ext_ptr->freq[0];
  tt = temp_ext;

  if (freq < fthresh)
    return (0.0);		// No recombination at frequencies lower than the threshold freq occur

  x = sigma_phot_topbase (cont_ext_ptr, freq);	//this is the cross-section
  integrand =
    x * freq * freq * exp (H_OVER_K * (fthresh - freq) / tt) * freq / fthresh;

  return (integrand);
}
