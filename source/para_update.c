/***********************************************************
                        University of Southampton

Synopsis: communicate_estimators_para

Arguments:		

Returns:
 
Description:	
	
Notes:
  This was originally done in python.c but I've moved here for more readable code.

History:
    JM Coded as part of fix to #132



**************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
#include "python.h"

int
communicate_estimators_para ()
{
  int mpi_i, mpi_j;
  double *maxfreqhelper, *maxfreqhelper2;
  /*NSH 131213 the next line introduces new helper arrays for the max and min frequencies in bands */
  double *maxbandfreqhelper, *maxbandfreqhelper2, *minbandfreqhelper, *minbandfreqhelper2;
  double *redhelper, *redhelper2;
  int *iredhelper, *iredhelper2;
  // int size_of_helpers;
  int plasma_double_helpers, plasma_int_helpers;

  /* The size of the helper array for doubles. We transmit 10 numbers 
     for each cell, plus three arrays, each of length NXBANDS */
  plasma_double_helpers = (14 + 3 * NXBANDS) * NPLASMA;

  /* The size of the helper array for integers. We transmit 6 numbers 
     for each cell, plus one array of length NXBANDS */
  plasma_int_helpers = (6 + NXBANDS) * NPLASMA;



  maxfreqhelper = calloc (sizeof (double), NPLASMA);
  maxfreqhelper2 = calloc (sizeof (double), NPLASMA);
  /* NSH 131213 - allocate memory for the band limited max and min frequencies */
  maxbandfreqhelper = calloc (sizeof (double), NPLASMA * NXBANDS);
  maxbandfreqhelper2 = calloc (sizeof (double), NPLASMA * NXBANDS);
  minbandfreqhelper = calloc (sizeof (double), NPLASMA * NXBANDS);
  minbandfreqhelper2 = calloc (sizeof (double), NPLASMA * NXBANDS);
  redhelper = calloc (sizeof (double), plasma_double_helpers);
  redhelper2 = calloc (sizeof (double), plasma_double_helpers);
  iredhelper = calloc (sizeof (int), plasma_int_helpers);
  iredhelper2 = calloc (sizeof (int), plasma_int_helpers);


  MPI_Barrier (MPI_COMM_WORLD);
  // the following blocks gather all the estimators to the zeroth (Master) thread


  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
    {
      maxfreqhelper[mpi_i] = plasmamain[mpi_i].max_freq;
      redhelper[mpi_i] = plasmamain[mpi_i].j / np_mpi_global;
      redhelper[mpi_i + NPLASMA] = plasmamain[mpi_i].ave_freq / np_mpi_global;
      redhelper[mpi_i + 2 * NPLASMA] = plasmamain[mpi_i].lum / np_mpi_global;
      redhelper[mpi_i + 3 * NPLASMA] = plasmamain[mpi_i].heat_tot / np_mpi_global;
      redhelper[mpi_i + 4 * NPLASMA] = plasmamain[mpi_i].heat_lines / np_mpi_global;
      redhelper[mpi_i + 5 * NPLASMA] = plasmamain[mpi_i].heat_ff / np_mpi_global;
      redhelper[mpi_i + 6 * NPLASMA] = plasmamain[mpi_i].heat_comp / np_mpi_global;
      redhelper[mpi_i + 7 * NPLASMA] = plasmamain[mpi_i].heat_ind_comp / np_mpi_global;
      redhelper[mpi_i + 8 * NPLASMA] = plasmamain[mpi_i].heat_photo / np_mpi_global;
      redhelper[mpi_i + 9 * NPLASMA] = plasmamain[mpi_i].ip / np_mpi_global;
      redhelper[mpi_i + 10 * NPLASMA] = plasmamain[mpi_i].j_direct / np_mpi_global;
      redhelper[mpi_i + 11 * NPLASMA] = plasmamain[mpi_i].j_scatt / np_mpi_global;
      redhelper[mpi_i + 12 * NPLASMA] = plasmamain[mpi_i].ip_direct / np_mpi_global;
      redhelper[mpi_i + 13 * NPLASMA] = plasmamain[mpi_i].ip_scatt / np_mpi_global;
      for (mpi_j = 0; mpi_j < NXBANDS; mpi_j++)
	{
	  redhelper[mpi_i + (14 + mpi_j) * NPLASMA] = plasmamain[mpi_i].xj[mpi_j] / np_mpi_global;
	  redhelper[mpi_i + (14 + NXBANDS + mpi_j) * NPLASMA] = plasmamain[mpi_i].xave_freq[mpi_j] / np_mpi_global;
	  redhelper[mpi_i + (14 + 2 * NXBANDS + mpi_j) * NPLASMA] = plasmamain[mpi_i].xsd_freq[mpi_j] / np_mpi_global;

	  /* 131213 NSH populate the band limited min and max frequency arrays */
	  maxbandfreqhelper[mpi_i * NXBANDS + mpi_j] = plasmamain[mpi_i].fmax[mpi_j];
	  minbandfreqhelper[mpi_i * NXBANDS + mpi_j] = plasmamain[mpi_i].fmin[mpi_j];

	}
    }
  /* 131213 NSH communiate the min and max band frequencies these use MPI_MIN or MPI_MAX */
  MPI_Reduce (minbandfreqhelper, minbandfreqhelper2, NPLASMA * NXBANDS, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce (maxbandfreqhelper, maxbandfreqhelper2, NPLASMA * NXBANDS, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce (maxfreqhelper, maxfreqhelper2, NPLASMA, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce (redhelper, redhelper2, plasma_double_helpers, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (rank_global == 0)
    {

      Log_parallel ("Zeroth thread successfully received the normalised estimators. About to broadcast.\n");
    }

  MPI_Bcast (redhelper2, plasma_double_helpers, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (maxfreqhelper2, NPLASMA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  /* 131213 NSH Send out the global min and max band limited frequencies to all threads */
  MPI_Bcast (minbandfreqhelper2, NPLASMA * NXBANDS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (maxbandfreqhelper2, NPLASMA * NXBANDS, MPI_DOUBLE, 0, MPI_COMM_WORLD);


  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
    {
      plasmamain[mpi_i].max_freq = maxfreqhelper2[mpi_i];
      plasmamain[mpi_i].j = redhelper2[mpi_i];
      plasmamain[mpi_i].ave_freq = redhelper2[mpi_i + NPLASMA];
      plasmamain[mpi_i].lum = redhelper2[mpi_i + 2 * NPLASMA];
      plasmamain[mpi_i].heat_tot = redhelper2[mpi_i + 3 * NPLASMA];
      plasmamain[mpi_i].heat_lines = redhelper2[mpi_i + 4 * NPLASMA];
      plasmamain[mpi_i].heat_ff = redhelper2[mpi_i + 5 * NPLASMA];
      plasmamain[mpi_i].heat_comp = redhelper2[mpi_i + 6 * NPLASMA];
      plasmamain[mpi_i].heat_ind_comp = redhelper2[mpi_i + 7 * NPLASMA];
      plasmamain[mpi_i].heat_photo = redhelper2[mpi_i + 8 * NPLASMA];
      plasmamain[mpi_i].ip = redhelper2[mpi_i + 9 * NPLASMA];
      plasmamain[mpi_i].j_direct = redhelper2[mpi_i + 10 * NPLASMA];
      plasmamain[mpi_i].j_scatt = redhelper2[mpi_i + 11 * NPLASMA];
      plasmamain[mpi_i].ip_direct = redhelper2[mpi_i + 12 * NPLASMA];
      plasmamain[mpi_i].ip_scatt = redhelper2[mpi_i + 13 * NPLASMA];
      for (mpi_j = 0; mpi_j < NXBANDS; mpi_j++)
	{
	  plasmamain[mpi_i].xj[mpi_j] = redhelper2[mpi_i + (14 + mpi_j) * NPLASMA];
	  plasmamain[mpi_i].xave_freq[mpi_j] = redhelper2[mpi_i + (14 + NXBANDS + mpi_j) * NPLASMA];
	  plasmamain[mpi_i].xsd_freq[mpi_j] = redhelper2[mpi_i + (14 + NXBANDS * 2 + mpi_j) * NPLASMA];

	  /* 131213 NSH And unpack the min and max banded frequencies to the plasma array */
	  plasmamain[mpi_i].fmax[mpi_j] = maxbandfreqhelper2[mpi_i * NXBANDS + mpi_j];
	  plasmamain[mpi_i].fmin[mpi_j] = minbandfreqhelper2[mpi_i * NXBANDS + mpi_j];
	}
    }
  Log_parallel ("Thread %d happy after broadcast.\n", rank_global);

  MPI_Barrier (MPI_COMM_WORLD);

  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
    {
      iredhelper[mpi_i] = plasmamain[mpi_i].ntot;
      iredhelper[mpi_i + NPLASMA] = plasmamain[mpi_i].ntot_star;
      iredhelper[mpi_i + 2 * NPLASMA] = plasmamain[mpi_i].ntot_bl;
      iredhelper[mpi_i + 3 * NPLASMA] = plasmamain[mpi_i].ntot_disk;
      iredhelper[mpi_i + 4 * NPLASMA] = plasmamain[mpi_i].ntot_wind;
      iredhelper[mpi_i + 5 * NPLASMA] = plasmamain[mpi_i].ntot_agn;
      for (mpi_j = 0; mpi_j < NXBANDS; mpi_j++)
	{
	  iredhelper[mpi_i + (6 + mpi_j) * NPLASMA] = plasmamain[mpi_i].nxtot[mpi_j];
	}
    }
  MPI_Reduce (iredhelper, iredhelper2, plasma_int_helpers, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (rank_global == 0)
    {
      Log_parallel ("Zeroth thread successfully received the integer sum. About to broadcast.\n");
    }

  MPI_Bcast (iredhelper2, plasma_int_helpers, MPI_INT, 0, MPI_COMM_WORLD);
  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
    {
      plasmamain[mpi_i].ntot = iredhelper[mpi_i];
      plasmamain[mpi_i].ntot_star = iredhelper[mpi_i + NPLASMA];
      plasmamain[mpi_i].ntot_bl = iredhelper[mpi_i + 2 * NPLASMA];
      plasmamain[mpi_i].ntot_disk = iredhelper[mpi_i + 3 * NPLASMA];
      plasmamain[mpi_i].ntot_wind = iredhelper[mpi_i + 4 * NPLASMA];
      plasmamain[mpi_i].ntot_agn = iredhelper[mpi_i + 5 * NPLASMA];
      for (mpi_j = 0; mpi_j < NXBANDS; mpi_j++)
	{
	  plasmamain[mpi_i].nxtot[mpi_j] = iredhelper2[mpi_i + (6 + mpi_j) * NPLASMA];
	}
    }

  free (maxfreqhelper);
  free (maxfreqhelper2);
  free (maxbandfreqhelper);
  free (maxbandfreqhelper2);
  free (minbandfreqhelper);
  free (minbandfreqhelper2);
  free (redhelper);
  free (redhelper2);
  free (iredhelper);
  free (iredhelper2);

  return (0);
}


/***********************************************************
                        University of Southampton

Synopsis: gather_spectra_para

Arguments:		

Returns:
 
Description:	
	
Notes:

History:
    JM Coded as part of fix to #132

**************************************************************/


int
gather_spectra_para (nspec_helper, nspecs)
     int nspec_helper;
     int nspecs;
{
  double *redhelper, *redhelper2;
  int mpi_i, mpi_j;

  redhelper = calloc (sizeof (double), nspec_helper);
  redhelper2 = calloc (sizeof (double), nspec_helper);


  for (mpi_i = 0; mpi_i < NWAVE; mpi_i++)
    {
      for (mpi_j = 0; mpi_j < nspecs; mpi_j++)
	{
	  redhelper[mpi_i * nspecs + mpi_j] = xxspec[mpi_j].f[mpi_i] / np_mpi_global;

	  if (geo.ioniz_or_extract)	// this is True in ionization cycles only, when we also have a log_spec_tot file
	    redhelper[mpi_i * nspecs + mpi_j + (NWAVE * nspecs)] = xxspec[mpi_j].lf[mpi_i] / np_mpi_global;
	}
    }

  MPI_Reduce (redhelper, redhelper2, nspec_helper, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast (redhelper2, nspec_helper, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  for (mpi_i = 0; mpi_i < NWAVE; mpi_i++)
    {
      for (mpi_j = 0; mpi_j < nspecs; mpi_j++)
	{
	  xxspec[mpi_j].f[mpi_i] = redhelper2[mpi_i * nspecs + mpi_j];

	  if (geo.ioniz_or_extract)	// this is True in ionization cycles only, when we also have a log_spec_tot file
	    xxspec[mpi_j].lf[mpi_i] = redhelper2[mpi_i * nspecs + mpi_j + (NWAVE * nspecs)];
	}
    }
  MPI_Barrier (MPI_COMM_WORLD);

  free (redhelper);
  free (redhelper2);

  return (0);
}



/***********************************************************
                        University of Southampton

Synopsis: communicate_matom_estimators_para

Arguments:		

Returns:
 
Description:	
	
Notes:

History:
    JM Coded as part of fix to #132


**************************************************************/


int
communicate_matom_estimators_para ()
{
  int n, mpi_i;
  double *gamma_helper, *alpha_helper;
  double *level_helper, *kpkt_helper, *jbar_helper;
  double *gamma_helper2, *alpha_helper2;
  double *level_helper2, *kpkt_helper2, *jbar_helper2;

  /* allocate helper arrays for the estimators we want to communicate */
  /* the sizes of these arrays should match the allocation in calloc_estimators in gridwind.c */
  /* we need two arrays for each set of variables. Note that we stick all estimators of
     the same size in the same helper array */
  jbar_helper = calloc (sizeof (double), NPLASMA * size_Jbar_est);
  gamma_helper = calloc (sizeof (double), NPLASMA * 4 * size_gamma_est);
  alpha_helper = calloc (sizeof (double), NPLASMA * 2 * size_alpha_est);
  level_helper = calloc (sizeof (double), NPLASMA * nlevels_macro);
  kpkt_helper = calloc (sizeof (double), NPLASMA);
  jbar_helper2 = calloc (sizeof (double), NPLASMA * size_Jbar_est);
  gamma_helper2 = calloc (sizeof (double), NPLASMA * 4 * size_gamma_est);
  alpha_helper2 = calloc (sizeof (double), NPLASMA * 2 * size_alpha_est);
  level_helper2 = calloc (sizeof (double), NPLASMA * nlevels_macro);
  kpkt_helper2 = calloc (sizeof (double), NPLASMA);

  /* set an mpi barrier before we start */
  MPI_Barrier (MPI_COMM_WORLD);

  /* now we loop through each cell and copy the values of our variables 
     into our helper arrays */
  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
    {
      /* one kpkt_abs quantity per cell */
      kpkt_helper[mpi_i] = plasmamain[mpi_i].kpkt_abs / np_mpi_global;

      for (n = 0; n < nlevels_macro; n++)
	{
	  level_helper[mpi_i + (n * NPLASMA)] = macromain[mpi_i].matom_abs[n] / np_mpi_global;
	}

      for (n = 0; n < size_Jbar_est; n++)
	{
	  jbar_helper[mpi_i + (n * NPLASMA)] = macromain[mpi_i].jbar[n] / np_mpi_global;
	}

      for (n = 0; n < size_gamma_est; n++)
	{
	  gamma_helper[mpi_i + (n * NPLASMA)] = macromain[mpi_i].alpha_st[n] / np_mpi_global;
	  gamma_helper[mpi_i + ((n + size_gamma_est) * NPLASMA)] = macromain[mpi_i].alpha_st_e[n] / np_mpi_global;
	  gamma_helper[mpi_i + ((n + 2 * size_gamma_est) * NPLASMA)] = macromain[mpi_i].gamma[n] / np_mpi_global;
	  gamma_helper[mpi_i + ((n + 3 * size_gamma_est) * NPLASMA)] = macromain[mpi_i].gamma_e[n] / np_mpi_global;
	}

      for (n = 0; n < size_alpha_est; n++)
	{
	  alpha_helper[mpi_i + (n * NPLASMA)] = macromain[mpi_i].recomb_sp[n] / np_mpi_global;
	  alpha_helper[mpi_i + ((n + size_alpha_est) * NPLASMA)] = macromain[mpi_i].recomb_sp_e[n] / np_mpi_global;
	}
    }

  /* because in the above loop we have already divided by number of processes, we can now do a sum
     with MPI_Reduce, passing it MPI_SUM as an argument. This will give us the mean across threads */
  MPI_Reduce (kpkt_helper, kpkt_helper2, NPLASMA, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (level_helper, level_helper2, NPLASMA * nlevels_macro, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (jbar_helper, jbar_helper2, NPLASMA * size_Jbar_est, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (gamma_helper, gamma_helper2, NPLASMA * 4 * size_gamma_est, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (alpha_helper, alpha_helper2, NPLASMA * 2 * size_alpha_est, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank_global == 0)
    {
      Log_parallel ("Zeroth thread successfully received the macro-atom estimators. About to broadcast.\n");
    }

  MPI_Bcast (kpkt_helper2, NPLASMA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (level_helper2, NPLASMA * nlevels_macro, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (jbar_helper2, NPLASMA * size_Jbar_est, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (gamma_helper2, NPLASMA * 4 * size_gamma_est, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (alpha_helper2, NPLASMA * 2 * size_alpha_est, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  /* We now need to copy these reduced variables to the plasma structure in each thread */
  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
    {

      plasmamain[mpi_i].kpkt_abs = kpkt_helper[mpi_i];

      for (n = 0; n < nlevels_macro; n++)
	{
	  macromain[mpi_i].matom_abs[n] = level_helper[mpi_i + (n * NPLASMA)];
	}

      for (n = 0; n < size_Jbar_est; n++)
	{
	  macromain[mpi_i].jbar[n] = jbar_helper[mpi_i + (n * NPLASMA)];
	}

      for (n = 0; n < size_gamma_est; n++)
	{
	  macromain[mpi_i].alpha_st[n] = gamma_helper[mpi_i + (n * NPLASMA)];
	  macromain[mpi_i].alpha_st_e[n] = gamma_helper[mpi_i + ((n + size_gamma_est) * NPLASMA)] / np_mpi_global;
	  macromain[mpi_i].gamma[n] = gamma_helper[mpi_i + ((n + 2 * size_gamma_est) * NPLASMA)];
	  macromain[mpi_i].gamma_e[n] = gamma_helper[mpi_i + ((n + 3 * size_gamma_est) * NPLASMA)];
	}

      for (n = 0; n < size_alpha_est; n++)
	{
	  macromain[mpi_i].recomb_sp[n] = alpha_helper[mpi_i + (n * NPLASMA)];
	  macromain[mpi_i].recomb_sp_e[n] = alpha_helper[mpi_i + ((n + size_alpha_est) * NPLASMA)];
	}
    }

  Log_parallel ("Thread %d happy after broadcast.\n", rank_global);

  MPI_Barrier (MPI_COMM_WORLD);

  free (kpkt_helper);
  free (level_helper);
  free (jbar_helper);
  free (gamma_helper);
  free (alpha_helper);

  free (kpkt_helper2);
  free (level_helper2);
  free (jbar_helper2);
  free (gamma_helper2);
  free (alpha_helper2);

  return (0);
}
