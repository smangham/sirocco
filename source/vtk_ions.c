/***********************************************************/
/** @file   vtk_ions.c
 * @author SWM
 * @date   July, 2015
 * @brief  Mini-functions for writing ions to disk
 *
 * File containing functions for outputting ions to VTK
 ***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
#include "python.h"

void
vtk_temperature_e (WindPtr w, FILE * fptr)
{
  fprintf (fptr, "%g\n", plasmamain[w->nplasma].t_e);
  return;
}

/**********************************************************/
/**
 * @brief      Read the an arbitray wind model in cylindrical
 *     coordinates
 *
 * @param [in] w
 * @param [in] fptr
 * @return     void
 *
 * @details
 * This routine just reads in the data and stores it in arrays
 *
 * ### Notes ###
 * Written by SWM 5/18
 **********************************************************/
void
vtk_ion_c4 (WindPtr w, FILE * fptr)
{
  int element = 6;
  int istate = 4;
  int nion = 0;
  while (nion < nions && !(ion[nion].z == element && ion[nion].istate == istate))
    nion++;

  int nelem = 0;
  while (nelem < nelements && ele[nelem].z != element)
    nelem++;

  double density = plasmamain[w->nplasma].density[nion];
  double nh_number_density = rho2nh * plasmamain[w->nplasma].rho;
  density /= (nh_number_density * ele[nelem].abun);

  fprintf (fptr, "%g\n", density);
  return;
}


/**********************************************************/
/**
 * @brief      Read the an arbitray wind model in cylindrical
 *     coordinates
 *
 * @param [in] w
 * @param [in] fptr
 * @return     void
 *
 * @details
 * This routine just reads in the data and stores it in arrays
 *
 * ### Notes ###
 * Written by SWM 5/18
 **********************************************************/
void
vtk_ion_h1 (WindPtr w, FILE * fptr)
{
  int element = 1;
  int istate = 1;
  int nion = 0;
  while (nion < nions && !(ion[nion].z == element && ion[nion].istate == istate))
    nion++;

  int nelem = 0;
  while (nelem < nelements && ele[nelem].z != element)
    nelem++;

  double density = plasmamain[w->nplasma].density[nion];
  double nh_number_density = rho2nh * plasmamain[w->nplasma].rho;
  density /= (nh_number_density * ele[nelem].abun);

  fprintf (fptr, "%g\n", density);
  return;
}
