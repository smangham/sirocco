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

void vtk_temperature_e(WindPtr w, FILE* fptr) {
  int n = w->nplasma;
  fprintf (fptr, "%g\n", plasmamain[n].t_e);
  return;
}
