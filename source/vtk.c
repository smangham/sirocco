/***********************************************************/
/** @file   vtk.c
 * @author SWM
 * @date   July, 2015
 * @brief  Output to VTK files
 *
 * File containing functions for outputting to VTK
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
#include "python.h"

/****************************************************************/
/**
 * @brief		Given trz index in wind, returns vtk data index.
 *
 * @param [in] i 	Theta index of cell.
 * @param [in] j 	Radius index of cell.
 * @param [in] k 	Height index of cell.
 * @param [in] i_top Whether the point is above or below the disk.
 * @param [in] resolution Radial bins for visualisation.
 * @return			Index for vtk poly.
 *
 * When given the position of a cell in the wind, returns the
 * corresponding index in the flat VTK data arrays. RZ uses the
 * uses the standard wind location, theta is set as defined in
 * geo. Used in wind_paths_output() only.
 *
 * ###Notes### Written by SWM 4/15.
 *****************************************************************/
int
vtk_point_index (int i, int j, int k, int i_top, DomainPtr dom, int resolution)
{
  int n;
  n = i * 2 * (resolution + 1) * dom->ndim + j * 2 * (resolution + 1) + k * 2 + i_top;
  return (n);
}

/****************************************************************/
/**
 * @brief		Given rtt index in wind, returns vtk data index.
 *
 * @param [in] i 	Theta index of cell.
 * @param [in] j 	Radius index of cell.
 * @param [in] k 	Height index of cell.
 * @param [in] resolution Radial bins for visualisation.
 * @return			Index for vtk poly.
 *
 * When given the position of a cell in the wind, returns the
 * corresponding index in the flat VTK data arrays. RZ uses the
 * uses the standard wind location, theta is set as defined in
 * geo. Used in wind_paths_output() only.
 *
 * ###Notes### Written by SWM 4/15.
 *****************************************************************/
int
vtk_sphere_point_index (int i, int j, int k, int resolution)
{
  int n;
  n = i * (resolution + 1) * (resolution + 1) + j * (resolution + 1) + k;
  return (n);
}

/**********************************************************/
/** @var *double reverb_path_bin
 * @brief	Resolution of cylindrical visualisation
 *
 * Set in reverb_init(), contains the boundaries for each
 * path bin.
 *
 * ###Notes###
 * 5/18	-	Written by SWM
***********************************************************/
int RESOLUTION = 100;

/****************************************************************/
/**
 * @brief		Outputs arbitrary wind information to vkt file
 *
 * @param [in] wind 		Pointer to wind array.
 * @param [in] c_file_in	Name of input file.
 * @param [in] property_pointers Array of pointers to functions.
 * @param [in] property_count Number of entries in that array.
 * @param [in] property_names Names of each property to plot
 * @return 					0
 *
 * When given a wind and an array of functions, calls those
 * functions on every cell of the wind so they can write to a
 * .vtk output file.
 *
 * ###Notes### Written by SWM 5/18.
*****************************************************************/
int
output_vtk (WindPtr wind, int ndom, int cycle,
            int property_count,
            void (*property_pointers[])(WindPtr, FILE*),
            char *property_names[])
{
  FILE *fopen (), *fptr;
  char c_file[LINELENGTH];
  int i, j, k, n, i_cells, i_points, i_property;
  double r_azi, r_inc, r_x, r_y, r_z;
  DomainPtr dom;

  int resolution = RESOLUTION;

  //Get output filename
  sprintf (c_file, "%s.%d.dom%d.wind_paths.vtk", files.root, cycle, ndom);

  if ((fptr = fopen (c_file, "w")) == NULL)
  {                             //If this file can't be opened, error out
    Error ("wind_paths_output_vtk: Unable to open %s for writing\n", c_file);
    exit (0);
  }
  Log ("Outputting wind path information to file '%s'.\n", c_file);

  dom = &zdom[ndom];

  //Setup the number of vertexes and cells within this mesh
  //Notably, we skip the outside cell as it's empty
  //For the spherical case, we generate theta and thi bins using the same resolution
  switch (dom->coord_type)
  {
  case SPHERICAL:
    i_cells = 1 * (dom->ndim - 1) * resolution * resolution;
    i_points = 1 * dom->ndim * (resolution + 1) * (resolution + 1);
    break;
  case CYLIND:
    i_cells = 2 * (dom->ndim - 1) * (dom->mdim - 1) * resolution;
    i_points = 2 * dom->ndim * dom->mdim * (resolution + 1);
    break;
  default:
    fclose (fptr);
    Error ("output_vtk: Mesh format not yet supported (CYLVAR/RTHETA)");
    return (0);
  }

  //Write out header
  fprintf (fptr, "# vtk DataFile Version 2.0\n");
  fprintf (fptr, "Wind file data\nASCII\n");

  //Write out positions of corners of each wind cell as vertexes
  fprintf (fptr, "DATASET UNSTRUCTURED_GRID\n");
  fprintf (fptr, "POINTS %d float\n", i_points);
  if (dom->coord_type == SPHERICAL)
  {
    for (i = 0; i < dom->ndim; i++)
    {
      n = dom->nstart + i;
      for (j = 0; j <= resolution; j++)
      {
        r_inc = j * (PI / (double) resolution) - PI / 2.0;

        for (k = 0; k <= resolution; k++)
        {
          r_azi = k * (PI / (double) resolution);
          r_x = wind[n].x[0] * sin (r_inc) * cos (r_azi);
          r_y = wind[n].x[0] * sin (r_inc) * sin (r_azi);
          r_z = wind[n].x[0] * cos (r_inc);
          fprintf (fptr, "%10.5g %10.5g %10.5g\n", r_x, r_y, r_z);
        }
      }
    }
  }
  else
  {
    for (i = 0; i < dom->ndim; i++)
    {
      for (j = 0; j < dom->mdim; j++)
      {
        wind_ij_to_n (ndom, i, j, &n);
        for (k = 0; k <= resolution; k++)
        {
          r_azi = k * (PI / (double) resolution);
          r_x = wind[n].x[0] * cos (r_azi);
          r_y = wind[n].x[0] * sin (r_azi);
          fprintf (fptr, "%10.5g %10.5g %10.5g\n", r_x, r_y, wind[n].x[2]);
          fprintf (fptr, "%10.5g %10.5g %10.5g\n", r_x, r_y, -wind[n].x[2]);
        }
      }
    }
  }
  fprintf (fptr, "\n");

  //Write out the vertexes comprising each cell
  fprintf (fptr, "CELLS %d %d\n", i_cells, 9 * i_cells);
  if (dom->coord_type == SPHERICAL)
  {
    for (i = 0; i < dom->ndim - 1; i++)
    {
      n = dom->nstart + i;
      for (j = 0; j < resolution; j++)
      {
        r_inc = j * (PI / (double) resolution) - PI / 2.0;

        for (k = 0; k < resolution; k++)
        {
          fprintf (fptr, "8 %d %d %d %d %d %d %d %d\n",
                   vtk_sphere_point_index (i, j, k, resolution),
                   vtk_sphere_point_index (i, j, k + 1, resolution),
                   vtk_sphere_point_index (i, j + 1, k + 1, resolution),
                   vtk_sphere_point_index (i, j + 1, k, resolution),
                   vtk_sphere_point_index (i + 1, j, k, resolution),
                   vtk_sphere_point_index (i + 1, j, k + 1, resolution),
                   vtk_sphere_point_index (i + 1, j + 1, k + 1, resolution),
                   vtk_sphere_point_index (i + 1, j + 1, k, resolution));
        }
      }
    }
  }
  else
  {
    for (i = 0; i < dom->ndim - 1; i++)
    {
      for (j = 0; j < dom->mdim - 1; j++)
      {
        for (k = 0; k < resolution; k++)
        {
          fprintf (fptr, "8 %d %d %d %d %d %d %d %d\n",
                   vtk_point_index (i, j, k, 1, dom, resolution),
                   vtk_point_index (i, j, k + 1, 1, dom, resolution),
                   vtk_point_index (i, j + 1, k + 1, 1, dom, resolution),
                   vtk_point_index (i, j + 1, k, 1, dom, resolution),
                   vtk_point_index (i + 1, j, k, 1, dom, resolution),
                   vtk_point_index (i + 1, j, k + 1, 1, dom, resolution),
                   vtk_point_index (i + 1, j + 1, k + 1, 1, dom, resolution),
                   vtk_point_index (i + 1, j + 1, k, 1, dom, resolution));
          fprintf (fptr, "8 %d %d %d %d %d %d %d %d\n",
                   vtk_point_index (i, j, k, 0, dom, resolution),
                   vtk_point_index (i, j, k + 1, 0, dom, resolution),
                   vtk_point_index (i, j + 1, k + 1, 0, dom, resolution),
                   vtk_point_index (i, j + 1, k, 0, dom, resolution),
                   vtk_point_index (i + 1, j, k, 0, dom, resolution),
                   vtk_point_index (i + 1, j, k + 1, 0, dom, resolution),
                   vtk_point_index (i + 1, j + 1, k + 1, 0, dom, resolution),
                   vtk_point_index (i + 1, j + 1, k, 0, dom, resolution));
        }
      }
    }
  }

  //Write the type for each cell (would be unnecessary if STRUCTURED_GRID used)
  //But this allows for more flexible expansion
  fprintf (fptr, "CELL_TYPES %d\n", i_cells);
  for (i = 0; i < i_cells; i++)
    fprintf (fptr, "12\n");
  fprintf (fptr, "\n");

  //Write out the arrays containing the various properties in the appropriate order
  fprintf (fptr, "CELL_DATA %d\n", i_cells);

  for (i_property = 0; i_property < property_count; i_property++)
  {
    fprintf (fptr, "SCALARS %s float 1\n", property_names[i_property]);
    fprintf (fptr, "LOOKUP_TABLE default\n");
    if (dom->coord_type == SPHERICAL)
    {
      for (i = 0; i < dom->ndim - 1; i++)
      {
        n = dom->nstart + i;
        for (j = 0; j < resolution; j++)
        {
          for (k = 0; k < resolution; k++)
          {
            property_pointers[i_property](&wind[n], fptr);
          }
        }
      }
    }
    else
    {
      for (i = 0; i < dom->ndim - 1; i++)
      {
        for (j = 0; j < dom->mdim - 1; j++)
        {
          wind_ij_to_n (ndom, i, j, &n);
          for (k = 0; k < resolution; k++)
          {
            property_pointers[i_property](&wind[n], fptr);
            property_pointers[i_property](&wind[n], fptr);
          }
        }
      }
    }
  }

  fclose (fptr);
  return (0);
}
