
/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	py_wind is a program which can be used to display various parameters of a wind
	as calculated by python.  This file contains many of the subroutiens
	used

Arguments:		


Returns:
 
Description:	
	
		
Notes:

History:
 	97jun	ksl	Coding on py_wind began.
	05jul	ksl	56d -- Moved all of the subroutines to py_wind_sub so that it was
			easier to make templates.

**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
#include "python.h"

#define LINELENGTH 132


/* This routine controls what is displayed on the screen.  
There are basically two options, determined by the variable
"determine"

If determine==1, then zoom sets up the variables py_wind_min,
py_wind_max and py_wind_delta so that the the entire wind is
displayed, but it is subsampled in the x direction.

If determine!=1, then a section of the wind is displayed at
full resolution

History:
	07jul	ksl	Made modifications to allow for the possibility that the
			wind has fewer than 10 elemetns in the x direction
*/
int
zoom (direction)
     int direction;
{
  int center;
  if (direction == 1)
    {				/* then unzoom */
      Log ("Showing selected positions throughout wind\n");
      py_wind_min = 0;
      py_wind_max = NDIM;
      py_wind_delta = NDIM / 10;
/* Allow for the possibility that the wind has an xdimension < 10 */
      if (py_wind_delta < 1)
	py_wind_delta = 1;
    }
  else
    {
      Log ("Select part of wind to display\n");
      center = 5;
      rdint ("Center_x", &center);
      py_wind_min = center - 5;
      if (py_wind_min < 0)
	{
	  Log
	    ("zoom: this choice of center is needlessly close to origin, adjusting py_wind_min to 0\n");
	  py_wind_min = 0;
	}
      py_wind_max = py_wind_min + 10;
      if (py_wind_max > NDIM)
	{
	  py_wind_max = NDIM;
	  Log
	    ("zoom: this choice of py_wind_max is lager than NDIM, adusting py_wind_max to NDIM");
	}
      py_wind_delta = 1;
    }
  return (0);
}


/* overview */
int
overview (w, rootname)
     WindPtr w;
     char rootname[];
{
  double lum, wind_luminosity ();
  int n;
  double heating, lines, ff, photo;

  heating = lines = ff = photo = 0.0;

  for (n = 0; n < NPLASMA; n++)
    {
      heating += plasmamain[n].heat_tot;
      lines += plasmamain[n].heat_lines;
      photo += plasmamain[n].heat_photo;
      ff += plasmamain[n].heat_ff;
    }
  lum = wind_luminosity (0., 1.e20);
  Log (" Total emission %8.2e heating %8.2e\n", lum, heating);
  Log ("    ff emission %8.2e heating %8.2e\n", geo.lum_ff, ff);
  Log ("    fb emission %8.2e heating %8.2e\n", geo.lum_fb, photo);
  Log ("  line emission %8.2e heating %8.2e\n", geo.lum_lines, lines);
  return (0);
}

/* Summary of everything at a given position */

int
position_summary (w)
     WindPtr w;
{
  double x[3], v[3];
  struct photon p;
  int n;
  int nplasma;

  x[0] = geo.wind_rmax / 5;
  x[1] = 0.0;
  x[2] = x[0];

a:Log ("Input x=0,y=0,z=0 to return to main routine\n");
  rddoub ("x", &x[0]);
  rddoub ("y", &x[1]);
  rddoub ("z", &x[2]);

  if (length (x) == 0.0)
    return (0);

  n = where_in_grid (x);
  nplasma = wmain[n].nplasma;

  Log ("Position %8.2e  %8.2e %8.2e  Cell %5d\n", x[0], x[1], x[2], n);

  Log ("Vertex position %8.2e  %8.2e %8.2e  rtheta %8.2e %8.2e \n", w[n].x[0],
       w[n].x[1], w[n].x[2], w[n].r, w[n].theta);
  Log ("Center position %8.2e  %8.2e %8.2e  rtheta %8.2e %8.2e \n",
       w[n].xcen[0], w[n].xcen[1], w[n].xcen[2], w[n].rcen, w[n].thetacen);

  Log ("Electron density: %8.2e  rho %8.2e\n", plasmamain[nplasma].ne,
       plasmamain[nplasma].rho);
  Log ("Vel cell: %8.2e  %8.2e %8.2e\n", w[n].v[0], w[n].v[1], w[n].v[2]);

  p.x[0] = x[0];
  p.x[1] = x[1];
  p.x[2] = x[2];

  vwind_xyz (&p, v);
  Log ("Velocity: %8.2e  %8.2e %8.2e\n", v[0], v[1], v[2]);


  goto a;

  return (0);
}


/* A summary of the energy absorbed in a cell */

int
abs_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  double x, xtot;
  char c;
  char name[LINELENGTH];
  char filename[LINELENGTH];
  int nplasma;

  c = 't';

  printf
    ("Absorption tot=t, lines=l,f=ff,b=fb,h=hydrogen,i=he1,j=he2,z=heavy elements\n");
  rdchar ("Choice", &c);
  switch (c)
    {
    case 't':			/* Total  Absorption */
      strcpy (name, "Total Absorbtion");
      break;
    case 'f':			/*ff */
      strcpy (name, "Free Free Absorbtion");
      break;
    case 'b':			/*Photoionization */
      strcpy (name, "Total photoionization Heating");
      break;
    case 'l':			/* Line  heating */
      strcpy (name, "Resonance Line Heating");
      break;
    case 'h':			/* H photoionization */
      strcpy (name, "H  Photoionization Heating");
      break;
    case 'i':			/* He1 photo */
      strcpy (name, "He I Photoionization Heating");
      break;
    case 'j':			/* He2 photo */
      strcpy (name, "He 2 Photoionization Heating");
      break;
    case 'z':			/* Metal photo */
      strcpy (name, "Metal Photoionization Heating");
      break;
    default:
      printf ("Not a valid choice\n");
      return (0);
    }


  strcpy (filename, rootname);

  for (n = 0; n < NDIM2; n++)
    {
      x = 0.0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  switch (c)
	    {
	    case 't':
	      {			/* Total heating */
		x = plasmamain[nplasma].heat_tot;
		strcat (filename, "heat_tot");
		break;
	    case 'f':		/*ff heating */
		x = plasmamain[nplasma].heat_ff;
		strcat (filename, "heat_ff");
		break;
	    case 'b':		/*photoionization heating */
		x = plasmamain[nplasma].heat_photo;
		strcat (filename, "heat_photo");
		break;
	    case 'l':		/* Line heating */
		x = plasmamain[nplasma].heat_lines;
		strcat (filename, "heat_lines");
		break;
	    case 'h':		/* H heating */
		x = plasmamain[nplasma].heat_ion[0];
		strcat (filename, "heat_h");
		break;
	    case 'i':		/* He1 heating */
		x = plasmamain[nplasma].heat_ion[2];
		strcat (filename, "heat_he1");
		break;
	    case 'j':		/* He2 heating */
		x = plasmamain[nplasma].heat_ion[3];
		strcat (filename, "heat_he2");
		break;
	    case 'z':		/* Line heating of high z elements */
		x = plasmamain[nplasma].heat_z;
		strcat (filename, "heat_z");
		break;
	    default:
		printf ("Not a valid choice\n");
		return (0);
	      }
	      xtot += x;

	    }
	  aaa[n] = x;
	}
    }

  display (name);
  Log ("Component heating %8.3e\n", xtot * py_wind_delta);

  if (ochoice)
    {
      write_array (filename, ochoice);
    }
  return (0);

}

/* A summary of adiabatic cooling */
int
adiabatic_cooling_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  double tot;
  double adiabatic_cooling ();
  char filename[LINELENGTH];
  double t_e;



  tot = 0.0;
  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  t_e = plasmamain[w[n].nplasma].t_e;
	  num_recomb (&plasmamain[w[n].nplasma], t_e);
	  tot += aaa[n] = adiabatic_cooling (&w[n], t_e);
	}
    }


  display ("Adiabatic cooling");
  Log ("The total adiabatic cooling is %8.2g\n", tot);

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".adiabatic");
      write_array (filename, ochoice);
    }
  return (0);

}


/* A summary of e lum of a cell */

int
lum_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  double x, xtot;
  char c;
  char name[LINELENGTH];
  char filename[LINELENGTH];
  int nplasma;

  c = 't';
  strcpy (filename, rootname);

  printf
    ("lum tot=t,r=rad, lines=l,f=ff,b=fb,h=hydrogen,i=he1,j=he2,z=heavy elements\n");
  rdchar ("Choice", &c);
  switch (c)
    {
    case 't':			/* Total luminosity */
      strcpy (name, "Total Luminosity");
      strcat (filename, ".tot_lum");
      break;
    case 'r':			/*Radiative energo loss total */
      strcpy (name, "Total Radiative Luminosity");
      strcat (filename, ".rad_lum");
      break;
    case 'f':			/*Radiative energo loss total */
      strcpy (name, "Free Free Luminosity");
      strcat (filename, ".ff_lum");
      break;
    case 'b':			/*Radiative energo loss total */
      strcpy (name, "Free Bound (Total Recombination) Luminosity");
      strcat (filename, ".fb_lum");
      break;
    case 'l':			/* Line lumosity */
      strcpy (name, "Line Luminosity");
      strcat (filename, ".line_lum");
      break;
    case 'h':			/* H lumosity */
      strcpy (name, "H  Recombination Luminosity");
      strcat (filename, ".h1_recomb_lum");
      break;
    case 'i':			/* Line lumosity */
      strcpy (name, "He I Recombination Luminosity");
      strcat (filename, ".he1_recom_lum");
      break;
    case 'j':			/* Line lumosity */
      strcpy (name, "He 2 Recombination Luminosoity");
      strcat (filename, ".he2_recomb_lum");
      break;
    case 'z':			/* Line lumosity */
      strcpy (name, "Metal RecombinationsLuminosity");
      strcat (filename, ".z_recomb_lum");
      break;
    default:
      printf ("Not a valid choice\n");
      return (0);
    }



  xtot = 0.0;


  for (n = 0; n < NDIM2; n++)
    {
      x = 0.0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  switch (c)
	    {
	    case 't':		/* Total luminosity */
	      x = plasmamain[nplasma].lum_rad;
	      break;
	    case 'r':		/*Radiative energo loss total */
	      x = plasmamain[nplasma].lum_rad;
	      break;
	    case 'f':		/*Radiative energo loss total */
	      x = plasmamain[nplasma].lum_ff;
	      break;
	    case 'b':		/*Radiative energo loss total */
	      x = plasmamain[nplasma].lum_fb;
	      break;
	    case 'l':		/* Line luminosity */
	      x = plasmamain[nplasma].lum_lines;
	      break;
	    case 'h':		/* H luminosity */
	      x = plasmamain[nplasma].lum_ion[0];
	      break;
	    case 'i':		/* Line luminosity */
	      x = plasmamain[nplasma].lum_ion[2];
	      break;
	    case 'j':		/* Line luminosity */
	      x = plasmamain[nplasma].lum_ion[3];
	      break;
	    case 'z':		/* Line luminosity */
	      x = plasmamain[nplasma].lum_z;
	      break;
	    default:
	      printf ("Not a valid choice\n");
	      return (0);
	    }
	  xtot += x;

	}
      aaa[n] = x;
    }

  display (name);

  Log ("Luminosity total %8.3e\n", xtot * py_wind_delta);

  if (ochoice)
    {
      write_array (filename, ochoice);
    }
  return (0);

}

/* A summary of the number of photoionizations in  a cell */

int
photo_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n, ion;
  char filename[LINELENGTH];
  int nplasma;

  ion = 0;

  rdint ("ion(<5)", &ion);
  if (!(0 <= ion && ion < NIONIZ))
    {
      Log ("ion out of range\n");
      return (0);
    }


  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] =
	    plasmamain[nplasma].ioniz[ion] * plasmamain[nplasma].density[ion];
	}


    }
  display ("No of ionizations per second in cell");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".photo_sum");
    }
  return (0);

}

/* A summary of the number of recombinations in  a cell */

int
recomb_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int ion;
  char filename[LINELENGTH];
  int nplasma;

  ion = 0;
  rdint ("ion(<5)", &ion);
  if (!(0 <= ion && ion < NIONIZ))
    {
      Log ("recomb_summary: ioniz out of range\n");
      return (0);
    }


  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  num_recomb (&plasmamain[nplasma], plasmamain[nplasma].t_e);
	  aaa[n] = plasmamain[nplasma].recomb[ion];
	}
    }

  display (" Number of recombinations per second in cell");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".recomb");
      write_array (filename, ochoice);
    }

  return (0);

}


/* A summary of electron densities */

int
electron_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;


  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0.0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].ne;
	}
    }
  display ("Electron densities");
  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".ne");
      write_array (filename, ochoice);
    }
  return (0);

}


/* A summary of the average frequency */

int
freq_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int nplasma;

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].ave_freq;
	}
    }
  display ("Average freqency");

  return (0);

}


/* A summary of the number of photons which passed through a cell */

int
nphot_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int nplasma;

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].ntot;
	}
    }
  display ("Nphot per cell");

  return (0);

}

/* A summary of the temperatures */

int
temp_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].t_e;
	}
    }
  display ("T_e");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".te");
      write_array (filename, ochoice);
    }

  return (0);

}


int
temp_rad (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;


  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].t_r;
	}
    }
  display ("T_rad");

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".tr");
      write_array (filename, ochoice);

    }

  return (0);

}

/* A summary of the radiative weights */

int
weight_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int nplasma;

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].w;
	}
    }
  display ("Radiative weights");

  return (0);

}

/* A summary of the velocities */

int
velocity_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  double x;
  char name[LINELENGTH];
  char filename[LINELENGTH];
  int ichoice;

  ichoice = 0;

  /* Note that EOF actually causes an exit via rdpar */
  while (rdint ("|v|=0;v_x=1,v_y=2,v_z=3,return=other", &ichoice) != EOF)
    {
      if (ichoice < 0 || ichoice > 3)
	return (0);

      if (ichoice == 0)
	strcpy (name, "Velocity");
      else if (ichoice == 1)
	strcpy (name, "V_x");
      else if (ichoice == 2)
	strcpy (name, "V_y");
      else
	strcpy (name, "V_z");

      for (n = 0; n < NDIM2; n++)
	{
	  x = 0;
	  if (w[n].vol > 0.0)
	    {
	      if (ichoice == 0)
		x =
		  sqrt (w[n].v[0] * w[n].v[0] + w[n].v[1] * w[n].v[1] +
			w[n].v[2] * w[n].v[2]);
	      else if (ichoice == 1)
		x = w[n].v[0];
	      else if (ichoice == 2)
		x = w[n].v[1];
	      else
		x = w[n].v[2];
	    }
	  aaa[n] = x;
	}
      display (name);


      if (ochoice)
	{
	  strcpy (filename, rootname);

	  if (ichoice == 0)
	    {
	      strcat (filename, ".vel");
	      write_array (filename, ochoice);
	    }
	  if (ichoice == 1)
	    {
	      strcat (filename, ".vrho");
	      write_array (filename, ochoice);
	    }
	  if (ichoice == 2)
	    {
	      strcat (filename, ".vtheta");
	      write_array (filename, ochoice);
	    }
	  if (ichoice == 3)
	    {
	      strcat (filename, ".vz");
	      write_array (filename, ochoice);
	    }

	}
    }

  return (0);

}


/* A summary of the velocities */

int
mo_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int ichoice;
  char name[LINELENGTH];
  double x;
  PlasmaPtr xplasma;


  ichoice = 0;

  /* Note that EOF actually causes an exit via rdpar */
  while (rdint
	 ("|F_rad|=0;Frac_x=1,Frad_y=2,Frad_z=3,return=other",
	  &ichoice) != EOF)
    {
      if (ichoice < 0 || ichoice > 3)
	return (0);

      if (ichoice == 0)
	strcpy (name, "|F_rad|");
      else if (ichoice == 1)
	strcpy (name, "F_rad_x");
      else if (ichoice == 2)
	strcpy (name, "F_rad_y");
      else
	strcpy (name, "F_rad_z");

      for (n = 0; n < NDIM2; n++)
	{
	  xplasma = &plasmamain[w[n].nplasma];
	  x = 0;
	  if (w[n].vol > 0.0)
	    {
	      if (ichoice == 0)
		x =
		  sqrt (xplasma->dmo_dt[0] * xplasma->dmo_dt[0] +
			xplasma->dmo_dt[1] * xplasma->dmo_dt[1] +
			xplasma->dmo_dt[2] * xplasma->dmo_dt[2]);
	      else if (ichoice == 1)
		x = xplasma->dmo_dt[0];
	      else if (ichoice == 2)
		x = xplasma->dmo_dt[1];
	      else
		x = xplasma->dmo_dt[2];
	    }
	  aaa[n] = x;
	}
      display (name);
    }

  return (0);

}

/* A summary of the volumes of each cell */

int
vol_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  aaa[n] = w[n].vol;
	}
    }
  display ("Volumes");

  return (0);

}

int
wind_element (w)
     WindPtr w;
{
  PlasmaPtr xplasma;
  int m, n, i, j, nn;
  int first, last;
  n = 50;
a:rdint ("element", &n);
  if (n < 0)
    goto b;
  wind_n_to_ij (n, &i, &j);
  xplasma = &plasmamain[w[n].nplasma];
  Log ("Element %d (%d,%d)  inwind %d ntot %d nioniz %d nrad %d\n", n, i,
       j, w[n].inwind, xplasma->ntot, xplasma->nioniz, xplasma->nrad);
  Log ("xyz %8.2e %8.2e %8.2e vel %8.2e %8.2e %8.2e\n", w[n].x[0],
       w[n].x[1], w[n].x[2], w[n].v[0], w[n].v[1], w[n].v[2]);
  Log ("nh %8.2e ne %8.2e t_r %8.2e t_e %8.2e w %8.2e vol %8.2e\n",
       xplasma->rho * rho2nh, xplasma->ne, xplasma->t_r, xplasma->t_e,
       xplasma->w, w[n].vol);
  Log
    ("t_e %8.2e lum_tot  %8.2e lum_lines  %8.2e lum_ff  %8.2e lum_fb     %8.2e %8.2e %8.2e %8.2e %8.2e\n",
     xplasma->t_e, xplasma->lum_rad, xplasma->lum_lines, xplasma->lum_ff,
     xplasma->lum_fb, xplasma->lum_ion[0], xplasma->lum_ion[2],
     xplasma->lum_ion[3], xplasma->lum_z);
  Log
    ("t_r %8.2e heat_tot %8.2e heat_lines %8.2e heat_ff %8.2e heat_photo %8.2e %8.2e %8.2e %8.2e %8.2e\n",
     xplasma->t_r, xplasma->heat_tot, xplasma->heat_lines, xplasma->heat_ff,
     xplasma->heat_photo, xplasma->heat_ion[0], xplasma->heat_ion[2],
     xplasma->heat_ion[3], xplasma->heat_z);
  Log ("The ratio of cooling to heating is %8.2f adiabatic cooling %8.2e\n",
       xplasma->lum_rad / xplasma->heat_tot, xplasma->lum_adiabatic);
  Log ("Number of ionizing photons in cell nioniz %d\n", xplasma->nioniz);
  Log ("ioniz %8.2e %8.2e %8.2e %8.2e %8.2e\n",
       xplasma->ioniz[0], xplasma->ioniz[1], xplasma->ioniz[2],
       xplasma->ioniz[3], xplasma->ioniz[4]);
  Log
    ("Convergence status: whole %d converging %d t_r %8.2e t_e %8.2e hc %8.2e \n",
     xplasma->converge_whole, xplasma->converging, xplasma->converge_t_r,
     xplasma->converge_t_e, xplasma->converge_hc);

  for (nn = 0; nn < 4; nn++)
    {
      first = ele[nn].firstion;
      last = first + ele[nn].nions;
      Log ("%-5s ", ele[nn].name);
      for (m = first; m < last; m++)
	Log (" %8.2e", xplasma->density[m]);
      Log ("\n");
    }

  goto a;

b:return (0);

}

int
tau_h_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int nplasma;

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = plasmamain[n].nplasma;
	  aaa[n] =
	    6.e-18 * plasmamain[nplasma].density[0] * pow (w[n].vol, 0.333);

	}
    }
  display ("tau_h summary");

  return (0);

}

/* A summary of the volumes of each cell */

int
coolheat_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  int nplasma;

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = plasmamain[n].nplasma;
	  aaa[n] = plasmamain[nplasma].lum_rad / plasmamain[nplasma].heat_tot;

	}
    }
  display ("Cooling over heating");

  return (0);

}

int
complete_file_summary (w, root, ochoice)
     WindPtr w;
     char root[];
     int ochoice;
{
  temp_summary (w, root, ochoice);
  temp_rad (w, root, ochoice);
  electron_summary (w, root, ochoice);
  ion_summary (w, 6, 3, 1, root, ochoice);
  ion_summary (w, 6, 4, 1, root, ochoice);
  ion_summary (w, 6, 5, 1, root, ochoice);
  ion_summary (w, 7, 4, 1, root, ochoice);
  ion_summary (w, 7, 5, 1, root, ochoice);
  ion_summary (w, 7, 6, 1, root, ochoice);
  ion_summary (w, 8, 4, 1, root, ochoice);
  ion_summary (w, 8, 5, 1, root, ochoice);
  ion_summary (w, 8, 6, 1, root, ochoice);
  ion_summary (w, 14, 3, 1, root, ochoice);
  ion_summary (w, 14, 4, 1, root, ochoice);
  ion_summary (w, 14, 5, 1, root, ochoice);
  return (0);
}

/* A summary of the regions in the wind */

int
wind_reg_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;


  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = -9;
      if (w[n].vol > 0.0)
	{
	  aaa[n] = w[n].inwind;

	}
    }
  display ("Regions of the wind");

  return (0);

}


/* A summary of the dvds_ave */

int
dvds_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  aaa[n] = w[n].dvds_ave;

	}
    }
  display ("Average dvds");

  return (0);
}
