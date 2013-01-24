
/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:

  Description:	
This program reads the file atom1.dat conatining permitted resonace lines compiled
by Verner, Verner and Ferland Atomic Data Nucl Data Tables in press
and reformats it so that it can be used with get_atomic.data.  

The results are written to lines.data.  

The elements selceted are controlled by the structure "element"   


  Arguments:		

  Returns:

  Notes:

  History:
	1/97	ksl 	Source for this modification Written and debugged. 
	01nov	ksl	I have lost the version of this program that
			was used to create lines.data.  Began to modify
			the only version of the program I could find to
			first reproduce lines.data as it exists. 
				(a) fixed element struct to produce same lines for
					same ions as previously
				(b) eliminated lines with unknown wavelength (which
					Verner had set to zero.
	01nov	ksl	Modify output to include uper and lower level energies of
			line.  This will enable get_atomicdata.c to be simplified,
			and a more accurate treatment of the lines in python.
	

 ************************************************************************/


#include <stdio.h>
#include <stdlib.h>

#define LINELENGTH 132
#define WORDLENGTH	25
#define WORDS		11
#define MULTUPLETS  0		/* if 0 get individual lines, other wise sum of multuplet */

#define H   				6.6262e-27
#define HC  				1.98587e-16
#define HEV				4.13620e-15	/* Planck's constant in eV */
#define C   				2.997925e10



int first[] = { 0, 5, 39, 47, 57, 73, 89, 92, 95, 104, 113 };

#define NELEM 14
struct element
{
  char name[5];
  int z;
  int ion_min, ion_max;
}
e[] =
{
"H", 1, 1, 2,
    "He", 2, 1, 2,
    "C", 6, 2, 6,
    "N", 7, 3, 7,
    "O", 8, 3, 7,
    "Ne", 10, 1, 8,
    "Na", 11, 2, 9,
    "Mg", 12, 2, 8,
    "Al", 13, 2, 8,
    "Si", 14, 2, 7,
    "S", 16, 2, 7, "Ar", 18, 1, 9, "Ca", 20, 1, 9, "Fe", 26, 2, 12};



int
main ()
{
  FILE *fptr, *fopen (), *optr;
  char line[LINELENGTH];
  int lengths[WORDS];
  char words[WORDS][WORDLENGTH];
  char element[WORDS], old_element[WORDS];
  int ion, old_ion;
  int nelem;

  int nrecords, nlines;
  int i, j, k, l, m, n, mmax;
  int z, ielec;
  float f, lambda;
  int g_i, g_k;
  float energy_i, energy_k;

  printf ("hello world\n");

  for (n = 0; n < WORDS - 1; n++)
    {
      lengths[n] = first[n + 1] - first[n] - 1;
    }
  lengths[WORDS - 1] = 2;

  fptr = fopen ("atom1.dat", "r");
  optr = fopen ("lines.data", "w");

  nrecords = nlines = 0;
  while (fgets (line, LINELENGTH, fptr) != NULL)
    {
      printf ("%s\n", line);
      nrecords++;
      for (n = 0; n < WORDS; n++)
	{
	  strncpy (&words[n][0], &line[first[n]], lengths[n]);
	}
//              for(n=0;n<WORDS;n++){
//                      printf("%d %s\n",n,words[n]);
//              }

      strcpy (element, "");
      if (sscanf (&words[0][0], "%s %d", element, &ion) != 2)
	{
	  printf ("Skipping record\n");
	  continue;
	}

      nelem = 0;
      while (nelem < NELEM && strncmp (element, e[nelem].name, 2) != 0)
	nelem++;
      if (nelem == NELEM)
	{
//                      printf("Could not match element %s\n",words[0]);
	  continue;
	}

      if (e[nelem].ion_min <= ion && ion <= e[nelem].ion_max)
	{
	  if (MULTUPLETS != 0)
	    {
	      sscanf (&words[3][0], "%f", &lambda);
	      sscanf (&words[9][0], "%f", &f);
	      printf ("Line %d %d %f %f\n", e[nelem].z, ion, lambda, f);
	      fprintf (optr, "Line %d %d %f %f\n", e[nelem].z, ion, lambda,
		       f);
	      nlines++;
	    }
	  else
	    {
	      sscanf (&words[10][0], "%d", &mmax);
	      printf ("Found %s %d %d \n", element, ion, mmax);
	      for (m = 0; m < mmax; m++)
		{
		  if (fgets (line, LINELENGTH, fptr) == NULL)
		    exit (0);
		  for (n = 0; n < WORDS; n++)
		    {
		      strncpy (&words[n][0], &line[first[n]], lengths[n]);
		    }

		  sscanf (&words[3][0], "%f", &lambda);
		  sscanf (&words[4][0], "%f", &energy_i);
		  sscanf (&words[5][0], "%f", &energy_k);
		  energy_i *= (HEV * C);
		  energy_k *= (HEV * C);
		  sscanf (&words[6][0], "%d", &g_i);
		  sscanf (&words[7][0], "%d", &g_k);
		  sscanf (&words[9][0], "%f", &f);
//                                      for(n=0;n<WORDS;n++){
//                                              printf("%d %s\n",n,words[n]);
//                                      }
		  if (lambda > 0)
		    {
		      printf (
// "Line %d %d %f %f %d %d %f %f\n",
			       "Line %2d %2d %11.6f %11.6f %3d %3d %11.6f %11.6f\n",
			       e[nelem].z, ion, lambda, f, g_i, g_k, energy_i,
			       energy_k);
		    }
		  if (lambda > 0)
		    {
		      fprintf (optr,
//"Line %d %d %f %f %d %d %f %f\n",
			       "Line %2d %2d %11.6f %11.6f %3d %3d %11.6f %11.6f\n",
			       e[nelem].z, ion, lambda, f, g_i, g_k, energy_i,
			       energy_k);
		    }
		  nlines++;
		}

	    }

	}
      else
	{
//                      printf("Failed ion match %d < %d < %d\n",
//                              e[nelem].ion_min,ion,e[nelem].ion_max);
	}


    }

  printf ("Grand total nlines %d\n", nlines);

  return EXIT_SUCCESS;
}
