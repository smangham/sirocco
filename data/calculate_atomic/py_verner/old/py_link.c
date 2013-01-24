
/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

py_link_lines_levs  reads the outputs of py_read_kurucz or py_read_verner
and creates a new file that associates a level number with the upper and 
lower level for each transition.

  Description:	

This program was created in an attempt to create the basic data for an 
enhanced version of python which would cleanly handle excited state lines.
It was intended to be a step toward developing a truly non-LTE version of
python, in which the occupation numbers of certain important lines were 
calculated in the code.

Nothing produced by this code is used in python...up to at least version python_32.

  Arguments:

  Returns:

  Notes:

  History:
	00oct	ksl	Coded as py_finish_kuruzc.c
	01nov	ksl	Generalized slightly so it was obvious that this
			would work with any appropriately formatted input
			files.


 ************************************************************************/

#include	<math.h>
#include	 <stdio.h>
#include        <strings.h>

#define HC	1.98587e-16
#define C	2.997925e10
#define HEV 	4.13620e-15	/* Planck's constant in eV */
#define BOLTZMANN 1.38062e-16



#define LINELENGTH 160
#define KLINELENGTH 160


struct kurucz
  {
    double freq, gf;
    int elem, ion;
    double e1, e2;
    double j1, j2;
  }
kz;

#define NIONS       500
#define NLEVELS   15000
#define NLINES    100000

struct ions
  {
    int z;
    int istate;
    int firstlevel;
    int nlevels;
  }
ion[NIONS];
int nion;

struct configurations
  {
    int z;
    int i;
    int nlev;
    double ex;			/*excitation energy of level */
    double g;			/* multiplicity for level */
  }
level[NLEVELS];
int nlevel;

struct lines
  {
    int z, i, g1, g2;
    double w, f, e1, e2;
    int n1, n2;
  }
xl[NLINES];
int nlines;




main (argc, argv)
     int argc;
     char *argv[];
{

  FILE *levptr, *linptr, *optr, *fopen ();
  char jline[KLINELENGTH];
  char word1[30], word2[30];
  int iielem, iion;
  char ichoice;
  int ilev, g;
  float elevel;
  int g1, g2;
  float wav, f, e1, e2;
  int ion_z, ion_state;
  int i, j;
  int jion, jlevel, lastlevel;
int kkk;



  char dummy[LINELENGTH];
  char root[LINELENGTH], input[LINELENGTH];
  char ionfile[LINELENGTH], linefile[LINELENGTH], levelfile[LINELENGTH];
  char olinefile[LINELENGTH];



  printf
    ("This program produces atomic data from the kurucz list of lines\n");

  /* Determine whether input data should be read from a file or from the terminal.  */
  if (argc == 2)
    {
      strcpy (dummy, argv[1]);
    }
  else
    {
      printf ("Input file (interactive=stdin):");
      fgets (dummy, LINELENGTH, stdin);
    }
  get_root (root, dummy);

  if (strncmp (root, "dummy", 5) == 0)
    {
      printf
        ("Proceeding to create rdpar file in dummy.pf, but will not run prog\n");
    }
  else if (strncmp (root, "stdin", 5) == 0 || strncmp (root, "rdpar", 5) == 0
           || root[0] == ' ' || strlen (root) == 0)
    {
      strcpy (root, "mod");
      printf
        ("Proceeding in interactive mode\n Output files will have rootname mod\n");
    }
  else
    {
      strcpy (input, root);
      strcat (input, ".pf");

      opar (input);
      printf ("Reading data from file %s\n", input);

    }


/* Initialize variables */

  strcpy (linefile, "lines.out");
  strcpy (levelfile, "levels.out");
  strcpy (olinefile, "lines_final.out");

/* Get inputs */

  rdstr ("ion_file", levelfile);
  rdstr ("input_line_file", linefile);
  rdstr ("output_line_file",olinefile);

// Got aall the inputs
  cpar ("py_link.pf");

/* Initialize variables */
  nion = 0;
  nlevel = 0;


  levptr = fopen (levelfile, "r");
  linptr = fopen (linefile, "r");
  optr = fopen (olinefile, "w");

// Read the level array 
  while (fgets (jline, KLINELENGTH, levptr) != NULL)
    {
      sscanf (jline, "%s", word1);
      ichoice = 'z';
      if (strncmp (word1, "Ion", 3) == 0)
	ichoice = 'i';
      if (strncmp (word1, "Level", 4) == 0)
	ichoice = 'l';

      switch (ichoice)
	{
	case 'i':
	  sscanf (jline, "%s %s %d %d", word1, word2, &iielem, &iion);
	  ion_z = ion[nion].z = iielem;
	  ion_state = ion[nion].istate = iion;
	  ion[nion].firstlevel = nlevel;
	  ion[nion].nlevels = 0;
	  nion++;
	  break;
	case 'l':
	  sscanf (jline, "%s %d %d %d %d %e",
		  word1, &iielem, &iion, &ilev, &g, &elevel);
kkk=0;
while((ion[kkk].z!=iielem||ion[kkk].istate!=iion) && kkk<nion) kkk++;
if(kkk<nion){
          if(ion[kkk].nlevels==0){
              ion[kkk].firstlevel=nlevel;
                }
          ion[kkk].nlevels += 1;

	  level[nlevel].z = iielem;
	  level[nlevel].i = iion;
	  level[nlevel].nlev = ilev;
	  level[nlevel].ex = elevel;	/*excitation energy of level */
	  level[nlevel].g = g;	/* multiplicity for level */
	  nlevel++;
}
	  break;
	default:
	  break;
	}
    }

// Read the line array
  while (fgets (jline, KLINELENGTH, linptr) != NULL)
    {
      sscanf (jline, "%s", word1);
      ichoice = 'z';
      if (strncmp (word1, "Line", 4) == 0)
	ichoice = 'x';

      switch (ichoice)
	{

	case 'x':
	  sscanf (jline, "%s %d %d %e %e %d %d %e %e",
		  word1, &iielem, &iion, &wav, &f, &g1, &g2, &e1, &e2);
	  xl[nlines].z = iielem;
	  xl[nlines].i = iion;
	  xl[nlines].w = wav;
	  xl[nlines].f = f;
	  xl[nlines].e1 = e1;
	  xl[nlines].g1 = g1;
	  xl[nlines].e2 = e2;
	  xl[nlines].g2 = g2;
	  if (e1 >= 0 && e2 >= 0)
	    nlines++;
	  break;
	default:
	  break;
	}
    }

  for (i = 0; i < nion; i++)
    printf ("%d %d %d %d \n",
	    ion[i].z, ion[i].istate, ion[i].firstlevel, ion[i].nlevels);

// Now attache the number of the level to the transition //

  for (i = 0; i < nlines; i++)
    {
      iielem = xl[i].z;
      iion = xl[i].i;
      e1 = xl[i].e1;
      e2 = xl[i].e2;
      g1 = xl[i].g1;
      g2 = xl[i].g2;

      // Find the right ion 
      jion = 0;
      while ((iielem != ion[jion].z || iion != ion[jion].istate) && jion < nion)
	jion++;
      if (jion == nion)
	{
	  printf ("Failed to find ion for this transition %d %d\n", iielem, iion);
	  exit (0);
	}
      jlevel = ion[jion].firstlevel;
      lastlevel = jlevel + ion[jion].nlevels;
//      while ((e1 != level[jlevel].ex || g1 != level[jlevel].g) && jlevel < lastlevel)
      while ((e1 != level[jlevel].ex ) && jlevel < lastlevel)
	jlevel++;
      if (jlevel == lastlevel)
	{
	  fclose (optr);
	  printf ("Failed to match e1 %d %d %10.6f %d\n", iielem, iion, e1,g1);
	  exit (0);
	}
      xl[i].n1 = jlevel - ion[jion].firstlevel;

      jlevel = ion[jion].firstlevel;
      lastlevel = jlevel + ion[jion].nlevels;
//      while ((e2 != level[jlevel].ex || g2 != level[jlevel].g) && jlevel < lastlevel)
      while ((e2 != level[jlevel].ex ) && jlevel < lastlevel)
	jlevel++;
      if (jlevel == lastlevel)
	{
	  fclose (optr);
	  printf ("Failed to match e2 %d %d %10.6f %d\n", iielem, iion, e2,g2);
	  exit (0);
	}
      xl[i].n2 = jlevel - ion[jion].firstlevel;

      fprintf (optr, "Line %2d %2d %11.6f %9.6f %3d %3d %4d %4d\n",
	       xl[i].z, xl[i].i, xl[i].w, xl[i].f, xl[i].g1, xl[i].g2, xl[i].n1, xl[i].n2);
    }


}
