


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  Description:	

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:

  This is roughly what the Makefile should look like. Uncomment the mv 
  statement to store binary in ~/bin

CC = gcc
CFLAGS = -c -g -I$$HOME/include -Wall
LDFLAGS= -L$$HOME/lib -lm -lkpar
BIN = $$HOME/bin

simple: simple.o
	gcc simple.o  $(LDFLAGS) -o simple
#	mv $@ $(BIN)

 

  History:
2004	ksl	Coded as better.c

 ************************************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


#include "log.h"

#define LINELENGTH 132
#define NPARS		10
#define NMODS		100

#define NWAVES		100

struct Model
{
  char name[LINELENGTH];
  float par[NPARS];
  float *w;
  float *f;
} mods[NMODS];


int
main (argc, argv)
     int argc;
     char *argv[];
{
  int nmod;
  float w[NWAVES], f[NWAVES];
  int n;

/* Make some models */
  for (nmod = 0; nmod < 100; nmod++)
    {
      make_fake_model (w, f);
      mods[nmod].w = calloc (sizeof (float), NWAVES);
      mods[nmod].f = calloc (sizeof (float), NWAVES);
      for (n = 0; n < NWAVES; n++)
	{
	  mods[nmod].w[n] = w[n];
	  mods[nmod].f[n] = f[n];
	}
    }
/* Now print something out */

  print_models ();



/* OK now end program and go on to the next important task. */
  exit (0);

}

int make_init = 0;

int
make_fake_model (w, f)
     float w[], f[];
{
  int n;
  n = 0;
  while (n < NWAVES)
    {
      w[n] = make_init * 100. + n;
      f[n] = make_init * 1000. + n;
      make_init++;
      n++;
    }
  return (0);
}


int
print_models ()
{
  int nmod;

  nmod = 0;
  while (nmod < NMODS)
    {

      printf ("hello %3d  %10.2f  %10.2f\n", nmod, mods[nmod].w[50],
	      mods[nmod].f[50]);
      nmod++;
    }

  return (0);
}
