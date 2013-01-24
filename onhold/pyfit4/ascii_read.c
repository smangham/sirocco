


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
#include "models.h"


int
ascii_read (filename)
     char *filename;
{

  FILE *fptr, *fopen ();


  /* Open the input file.  Exit if it is not opened successfully */
  if ((fptr = fopen (filename, "r")) == NULL)
    {
      Error ("Failed to open file %s\n", filename);
      exit (0);
    }

/* Read and print the input file */
  ascii_nlines = 0;
  while (fgets (ascii_store[ascii_nlines], LINELEN, fptr) != NULL
	 && ascii_nlines < MAX_LINES)
    {
      ascii_nlines++;
    }


  Log ("ascii_read: read %d from %s\n", ascii_nlines, filename);

/* OK now end program and go on to the next important task. */
  return (ascii_nlines);

}

int
check_comment (xstring, answer)
     char *xstring;
     char *answer;
{
  int nline;
  char *line, format[LINELEN];
  char word[LINELEN];


  if (ascii_nlines == 0)
    {
      Error ("check comment: no file has been read:");
      strcpy (answer, "Unknown");
      return (1);
    }

  strcpy (format, xstring);
  strcat (format, " %s");

  nline = 0;
  while (nline < ascii_nlines)
    {
      line = ascii_store[nline];
      // printf ("Parsing %s\n", line);

      if ((line = strstr (line, "#")) != NULL
	  && (line = strstr (line, xstring)) != NULL)
	{
	  //line++;             // To push past #
	  // printf ("Got here  %s\n", line);
	  if (sscanf (line, format, word) == 1)
	    {
	      strcpy (answer, word);
	      return (0);
	    }
	}
      nline++;
    }
  strcpy (answer, "Unknown");
  return (1);
}

/* Return a comment line with the entire string for that line 
 * which the reading program can parse
 *
 * Note that this always stars from the beginning of the
 * file so the idenfifier xstring must be unique
 */

int
get_comment_line (xstring, answer)
     char *xstring;
     char *answer;
{
  int nline;
  char *line, format[LINELEN];


  if (ascii_nlines == 0)
    {
      Error ("check comment: no file has been read:");
      strcpy (answer, "Unknown");
      return (1);
    }

  strcpy (format, xstring);
  strcat (format, " %s");

  nline = 0;
  while (nline < ascii_nlines)
    {
      line = ascii_store[nline];
      // printf ("Parsing %s\n", line);

      if ((line = strstr (line, "#")) != NULL
	  && (line = strstr (line, xstring)) != NULL)
	{
	  //line++;             // To push past #
	  // printf ("Got here  %s\n", line);
	  strcpy (answer, line);
	  return (0);
	}
      nline++;
    }

  strcpy (answer, "Unknown");
  return (1);
}
