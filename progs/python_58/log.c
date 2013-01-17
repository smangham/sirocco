/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	These are a simple series of routines designed to store comments and errors
	in a diagnostic file.
		
		int Log_init(filename)				Open a logfile

		int Log ( char *format, ...)			Send a message to the screen and logfile
		int Log_silent ( char *format, ...)		Send a message to the logfile
		int Error ( char *format, ...)			Send a message prefaced by the word "Error:" to
										the screen and to the logfile.
		int Error_silent ( char *format, ...)		Send a message prefaced by the word "Error:" to
										the logfile
		
		int Log_close()						Close the current logfile



Arguments:		


Returns:
 
Description:	
	
	Normally, one would begin by issuing a Log_init command.  This will open a
	file which will be used for logging.  If the user does not issue a Log_init command
	but attempts to use one of the other routines in the file, then Log_init will
	be issued internally but the and the which will used for logging will be called
	"logfile".
	
	All of the Log... and Error... allow one to send what is essentially fprintf and/if desired
	printf commands.  They are designed to handle variable numbers of arguments.  So
	for example if you want to send a message to the log file only, which includes the
	variables i, and j, one would say:
		
		Log_silent("This message writes i (%d) and j (%d) to the screen\n",i,j);
		
	Log_close will close the existing log file.  If one actually uses Log_close, one should
	be sure not to try and reopen the same file for logging because this will overwrite
	the original file.
	
	

Notes:
	These could be a used more generally than in python.  The only reason python specific
	headers are called is due to the fact that for working on the MAC, I have macros that
	replace fopen so that the files will be located in a special directory.

History:
  	98feb	ksl	Coding of these subroutines began.
	99dec	ksl	Rewrote sane_check for linux using routine finate

 
**************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "atomic.h"
#include "python.h"		/* Needed so that fopen will be redefined on Mac */


FILE *diagptr;
int init_log = 0;

int
Log_init (filename)
     char *filename;
{
  FILE *fopen ();

  if ((diagptr = fopen (filename, "w")) == NULL)
    {
      printf ("Yikes: could not even open log file %s\n", filename);
      exit (0);
    }
  init_log = 1;
  return (0);
}

int
Log_close ()
{
  fclose (diagptr);
  init_log = 0;
  return (0);
}

int
Log (char *format, ...)
{
  va_list ap;
  int result;

  if (init_log == 0)
    Log_init ("logfile");

  va_start (ap, format);
  result = vprintf (format, ap);
  result = vfprintf (diagptr, format, ap);
  va_end (ap);
  return (result);
}

int
Log_silent (char *format, ...)
{
  va_list ap;
  int result;

  if (init_log == 0)
    Log_init ("logfile");

  va_start (ap, format);
  result = vfprintf (diagptr, format, ap);
  va_end (ap);
  return (result);
}

int
Error (char *format, ...)
{
  va_list ap;
  int result;

  if (init_log == 0)
    Log_init ("logfile");

  printf ("Error: ");
  va_start (ap, format);
  result = vprintf (format, ap);
  fprintf (diagptr, "Error: ");
  result = vfprintf (diagptr, format, ap);
  va_end (ap);
  return (result);
}

int
Error_silent (char *format, ...)
{
  va_list ap;
  int result;

  if (init_log == 0)
    Log_init ("logfile");

  fprintf (diagptr, "Error: ");
  va_start (ap, format);
  result = vfprintf (diagptr, format, ap);
  va_end (ap);
  return (result);
}


int
sane_check (x)
     double x;
{
  int i;
  if ((i = finite (x)) == 0)
    {
      Error ("sane_check: %d %e\n", i, x);
      return (-1);
    }
  return (0);
}
