#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "log.h"
#define LINELEN		256
#define	OLD		100
#define	NORMAL		1
#define REISSUE		-199


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
 
This file contains various routines which implement a relatively straightforward
interface to get data either from the user or a command line.

The basic routines are as follows:
	rdstr(question,answer)	gets a single contiguous string 
	rdchar(question,answer)	gets a single character
	rdint(question,answer)	gets a single integer
	rdflo(question,answer)	gets a single precision floating point number
	rddoub(question,answer)	gets a double precision floating point number
	rdline(question,answer)	gets an entire line
	rdpar_comment(string)   Adds a comment to any output.pf files



The file also contains routines to control how rdpar operates:

opar(infilename)   --- allows one to read from a file instead of from the
	command line.  If opar is not called then interactive mode is 
	assumed.  If the filedoes not exist, the routines drop into interactive
	mode. 
cpar(outfilename)  ---- causes the temporary file being created by the various
	routines to be copied to outfilename.  If there was a file with
	that name previously it will be copied to filename.old.   If cpar
	is not called, then the datat which is aquired with the readpar
	routines will be found in "tmp.rdpar"
int get_root(root,total) -- strips off the trailing characters in a file
	name.
int rdpar_set_mpi_rank(rank) -- passes the rank of the parallel process to rdpar

The file also contains a little subroutine to put out a message of
one or more lines.

message(string)

Arguments:

	For the io routines the arguments are

	char question[];	The string printed to the screen to prompt the user
					or read from the file and matched
	*answer			The value returned; the type depends on the call.
		

Returns:
	The individual io routines return 
		1 or more for a successful read
		EOF if receive done or EOF 
 
Description:	
	
	In interactive mode, the programs print out the question and
	expects an answer.  Only one integer (or whatever) is read at a time.
	The calls are generally of the form
		rdxxx(question,whatever)
		char question[]
		whatevertype *whatever
	To keep the current value, respond with \n
	In interactive mode, type ^D to exit the program immedieately
	
	In noninteractive mode, the program reads a line of a file.  It uses
	the first uninterrupted string to assure itself that it has read the correct
	line of the file (it expects the first string to be the question)
	and then the second uninterrupted string as the value
	of answer.  (An exception is rdline in which case it returns the entire
	rest of the line).
	
	In noninteractive mode, if the second string in the file is "$" (or 
	there is no second string, it will
	ask the user that question interactively and then continue to read
	other answers from the file.  (Thus if one wants to run the program
	a bunch of times, changind only one or two of the values in the
	program, one can include most of the values in the file, but leave
	the values which will change as $)

	In noninteractive mode, if the string in the file and the string in the
	question do not match, it will read through the file until one does
	match or until it reaches the EOF in which case it will switch to 
	interactive mode and try that.


	When used as an interactive interface, this
	version generates a file which can be used to run the program again.  

	That file can be edited, but one needs to stick fairly closely to the file
	which is produced.  !!!!  	
		
Notes:

	The original rdpar was written at SAO for handling the reduction programs
	assocated with Einstein data.  C and Fortran versions of the code exist.
	The SAO version of rdpar always works from parameter (".pf" files".  
	When I began working with c programs on Macs (and Atari's), I wrote
	a set of routines which I called "irdpar" which were intended to provide
	a simple interactive version of rdpar, i.e. a very readable way to get
	single input values into a program.  Eventually, it became clear that I
	was writing programs that would benefit from a file i/o interface as
	well and so this version of rdpar has been written.  
	
	This version of rdpar does not have all of the advantages of the old
	rdpar.  It does allow one to use $ in file to switch to interactive mode
	for that particular input.  But this version of rdpar cannot read more 
	than one variable from a line.  Finally, unlike rdpar this routine expects 
	variables to be in the file in the prescribed manner.  It bails out
	if the variables are not in the prescibed order.
	
	Note that you can have extra lines in the input file.  These will simply
	be copied to the output file.  But if you are thinking of using several
	programs which have some of the same keywords, be sure to put them
	in the same order in both programs.  It is OK to have 
		key1
		newkey
		key2
	
	where newkey is used by the second program but not the first.  But
	it is not a good idea to look for key2 before key1 in the second program.
		
	
	If you are writing programs using rdpar, and you seem to be having
	problems with the inputs, the usual problem is that you did not pass
	a pointer to rdxxx.
	

History:
	96dec     ksl	Began work to modifiy old interactive version of rdpar 
 				("irdpar") to function as a simpler version of "rdpar". 
 				Coded and debugged as part of Monte and Python efforts.  
	97mar   ksl  	Modified to allow one to use $ to switch to interactive mode
				for one line.
	97apr  ksl  	Modified to cause rdpar to continue reading through a 
				file until it finds the correct question keyword.  
				Thus keywords can be skipped.  Note that it
                		does not read through the file multiple times as 
                		the original rdpar would.
	97aug	ksl	Made changes needed so that rdpar would operate on the
         			Mac in a directory other than the directory in which the
         			program had been compiled.  (If on the Mac,  then rdpar
         			must be compled with the routines in "mac_file" and MAC
         			should be set to !0.)
	98feb	ksl	Made changes so that extra lines and comments would all go
         			to output file.  
	99jul	ksl	Eliminated all calls to gets in favor of fgets.  Note that fgets
				handles \n differently than does gets.  Also added a routine
				get_root to strip off the trailing characters in a file name.
				Removed special commands so would work on MAC (with Symentec C).
				These commands had only to do with the way files were opened since
				MAC assumes files are local to the location of the program rather than
				the place where the program is being executed.  I have also added restart_pf
                                which is designed to allow one to read the file multiple times.  
				There are interactions between cpar and 
                                restart_par which still need to be worked out better.  At present cpar needs
				to come after restart_par.
	00dec	ksl	Found small problem in rdchar.  It would not read properly from a file.  
			Hopefully it still works in the interactive situation.
	02jun	ksl	Modified behavior so only compares the input string to the first (.  This is
			to help make it possible to add new options to programs without forcing one
			to recreate parameter files.  
	04may	ksl	Corrected a problem where the first word ended up to have length 0.
	05jan	ksl	Corrected a problem where the first word matches a desired variable, but there
			was no second variable.  The solution is to ask for the variable on the command
			line.
        06jul	ksl	Increased LINELEN to 256
	07may	ksl	Made changes to make robust on bluegill (MacOSX)
	09mar	ksl	Added a capability to record all of the accepted values and write them out
			in a specified format to an already open file.  The purpose of this was
			to allow me to record the rdpar information in an output file, e.g a spectrum
	13jul	jm	Modified so in parallel mode only print to screen if thread 0 when reading from
			pf file. Also will not print to screen in a certain verbosity setting.
	16sept	ksl	Split string and file processing into separate subroutines, and modified
			so that the rdpar file is read in.  At present, if all of the keywords are
			contained in the .pf file, they will be read even if they are in the wrong
			order.  A line will not be read multiple times.  Also it is now possible
			to add a comment to a rdpar output file.  The behavior of the routine is also
			changed in the sense that if a variable is missing from the .pf file, one
			will be asked for the value of this variable only, and then the program will
			go on parsing the .pf file. (Note - removed restart_par as it is no longer needed)

**************************************************************/


FILE *rdin_ptr, *rdout_ptr;     /* Pointers to the input and output files */
int rdpar_stat = 0;             /*rdpar_stat=0 initially, 
                                   rd_stat=1 while reading from stdin,
                                   rdpar_stat=2 if reading from file */
int verbose = 1;                /* verbose = 1 -> printout all errors
                                   verbose=2 -->suppresses some error messages */

char current_filename[LINELEN];

int rd_rank = 0;                // rank of mpi process, set to zero

/* Array to record the actual values of rdpar accepted by rdpar.  All values are 
 * stored as strings */
#define  MAX_RECORDS   500

int rdpar_nrec = 0;

struct rdpar
{
  char name[LINELEN];
  char value[LINELEN];
}
rdpar_record[MAX_RECORDS];

/* Struture to hold the input rdpar file */

#define UNUSED 0
#define USED  1

int rdpar_cursor = 0;
int rdpar_ntot = 0;             // total number of raw input lines
struct rdpar_raw
{
  char line[LINELEN];           // Raw input line 
  int icheck;                   // Check off to indicate a line has been used
}
input[MAX_RECORDS];



/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

	opar sets rdpar up to read from a file.  If it fails you proceed interactively


  Description:	


  Returns:

  Notes:

 

  History:
  	16sept	ksl	Standard header included

 ************************************************************************/

int
opar (filename)
     char filename[];
{
  FILE *fopen (), *tmp_ptr;
  int rdpar_init ();

  /* Check that an input file is not currently open */
  if (rdpar_stat == 2)
  {
    printf ("Error: opar: A file for reading input is already open.\n");
    return (rdpar_stat);
  }

  /*Open a temporary outbut file */
  rdpar_init ();
  if ((tmp_ptr = fopen (filename, "r")) == NULL)
  {
    printf ("Error: opar: Could not open filename %s\n", filename);
    printf ("               Proceeding in interactive mode\n");
    return (rdpar_stat);
  }
  else
  {
    rdin_ptr = tmp_ptr;
    rdpar_stat = 2;             /* implies we are now trying to read from a file */
    while (fgets (input[rdpar_ntot].line, LINELEN, rdin_ptr) != NULL)
    {
      input[rdpar_ntot].icheck = UNUSED;
      rdpar_ntot++;
    }
    fclose (rdin_ptr);
  }
  strcpy (current_filename, filename);


  return (rdpar_stat);
}



/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

	add_par routine appends data from another parameter file to the existing values.


  Description:	


  Returns:

  Notes:
	This is intended for the case where one would like to use a second parameter
	file for extra parameters or to add a parameter to an old file. 

	Note that if cpar is called, it will write all of the parameters in the order
	they were used.
 

  History:
  	16sept	ksl	Coded                     

 ************************************************************************/


int
add_par (filename)
     char filename[];
{
  int rdpar_init ();

  /* Check that an input file is not currently open */
  if (rdpar_stat != 2)
  {
    printf ("Error: add_par: opar needs to be called before add_par\n");
    return (rdpar_stat);
  }

  if ((rdin_ptr = fopen (filename, "r")) == NULL)
  {
    printf ("Error: add_par: Could not additional file %s\n", filename);
    return (rdpar_stat);
  }

  while (fgets (input[rdpar_ntot].line, LINELEN, rdin_ptr) != NULL)
  {
    input[rdpar_ntot].icheck = UNUSED;
    rdpar_ntot++;
  }
  fclose (rdin_ptr);


  return (rdpar_stat);
}


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

	This routine should be called to copy the tmp.rdpar to filename.  
   	The old parameters of the input file are stored in filename.old.  Also has the 
   	effect of closing all the files associated with rdpar

  Description:	


  Returns:

  Notes:

 

  History:
  	16sept	ksl	Standard header included

 ************************************************************************/


int
cpar (filename)
     char filename[];
{
  char old_filename[LINELEN];

  if (rd_rank != 0)
    return 0;

  fclose (rdout_ptr);
  rdpar_stat = 0;


  strcpy (old_filename, filename);
  strcat (old_filename, ".old");

  if (remove (old_filename) != 0 && verbose)
    printf ("Old file %s did not exist\n", old_filename);

  if (rename (filename, old_filename) != 0 && verbose)
    printf ("Could not rename %s to %s\n", filename, old_filename);

  if (rename ("tmp.rdpar", filename) != 0 && verbose)
    printf ("Could not rename %s to %s", "tmp.rdpar\n", filename);

  printf ("A new .pf file %s has been written\n", filename);
  return (NORMAL);
}

/* Despite its name, this routine is internal to the other rdpar parts of rdpar. It
   does however set put rdpar into a known mode and set up the output file */

int
rdpar_init ()
{
  FILE *fopen ();
  rdin_ptr = stdin;             /* Initialize rdin_ptr to standard input */
  if ((rdout_ptr = fopen ("tmp.rdpar", "w")) == NULL)
  {
    printf ("Error: rdpar_init: Problem opening tmp.rdpar\n");
    exit (0);
  }
  rdpar_stat = 1;
  strcpy (current_filename, "tmp.rdpar.out");
  return (NORMAL);
}

/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

	string_process is the routine called by all of the individual
	routines, like rdpar for processing a line of input whether
	it comes from the command line or from a .pf file
	
 
  Description:	
  	At present this is merely a steering routine


  Returns:
	When the routine returns NORMAL or OLD the input has been 
	successfully processed

  Notes:

  History:
  	16sept	ksl	Moved the command line processing and the file
			processing into separate routines below

 ************************************************************************/

int
string_process (question, dummy)
     char question[], dummy[];
{

  if (rdpar_stat == 0)
    rdpar_init ();              /* Set rdin_ptr to stdin, and rdout_ptr to file tmp.rdpar */


  if (rdpar_stat == 1)
  {                             /* Then we are reading from stdin */
    return (string_process_from_command_line (question, dummy));
  }


  else
  {                             /*Then read from the file */

    return (string_process_from_file (question, dummy));

  }
}

/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

	string_process_from the command line, handles the case
	where one is processing a line from the command line
	
 
  Description:	


  Returns:
	When the routine returns NORMAL or OLD the input has been 
	successfully processed

  Notes:

  History:
  	16sept	ksl	Moved from string_proces above.

 ************************************************************************/


int
string_process_from_command_line (question, dummy)
     char question[], dummy[];
{
  char tdummy[LINELEN];
  fprintf (stderr, "%s (%s) :", question, dummy);
  fflush (stderr);
  strcpy (tdummy, "");
  if (fgets (tdummy, LINELEN, stdin) == NULL)
  {
    printf ("Exiting since rdpar got EOF in interactive mode\n");
    exit (0);
  }
  else if (tdummy[0] == '\n')
  {                             //Use the current value

    printf ("%s %s\n", question, dummy);
    fprintf (rdout_ptr, "%-30s %20s\n", question, dummy);
    rdpar_store_record (question, dummy);
    return (OLD);
  }
  else if (tdummy[0] == '!')
  {
    system (&tdummy[1]);        /* Send a command to the system */
    return (REISSUE);
  }
  else if (strncmp (tdummy, "done", 4) == 0)
    return (EOF);
  else
  {                             //Use the value input from the command line

    strcpy (dummy, tdummy);
    fprintf (rdout_ptr, "%-30s %20s", question, dummy);
    rdpar_store_record (question, dummy);
    return (NORMAL);
  }
}


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

	string_process_from_file handles the case
	where one is processing a line from the a file
	
 
  Description:	

  	Basically this routine reads information from the structure
	called input, which was filed by calls now contained in
	opar or rdpar.  


  Returns:
	When the routine returns NORMAL or OLD the input has been 
	successfully processed

  Notes:

  	When a keyword is missing the routine invokes the inteactive version
	of rdpar, but it does not force all future inputs to come from
	the command line version, as was the case for version of rdpar
	prior to 16sept.

  History:
  	16sept	ksl	Moved from string_proces above, and significantly
			modified.
	16oct	ksl	Added a call to check_synonyms, which looks for
			names that we have changed.

 ************************************************************************/

int
string_process_from_file (question, dummy)
     char question[], dummy[];
{

  char firstword[LINELEN], secondword[LINELEN];
  char *line;
  char *ccc, *index (), *fgets ();
  int nwords, wordlength;
  char old_question[LINELEN];

  for (rdpar_cursor = 0; rdpar_cursor < rdpar_ntot; rdpar_cursor++)
  {
    if (input[rdpar_cursor].icheck == USED)
    {
      continue;
    }

    line = input[rdpar_cursor].line;

    // Parse the first two words of the input line

    strcpy (firstword, "");
    strcpy (secondword, "");
    nwords = sscanf (line, "%s %s", firstword, secondword);

    wordlength = strlen (firstword);
    if (nwords < 2 || wordlength == 0)
    {
      continue;                 // Read the next line in the file
    }

    // Strip off everthing in paren
    if ((ccc = index (firstword, '(')) != NULL)
    {
      wordlength = (int) (ccc - firstword);
      if (wordlength == 0)
      {
        continue;
      }
    }


    // If there is a match go to the section that handles the various possibilities

    if (strncmp (question, firstword, wordlength) == 0)
    {
      break;                    // We have matched the keywords
    }

    /* Check for synonyms - That is look for keywords in the file we are reading 
     * that have been replaced by a new keyword
     */

    if (check_synonyms (question, old_question) == 1 && strncmp (old_question, firstword, wordlength) == 0)
    {
      break;
    }

    // Print a warning if the lines do not match
    if (verbose > 1 && rd_rank == 0)
      printf ("Warning:  Question (%s) does not match word (%s) in file. Continuing!\n", question, firstword);
  }

  // Handle the EOF since we were not successful in identifying the keyword
  if (rdpar_cursor == rdpar_ntot)
  {
    // printf ("Error: string_proces: Unexpectedly reached EOF\n");
    printf ("Switching to interactive mode for this variable\n");
    return (string_process_from_command_line (question, dummy));
  }
  else
  {
    input[rdpar_cursor].icheck = USED;
    rdpar_cursor++;             // This is needed because we have already processed the earlier line
  }


  // At this point we know we have the correct line in the input file
  if (strncmp (secondword, "$", 1) == 0 || nwords == 1)
  {                             // This if for the case where one wants to read this value only from the command line
    fprintf (stderr, "%s (%s) :", question, dummy);
    fflush (stderr);

    strcpy (secondword, dummy);

    fgets (dummy, LINELEN, stdin);

    if (strcmp (dummy, "\n") == 0)
    {                           /* Store the provided value since \n */
      rdpar_store_record (question, secondword);
      fprintf (rdout_ptr, "%-40s   %s\n", question, secondword);
    }
    else
    {                           /* Store the value received via the command line */

      rdpar_store_record (question, dummy);
      fprintf (rdout_ptr, "%-40s   %s\n", question, dummy);
    }

    return (NORMAL);
  }
  else                          // This handles the situation where the variable is actually read from the rdpar file
    strcpy (dummy, secondword);
  rdpar_store_record (question, secondword);
  fprintf (rdout_ptr, "%-40s   %s\n", question, secondword);
  return (NORMAL);
}




/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  	rdpar_store_record stores the parsed version of an input into 
	a structure so that it can be written out to a file, generally
	the header for one of the ascii files produced by python

	
 
  Description:	

	rdpar_store_record(name,value) adds a line to the structrue rdpar_record
	which records the accepted values for all of the variables.



  Returns:

  Notes:

	rdpar_save is the routine that writes out the header

  History:
	0903	ksl	Routine added so that one could more easily copy values
		read in through rdpar to an output file

 ************************************************************************/


int
rdpar_store_record (name, value)
     char *name, *value;
{
  strcpy (rdpar_record[rdpar_nrec].name, name);
  strcpy (rdpar_record[rdpar_nrec].value, value);
  rdpar_nrec++;

  return (0);
}



/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

	rdpar_save(file_ptr)  savve the actual values input to rdpar to a file 

 
  Description:	

	saves the stored parameter file vectors to an already open 
	file, e.g. a spectrum file.  It allows one to attach the
	inputs variables to the file for future reference.


  Returns:

  Notes:
  	The input is a file_ptr because genearlly this is intended simply
	to document the inputs that were used to produce and output file, 
	in the case of python, a spectrum


  History:
	16sept	ksl	Standard header added

 ************************************************************************/


int
rdpar_save (file_ptr)
     FILE *file_ptr;
{
  int i;


  for (i = 0; i < rdpar_nrec; i++)
  {
    fprintf (file_ptr, "# Var	 %-40s	%s\n", rdpar_record[i].name, rdpar_record[i].value);
  }
  return (0);
}


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

	Add a comment to the header of the output file

 
  Description:	



  Returns:

  Notes:


  History:
	16sept	ksl	Standard header added

 ************************************************************************/

int
rdpar_comment (char *format, ...)
{
  va_list ap, ap2;
  int result = 0;

  if (rdpar_stat == 0)
    rdpar_init ();              /* Set rdin_ptr to stdin, and rdout_ptr to file tmp.rdpar */


  va_start (ap, format);
  va_copy (ap2, ap);            /* ap is not necessarily preserved by vprintf */

  if (rd_rank == 0)
  {
    fprintf (rdout_ptr, "\n### ");
    result = vfprintf (rdout_ptr, format, ap2);
    fprintf (rdout_ptr, "\n");
  }
  va_end (ap);
  return (result);
}

/*Add a general purpose message line */
int
message (string)
     char string[];
{
  fprintf (stderr, "%s\n", string);
  fflush (stderr);
  return (NORMAL);
}



/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  	rdstr
	rddhar
  	rdint
	rdint2  - will read two integers from the command line (but this is
		not integrated with the rest of rdpar
	rdflo
	rddoub
	rdline - will read an entire line, but is unclear that this is integrated
		with the rest of rdpar


	These are the specific routines that are called to read in a variable. In
	general these routiens have the form

	rdxxx (question,answer)

 
  Description:	

  	All of these routines are very similar, they are intended for different
	types of inpurs



  Returns:

  Notes:


  History:
	16sept	ksl	Standard header added

 ************************************************************************/


int
rdstr (question, answer)
     char question[], answer[];
{
  int query;
  char dummy[LINELEN];
  query = REISSUE;
  if (rdpar_stat == 0)
    rdpar_init ();              /* Set rdin_ptr to stdin, and rdout_ptr to file tmp.rdpar */
  while (query == REISSUE)
  {
    sprintf (dummy, "%s", answer);
    query = string_process (question, dummy);
    if (query == NORMAL)
    {
      if (sscanf (dummy, "%s", answer) != 1)
      {
        printf ("Could not convert input (%s) to string. Try again\n", dummy);
        query = REISSUE;
      }
      if (rd_rank == 0 && verbose == 1)
        printf ("%s   %s\n", question, answer);
    }
  }
  return (query);
}


int
rdchar (question, answer)
     char question[];
     char *answer;
{
  int query;
  char dummy[LINELEN];
  query = REISSUE;
  if (rdpar_stat == 0)
    rdpar_init ();              /* Set rdin_ptr to stdin, and rdout_ptr to file tmp.rdpar */
  while (query == REISSUE)
  {
    sprintf (dummy, "%c", *answer);
    query = string_process (question, dummy);
    if (query == NORMAL)
    {
      if (sscanf (dummy, " %1c", answer) != 1)
      {
        printf ("Could not convert input (%s) to character. Try again\n", dummy);
        query = REISSUE;
      }

      if (rd_rank == 0 && verbose == 1)
        printf ("%s   %1c\n", question, *answer);
    }
  }
  return (query);
}

int
rdint (question, answer)
     char question[];
     int *answer;
{
  int query;
  char dummy[LINELEN];
  query = REISSUE;
  if (rdpar_stat == 0)
    rdpar_init ();              /* Set rdin_ptr to stdin, and rdout_ptr to file tmp.rdpar */
  while (query == REISSUE)
  {
    sprintf (dummy, "%d", *answer);
    query = string_process (question, dummy);
    if (query == NORMAL)
    {
      if (sscanf (dummy, "%d", answer) != 1)
      {
        printf ("Could not convert input (%s) to integer. Try again\n", dummy);
        query = REISSUE;
      }
      if (rd_rank == 0 && verbose == 1)
        printf ("%s	  %d\n", question, *answer);
    }
  }
  return (query);
}

// XXX it is not obvious the next routine works.  It should be tested. 
// It's probably not actually used in Python

int
rdint2 (question, answer1, answer2)
     char question[];
     int *answer1, *answer2;
{
  int query, found = 0;
  char dummy[LINELEN];
  query = REISSUE;
  if (rdpar_stat == 0)
    rdpar_init ();              /* Set rdin_ptr to stdin, and rdout_ptr to file tmp.rdpar */
  while (query == REISSUE)
  {
    if (found == 0)
      sprintf (dummy, "%d", *answer1);
    else if (found == 1)
      sprintf (dummy, "%d", *answer2);
    query = string_process (question, dummy);
    if (query == NORMAL)
    {
      if (found == 0)
        sscanf (dummy, "%d", answer1);
      else if (found == 1)
      {
        sscanf (dummy, "%d", answer2);
        if (rd_rank == 0 && verbose == 1)
          printf ("%s	  %d %d\n", question, *answer1, *answer2);
      }
    }
  }
  return (query);
}


int
rdflo (question, answer)
     char question[];
     float *answer;
{
  int query;
  char dummy[LINELEN];
  query = REISSUE;
  if (rdpar_stat == 0)
    rdpar_init ();              /* Set rdin_ptr to stdin, and rdout_ptr to file tmp.rdpar */
  while (query == REISSUE)
  {
    sprintf (dummy, "%g", *answer);
    query = string_process (question, dummy);
    if (query == NORMAL)
    {
      if (sscanf (dummy, "%e", answer) != 1)
      {
        printf ("Could not convert input (%s) to float. Try again\n", dummy);
        query = REISSUE;
      }
      if (rd_rank == 0 && verbose == 1)
        printf ("%s	  %e\n", question, *answer);
    }
  }
  return (query);
}

int
rddoub (question, answer)
     char question[];
     double *answer;
{
  int query;
  char dummy[LINELEN];
  query = REISSUE;
  if (rdpar_stat == 0)
    rdpar_init ();              /* Set rdin_ptr to stdin, and rdout_ptr to file tmp.rdpar */
  while (query == REISSUE)
  {
    sprintf (dummy, "%g", *answer);
    query = string_process (question, dummy);
    if (query == NORMAL)
    {
      if (sscanf (dummy, "%le", answer) != 1)
      {
        printf ("Could not convert input (%s) to double. Try again\n", dummy);
        query = REISSUE;
      };
      if (rd_rank == 0 && verbose == 1)
        printf ("%s	  %e\n", question, *answer);
    }
  }
  return (query);
}

int
rdline (question, answer)
     char question[];
     char answer[];
{
  int query;
  char dummy[LINELEN];
  query = REISSUE;
  if (rdpar_stat == 0)
    rdpar_init ();              /* Set rdin_ptr to stdin, and rdout_ptr to file tmp.rdpar */
  while (query == REISSUE)
  {
    sprintf (dummy, "%s", answer);
    query = string_process (question, dummy);
    if (query == NORMAL)
    {
      strcpy (answer, dummy);
      if (rd_rank == 0 && verbose == 1)
        printf ("%s	  %s\n", question, answer);
    }
  }
  return (query);
}

/* This is the end of the various routines which parse different kinds of inputs */




/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

	get_root takes the string in "total" and constructs the "root" name for 
	a file keeping all of the charcthers up to .pf.  If the file does not
    include a .pf, the root is taken to be same as "total"


  Description:	


  Returns:

  Notes:

	The order of the arguments may be some what odd, but it resembles
	routines like strcpy in which the output is the first string

  History:
	99jul	ksl	Coded 
    16dec   ksl Adjusted so that parameter files can have periods in them.
                Previously get_root took the root 

 ************************************************************************/

int
get_root (root, total)
     char root[], total[];
{
  int j;
  char *pf;
  int position;

  /* Check whether total is an empty string */

  j = strcspn (total, "\n");
  if (j == 0)
  {
    strcpy (root, "rdpar");
    return (0);
  }

  /* Check for .pf at the end of the string 
   * Note that there is no easy way to actually
   * check for the last occurence.  Here we
   * assume there is only one occurrcne of .pf
   */

  pf = strstr (total, ".pf");
  if (pf != NULL)
  {
    position = pf - total;
    strncpy (root, total, position);
    root[position] = '\0';
    printf ("xxx %s\n", root);
    return (0);
  }

  strncpy (root, total, j);
  root[j] = '\0';


  return (0);
}





/**************************************************************************
                    Southampton University            


  Synopsis:  

	the next set of routines are designed to clean up the code 
	in parallel mode


  Description:	



  Returns:

  Notes:
	The next routine simply sets the rank of the process 
if not in parallel mode then we set rd_rank to zero

 

  History:
  	13Jul	jm	Coded

 ************************************************************************/


int
rdpar_set_mpi_rank (rank)
     int rank;
{
  rd_rank = rank;
  return (0);
}


/* sets the verbosity level
 */

int
rdpar_set_verbose (vlevel)
     int vlevel;
{
  if (vlevel < 2)
    verbose = 0;
  return (0);
}
