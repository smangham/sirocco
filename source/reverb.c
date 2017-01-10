/***********************************************************/
/** @file   reverb.c
 * @Author SWM
 * @date   July, 2015
 * @brief  Reverberation mapping functions.
 *
 * File containing reverberation mapping functions.
 ***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sqlite3.h>
#include "atomic.h"
#include "python.h"


sqlite3 *db;
int subzero;
int seconds_per_day = 86400;

/**********************************************************/
/** @name 	delay_to_observer
 * @brief	Calculates the delay to the observer plane
 *
 * @param [in] pp			Pointer to test photon
 * @return 					Distance from photon to plane
 *
 * Sets up the observer plane, then calculates the distance
 * from the current photon to it. Uses ds_to_plane() for
 * reliability.
 *
 * @notes
 * 9/14	-	Written by SWM
***********************************************************/
double
delay_to_observer (PhotPtr pp)
{
  plane_dummy observer;
  observer.x[0] = pp->lmn[0] * geo.rmax;
  observer.x[1] = pp->lmn[1] * geo.rmax;
  observer.x[2] = pp->lmn[2] * geo.rmax;
  observer.lmn[0] = pp->lmn[0];
  observer.lmn[1] = pp->lmn[1];
  observer.lmn[2] = pp->lmn[2];
  return (pp->path + ds_to_plane (&observer, pp));
}

/**********************************************************/
/** @name 	delay_dump_prep
 * @brief	Prepares delay dump output file
 *
 * @param [in] filename		File root for run
 * @param [in] restart_stat If this is a restart run
 * @param [in] i_rank		Parallel rank
 * @return 					0
 *
 * Sets up filenames, allocates bank for temporary storage
 * of photons that are to be dumped, and if it's not a resume
 * run prints out a header for the file. Filename is set per
 * thread. The file is then built up in batches using
 * delay_dump() in increments of #delay_dump_bank_size.
 *
 * @notes
 * 9/14	-	Written by SWM
***********************************************************/
int
delay_dump_prep_master (int restart_stat)
{
  char c_file[LINELENGTH], *err_msg;
  subzero = 0;

  // Open SQLite database
  strcpy (c_file, files.root);  //Copy filename to new string
  strcat (c_file, ".db");
  printf("Preparing to open database '%s'",c_file);
  if (SQLITE_OK != sqlite3_open(c_file, &db)) {
    fprintf(stderr, "Master cannot open database: %s\n", sqlite3_errmsg(db));
    sqlite3_close(db);
    exit(0);
  } 
  // If this is not a restart run, purge the table and create a new one.
  if(!restart_stat) {
	  if (SQLITE_OK != sqlite3_exec(db, "DROP TABLE IF EXISTS Photons;", 0, 0, &err_msg)) {
	    fprintf(stderr, "SQL error: %s\n", err_msg);
	    sqlite3_free(err_msg);        
	    sqlite3_close(db);
	    exit(0);
	  } 
	  if (SQLITE_OK != sqlite3_exec(db, "CREATE TABLE Photons("
	  		"id INTEGER PRIMARY KEY, Wavelength FLOAT, Weight FLOAT, X FLOAT, Y FLOAT, Z FLOAT,"
	  		"ContinuumScatters INTEGER, ResonantScatters INTEGER, Delay FLOAT, Spectrum INTEGER," 
		 		"Origin INTEGER, Resonance INTEGER, Origin_matom BOOLEAN);", 0, 0, &err_msg)) {
	    fprintf(stderr, "SQL error: %s\n", err_msg);
	    sqlite3_free(err_msg);        
	    sqlite3_close(db);
	    exit(0);
	  } 
  }
  return (0);
}
int delay_dump_prep_slave() {
  char c_file[LINELENGTH], err_msg;
  subzero = 0;

  // Open SQLite database
  strcpy (c_file, files.root);  //Copy filename to new string
  strcat (c_file, ".db");
  if (SQLITE_OK != sqlite3_open(c_file, &db)) {
    fprintf(stderr, "Slave cannot open database: %s\n", sqlite3_errmsg(db));
    sqlite3_close(db);
    exit(0);
  } 
  return(0);
}

/**********************************************************/
/** @name 	delay_dump_finish
 * @brief	Finishes dumping tracked photons to file
 *
 * @return 					0
 *
 * Dumps the remaining tracked photons to file, frees memory.
 *
 * @notes
 * 6/15	-	Written by SWM
***********************************************************/
int
delay_dump_finish (void)
{
	sqlite3_close(db);
  return (0);
}
/**********************************************************/
/** @name 	delay_dump_single
 * @brief	Preps a single photon to be dumped
 *
 * @param [in] pp			Pointer to extracted photon
 * @param [in] iExtracted	Is p an extract photon
 * @return 					0
 *
 * Takes a photon and copies it to the staging arrays for 
 * delay dumping, to be output to file later.
 *
 * @notes
 * 6/15	-	Written by SWM
***********************************************************/
int
delay_dump (PhotPtr pp, int i_spec)
{
	char *err_msg;
  
 	//If we're filtering out continuum photons and this is a continuum photon, throw it away.
  if(geo.reverb_filter_lines == -1 && pp->nres == -1) return (1);
	char c_sql[LINELENGTH];
	int b_matom = 0;
	int i_origin = pp->origin;
	if(pp->origin > 10) {
		b_matom = 1;
		i_origin -= 10;
	}
	double r_delay = (delay_to_observer (pp) - geo.rmax) / (C * seconds_per_day);
	if(r_delay<0) subzero++;

  sprintf(c_sql,"INSERT INTO Photons VALUES(%g, %g, %g, %g, %g, %d, %d, %g, %d, %d, %d, %d) ", 
  				C * 1e8 /pp->freq, pp->w, pp->x[0], pp->x[1], pp->x[2], 
  				pp->nscat-pp->nrscat, pp->nrscat, r_delay, i_spec, i_origin, pp->nres, b_matom);

	if (SQLITE_OK != sqlite3_exec(db, c_sql, 0, 0, &err_msg)) {
	      fprintf(stderr, "SQL error: %s\n", err_msg);
	      sqlite3_free(err_msg);        
	      sqlite3_close(db);
	      exit(0);
	} 
  return (0);
}
