/********************************************************//**
 * @file   reverb.c
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
#include "atomic.h"
#include "python.h"


char	delay_dump_file[LINELENGTH];
int		delay_dump_bank_size = 65535, delay_dump_bank_curr = 0;
PhotPtr	delay_dump_bank;
int     *delay_dump_bank_ex;

/********************************************************//*
 * @name 	delay_to_observer
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
delay_to_observer(PhotPtr pp)
{
	plane_dummy	observer;
	observer.x[0] = pp->lmn[0] * geo.rmax;
	observer.x[1] = pp->lmn[1] * geo.rmax;
	observer.x[2] = pp->lmn[2] * geo.rmax;
	observer.lmn[0] = pp->lmn[0];
	observer.lmn[1] = pp->lmn[1];
	observer.lmn[2] = pp->lmn[2];
	return (pp->path + ds_to_plane(&observer, pp));
}

/********************************************************//*
 * @name 	delay_dump_prep
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
delay_dump_prep(char filename[], int restart_stat, int i_rank)
{
	FILE *fopen(), *fptr;
	char string[LINELENGTH], c_file[LINELENGTH], c_rank[LINELENGTH];
	int	i;

	//Get output filename
	strcpy(c_file, filename);		//Copy filename to new string
	strcat(c_file, ".delay_dump");
	if (i_rank > 0) 
	{
		sprintf(c_rank, "%i", i_rank);	//Write thread to string
		strcat(c_file, c_rank);			//Append thread to filename
	}
	strcpy(delay_dump_file, c_file);	//Store modified filename for later

	//Allocate and zero dump files and set extract status
	delay_dump_bank = (PhotPtr) calloc(sizeof(p_dummy), delay_dump_bank_size);
	delay_dump_bank_ex = (int *)calloc(sizeof(int), delay_dump_bank_size);
	for (i = 0; i < delay_dump_bank_size; i++)
		delay_dump_bank_ex[i] = 0;

	if (restart_stat == 1) 
	{	//Check whether the output file already has a header
		Log("delay_dump_prep: Resume run, skipping writeout\n");
		return (0);
	}

	if ((fptr = fopen(delay_dump_file, "w")) == NULL) 
	{	//If this isn't a continue run, prep the output file
	
		Error("delay_dump_prep: Unable to open %s for writing\n",
		      delay_dump_file);
		exit(0);
	}

	if (i_rank > 0) 
	{
		fprintf(fptr, "# Delay dump file for slave process %d\n", i_rank);
	}
	else
	{	// Construct and write a header string for the output file
		fprintf(fptr, "# Python Version %s\n", VERSION);
		get_time(string);
		fprintf(fptr, "# Date	%s\n#  \n", string);
		fprintf(fptr, "# \n# Freq      Wavelength  Weight   "
			" Last X     Last Y     Last Z    "
		    " Scatters   RScatter   Delay      Extracted  Spectrum   Origin   Last_Res  \n");
	}
	fclose(fptr);
	return (0);
}

/********************************************************//*
 * @name 	delay_dump_finish
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
	if (delay_dump_bank_curr > 0) 
	{
		delay_dump(delay_dump_bank, delay_dump_bank_curr - 1, 1);
	}
	free(delay_dump_bank);
	free(delay_dump_bank_ex);
	return (0);
}

/********************************************************//*
 * @name 	delay_dump_combine
 * @brief	Prepares delay dump output file
 *
 * @param [in] i_ranks		Number of parallel processes
 * @return 					0
 *
 * Collects all the delay dump files together at the end. 
 * Called by the master thread. Uses 'cat' for simplicity.
 *
 * @notes
 * 6/15	-	Written by SWM
***********************************************************/
int
delay_dump_combine(int i_ranks)
{
	FILE *fopen();//, *f_base, *f_cat;
	char c_call[LINELENGTH]; //, c_cat[LINELENGTH], c_char;
	//int i;
/*
	f_base = fopen(delay_dump_file, 'a');
	for(i=1;i<iRanks;i++)
	{
		sprintf(c_cat, "%s%d", c_cat, i);
		if((f_cat = fopen(c_cat, 'r')) == NULL)
			Error("delay_dump_combine: Missing file %s%d", c_cat,i);
		else
			while((c_char = fgetc(f_cat)) != EOF) fputc(c_char, f_base);
		fclose(f_cat);
	}
	fclose(f_base);
*/
	//Yes this is done as a system call and won 't work on Windows machines. Lazy solution!
	sprintf(c_call, "cat %s[0-9]* >> %s", delay_dump_file, delay_dump_file);
	if (system(c_call) < 0)
	{
		Error("delay_dump_combine: Error calling system command '%s'", c_call);
	}
	else
	{
		sprintf(c_call,"rm %s[0-9]*");
	}
	return (0);
}

/********************************************************//*
 * @name 	delay_dump
 * @brief	Dumps tracked photons to file
 *
 * @param [in] np			Pointer to photon array tp dump
 * @param [in] np			Number photons in p
 * @param [in] iExtracted	Array of extract flags for p
 * @return 					0
 *
 * Sifts through the photons in the delay_dump file, checking
 * if they've scattered or were generated in the wind and so
 * contribute to the delay map. Uses the same filters as 
 * the spectra_create() function for scatters & angles.
 *
 * @notes
 * 6/15	-	Written by SWM
***********************************************************/
int
delay_dump(PhotPtr p, int np, int iExtracted)
{
	FILE *fopen(), *fptr;
	int	 nphot, mscat, mtopbot, i, subzero;
	double zangle, delay;
	subzero=0;

	printf("delay_dump: Dumping %d photons\n",np);
	/*
	 * Open a file for writing the spectrum
	 */
	if ((fptr = fopen(delay_dump_file, "a")) == NULL) 
	{
		Error("delay_dump: Unable to reopen %s for writing\n",
		      delay_dump_file);
		exit(0);
	}
	for (nphot = 0; nphot < np; nphot++) 
	{
		if (iExtracted  || 
			(p[nphot].istat == P_ESCAPE && 
			(p[nphot].nscat > 0 || p[nphot].origin == PTYPE_WIND || p[nphot].origin == PTYPE_WIND_MATOM)))
		{
			zangle = fabs(p[nphot].lmn[2]);
			/*
			 * Complicated if statement to allow one to choose
			 * whether to construct the spectrum from all photons
			 * or just from photons which have scattered a
			 * specific number of times.  01apr13--ksl-Modified
			 * if statement to change behavior on negative
			 * numbers to say that a negative number for mscat
			 * implies that you accept any photon with |mscat| or
			 * more scatters
			 */
			for (i = MSPEC; i < nspectra; i++) 
			{
		  		if (((mscat = xxspec[i].nscat) > 999 ||
		       		p[nphot].nscat == mscat ||
		       		(mscat < 0 && p[nphot].nscat >= (-mscat)))
		      		&& ((mtopbot = xxspec[i].top_bot) == 0
			  		|| (mtopbot * p[nphot].x[2]) > 0))
				{
					if (iExtracted ||
						(xxspec[i].mmin < zangle && zangle < xxspec[i].mmax))
					{
						delay = (delay_to_observer(&p[nphot]) - geo.rmax) / C;
						if(delay<0)subzero++;

						fprintf(fptr,
							"%10.5g %10.5g %10.5g %+10.5g %+10.5g %+10.5g %3d     %3d     %10.5g %5d %5d %5d %10d\n",
							p[nphot].freq, C * 1e8 / p[nphot].freq, p[nphot].w, 
							p[nphot].x[0], p[nphot].x[1], p[nphot].x[2],
							p[nphot].nscat, p[nphot].nrscat, delay,
							(iExtracted ? delay_dump_bank_ex[nphot] : 0),
							i - MSPEC, p[nphot].origin, 
							p[nphot].nres);
					}
				}
			}
		}
	}

	if(subzero>0)
	{
		Error ("delay_dump: %d photons with <0 delay found! Increase path bin resolution to minimise this error.", subzero);
	}
	fclose(fptr);
	return (0);
}

/********************************************************//*
 * @name 	delay_dump_single
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
delay_dump_single(PhotPtr pp, int extract_phot)
{
	stuff_phot(pp, &delay_dump_bank[delay_dump_bank_curr]); 		//Bank single photon in temp array
		delay_dump_bank_ex[delay_dump_bank_curr] = extract_phot; 	//Record if it's extract photon
	if (delay_dump_bank_curr == delay_dump_bank_size - 1)			//If temp array is full
	{
		delay_dump(delay_dump_bank, delay_dump_bank_size, 1);
		delay_dump_bank_curr = 0;										//Dump to file, zero array position
		int		i;
		for (i = 0; i < delay_dump_bank_size; i++)
			delay_dump_bank_ex[i] = 0;									//Zero extract status of single array
	}
	else
	{
		delay_dump_bank_curr++;
	}
	return (0);
}