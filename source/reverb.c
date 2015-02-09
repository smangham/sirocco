/* 
 * The subroutines in this file handle outputting spectrum array reverb data
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"

#include "python.h"

/***********************************************************
                                    Space Telescope Science Institute

Synopsis:
	delay_spectrum_summary(filename,mode,nspecmin,nspecmax,select_spectype,renorm,loglin)  
		writes out the spectrum to a file 

Arguments:		

	char filename[]		The name of the file to write
	char mode[];		The mode in which the file should be opened, usually 'w', but
					it could be 'a'
	int nspecmin,nspecmax	These two numbers define the spectra you want to write.    
	int select_spectype	The type of spectral file you want to create, 
					0 = raw, 1= flambda,2=fnu 
        double	renorm		This is renormalization which incrementally decreases to
				one as the detailed spectral calculation goes forward.  It
				was added to allow one to print out the spectrum at the
				end of each cycle, rather than the end of the entire
				calculation.
	char loglin[]		Are we outputting a log or a linear spectrum
				

Returns:
  
Description:	
	This simple routine simply writes the spectra to a file in an easily interpretable
	ascii format. Normally one would write all of the spectra in one go, but  one can use 
	spectrum summary to write various spectra to various files by using the variables
	nspecmin and nspecmax.. 
			
Notes:

History:
	14aug	-	Cloned from spectrum_summary by SWM
**************************************************************/
int 
delay_spectrum_summary (
    char filename[],
    char mode[],
    int nspecmin,
    int nspecmax,
    int select_spectype,
    double renorm,				/* parameter used to rescale spectrum as it is building up */
    int loglin				/* switch to tell the code if we are utputting a log or a lin */
)
{
	FILE *fopen(), *fptr;
	int i, n;
	char string[LINELENGTH];
	double freq, freqmin, dfreq, freq1;
	double lfreqmin, lfreqmax, ldfreq;
	double x, dd;



	/* 
	 * Open or reopen a file for writing the spectrum 
	 */
	if ((fptr = fopen(filename, "w")) == NULL)
	{
		Error("delay_spectrum_summary: Unable to open %s for writing\n", filename);
		exit(0);
	}

	/* 
	 * Check that nspecmin and nspecmax are reasonable 
	 */
	if (nspecmin < 0 || nspecmax < 0 || nspecmin > nspecmax)
	{
		Error("delay_spectrum_summary: nspecmin %d or nspecmax %d not reasonable \n", nspecmin, nspecmax);
		exit(0);
	}

	/* 
	 * Construct and write a header string for the output file 
	 */
	fprintf(fptr, "# Python Version %s\n", VERSION);

	get_time(string);
	fprintf(fptr, "# Date	%s\n#  \n", string);

	/* 
	 * Write the rest of the header for the spectrum file 
	 */

	fprintf(fptr, "# \n# Freq.        Lambda");
	for (n = nspecmin; n <= nspecmax; n++)
	{
		fprintf(fptr, "   %8s", xxspec[n].name);
	}
	fprintf(fptr, " Weight\n");

	/* 
	 * Don't print out the end bins because they include all photons outside the frequency range and there may be some as a result of the fact that the bb function generate some IR photons 
	 */
	dd = 4. * PI * (100. * PC) * (100. * PC);

	if (loglin == 0)
	{
		freqmin = xxspec[nspecmin].freqmin;
		dfreq = (xxspec[nspecmin].freqmax - freqmin) / NWAVE;
		for (i = 1; i < NWAVE - 1; i++)
		{
			freq = freqmin + i * dfreq;
			fprintf(fptr, "%-8e %.3f ", freq, C * 1e8 / freq);
			for (n = nspecmin; n <= nspecmax; n++)
			{
				if(xxspec[n].delay[i] == 0.) x = 0;
				else x = xxspec[n].delay[i] / xxspec[n].delay_weight[i];
				fprintf(fptr, " %10.5g %10.5g", x * renorm, xxspec[n].f[i]);
			}
			fprintf(fptr, "\n");
		}
	}
	else if (loglin == 1)
	{
		lfreqmin = log10(xxspec[nspecmin].freqmin);
		freq1 = lfreqmin;
		lfreqmax = log10(xxspec[nspecmin].freqmax);
		ldfreq = (lfreqmax - lfreqmin) / NWAVE;

		for (i = 1; i < NWAVE - 1; i++)
		{
			freq = pow(10., (lfreqmin + i * ldfreq));
			dfreq = freq - freq1;
			fprintf(fptr, "%-8e %.3f ", freq, C * 1e8 / freq);
			for (n = nspecmin; n <= nspecmax; n++)
			{
				if(xxspec[n].delay[i] == 0.) x = 0;
				else x = xxspec[n].delay[i] / xxspec[n].delay_weight[i];
				fprintf(fptr, " %10.5g", x * renorm);	/* this really shouldn't get called if we are outputting log data */
			}


			fprintf(fptr, "\n");
			freq1 = freq;
		}
	}
	fclose(fptr);
	return (0);
}


/***********************************************************
Synopsis:
	delay_dump(filename, p, f1, f2, nspec)  
		Calculates the delay for the photons. Cloned from spectra.c

Arguments:		
	filename		File to output to
	PhotPtr p		Array of photons (pass global p)
	int nspec		Observer spectrum to dump for

Returns:
	To file only
  
Description:	
	Prep wipes clean the output file and writes the header, dump then
	builds up the output file in batches.

Notes:

History:
	14aug	-	Written by SWM
***********************************************************/
char delay_dump_file[LINELENGTH];
int delay_dump_bank_size=65535, delay_dump_bank_curr=0,delay_dump_spec;
PhotPtr delay_dump_bank;
int *delay_dump_bank_ex;

double
delay_to_observer(PhotPtr pp)
{
	plane_dummy observer;
	observer.x[0]   = pp->lmn[0]*geo.rmax;
	observer.x[1]   = pp->lmn[1]*geo.rmax;
	observer.x[2]   = pp->lmn[2]*geo.rmax;
	observer.lmn[0] = pp->lmn[0];
	observer.lmn[1] = pp->lmn[1];
	observer.lmn[2] = pp->lmn[2];
	return(pp->path + ds_to_plane(&observer,pp));
}

int 
delay_dump_prep (char filename[], int nspec, int restart_stat)
{
	FILE *fopen(), *fptr;
	char string[LINELENGTH];
	int i;
	
	//Allocate and zero dump files and set extract status
	delay_dump_spec = nspec;
	delay_dump_bank = (PhotPtr) calloc(sizeof(p_dummy), delay_dump_bank_size);
	delay_dump_bank_ex = (int*) calloc(sizeof(int),		delay_dump_bank_size);
	for(i=0;i<delay_dump_bank_size;i++)	delay_dump_bank_ex[i] = 0;


	//Get output filename
	strcpy(delay_dump_file,filename);
	if(restart_stat) return(0);

	/* If this isn't a continue run, prep the output file */
	if ((fptr = fopen(filename, "w")) == NULL)
	{
		Error("delay_dump: Unable to open %s for writing\n", filename);
		exit(0);
	}

	/* Check that nspec is reasonable */
	if (nspec < MSPEC)
	{
		Error("delay_dump: nspec %d below MSPEC value \n", nspec);
		exit(0);
	}

	/* 
	 * Construct and write a header string for the output file 
	 */
	fprintf(fptr, "# Python Version %s\n", VERSION);

	get_time(string);
	fprintf(fptr, "# Date	%s\n#  \n", string);
	fprintf(fptr, "# Spectrum: %s\n", xxspec[nspec].name);

	/* 
	 * Write the rest of the header for the spectrum file 
	 */
	fprintf(fptr, "# \n# Freq      Wavelength  Weight   "
			"  Last X       Last Y       Last Z     "
			"  Last L       Last M       Last N     "
			"Scatters RScatter Delay Extracted\n");	
	fclose(fptr);
	return(0);
}

int 
delay_dump_finish()
{
	if(delay_dump_bank_curr>0)
	{
		delay_dump(delay_dump_bank,delay_dump_bank_curr-1,delay_dump_spec, 1);
	}
	free(delay_dump_bank);
	free(delay_dump_bank_ex);
	return(0);
}

int 
delay_dump (PhotPtr p, int np, int nspec, int iExtracted)
{
	FILE *fopen(), *fptr;
	int nphot, mscat, mtopbot, i;
	double zangle;

	/* 
	 * Open a file for writing the spectrum 
	 */	
	if ((fptr = fopen(delay_dump_file, "a")) == NULL)
	{
		Error("delay_dump: Unable to reopen %s for writing\n", delay_dump_file);
		exit(0);
	}
	
	for (nphot = 0; nphot < np; nphot++)
	{
		if (p[nphot].istat == P_ESCAPE && p[nphot].nscat>0)
		{
			zangle = fabs(p[nphot].lmn[2]);
			/* 
			 * Complicated if statement to allow one to choose whether to construct the spectrum from all photons or just from photons which have scattered a specific number of times.  01apr13--ksl-Modified if statement to change behavior on negative numbers to say that a negative number for mscat implies that you accept any photon with |mscat| or more scatters 
			 */
			if (((mscat = xxspec[nspec].nscat) > 999 || p[nphot].nscat == mscat 
				|| (mscat < 0 && p[nphot].nscat >= (-mscat)))
				&& ((mtopbot = xxspec[nspec].top_bot) == 0 || (mtopbot * p[nphot].x[2]) > 0))
			{
				if (xxspec[nspec].mmin < zangle 
					&& zangle < xxspec[nspec].mmax)
				{	/* SWM 15/8/14 - Added path delay in comparison to photon heading straight from origin to rmax*/
					fprintf(fptr, "%10.5g %10.7f %10.5g %+10.5e %+10.5e %+10.5e %+10.5g %+10.5g %+10.5g %3d     %3d     %10.5g %i\n", 
						p[nphot].freq, C * 1e8 / p[nphot].freq, p[nphot].w, 
						p[nphot].x[0], p[nphot].x[1], p[nphot].x[2], 
						p[nphot].lmn[0], p[nphot].lmn[1], p[nphot].lmn[2], 
						p[nphot].nscat, p[nphot].nrscat, delay_to_observer(&p[nphot]),
						(iExtracted ? delay_dump_bank_ex[nphot] : 0));
				}
			}
		}
	}
	fclose(fptr);
	return (0);
}

int 
delay_dump_single (PhotPtr pp, int extract_phot)
{
	stuff_phot(pp, &delay_dump_bank[delay_dump_bank_curr]);			//Bank single photon in temp array
	delay_dump_bank_ex[delay_dump_bank_curr] = extract_phot;		//Record if it's extract photon
	if(delay_dump_bank_curr == delay_dump_bank_size-1)				//If temp array is full
	{																	//Dump to file, zero array position
		delay_dump(delay_dump_bank,delay_dump_bank_size,delay_dump_spec,1);
		delay_dump_bank_curr=0;
		int i;
		for(i=0;i<delay_dump_bank_size;i++) delay_dump_bank_ex[i]=0;	//Zero extract status of single array
	}
	else
	{
		delay_dump_bank_curr++;
	}
	return(0);
}
