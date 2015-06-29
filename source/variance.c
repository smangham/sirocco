#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
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
importance_map_iterate (
	WindPtr w)
{
	int n, iCells;
	double rFluxMax=0.0, rFluxMin=1e99, rFlux, rImpMin=1e99, rImpMax=0.0, rImp;

	if(geo.coord_type != SPHERICAL)
			iCells = NDIM*MDIM;
	else	iCells = NDIM;

	for (n = 0; n < iCells; n++)
    {
    	if(w[n].inwind == W_ALL_INWIND)
    	{
			if(plasmamain[w[n].nplasma].ntot > 0)
			{
				rFlux = plasmamain[w[n].nplasma].j;
				if(rFlux > rFluxMax) rFluxMax = rFlux;
				if(rFlux < rFluxMin) rFluxMin = rFlux;
				rImp = w[n].imp_inc;
				if(rImp > rImpMax) rImpMax = rImp;
				if(rImp < rImpMin) rImpMin = rImp;
			}
		}
	}

	if(rFluxMax <= 0 || rFluxMin >= 1e99)
	{
		Error ("importance_map_iterate: Could not find max/min photon flux!\n");
	}

	for (n = 0; n < iCells; n++)
    {
    	printf("IMPM  %d %g %g %g %g %d\n",geo.wcycle, w[n].x[0],w[n].x[1],w[n].x[2],w[n].imp_inc,w[n].inwind);
    	if(w[n].inwind == W_ALL_INWIND)
    	{
    		if(plasmamain[w[n].nplasma].ntot > 0) w[n].importance = 
    			pow(rImpMax / w[n].imp_inc, 1.0);
    		//if(plasmamain[w[n].nplasma].ntot > 0) w[n].importance = 
    		//	pow(rFluxMax / plasmamain[w[n].nplasma].j, 1.0); 
    	}
    	w[n].imp_inc = 0.0;
	}
	return(0);
}

int
importance_map_2d_flat (
	WindPtr w)
{
	int n;

	for (n = 0; n < NDIM*MDIM; n++)
    {
   		w[n].importance = 1.0;
	}
	return(0);
}

int
importance_map_sphere (
	WindPtr w)
{
	int n, iShell;
	iShell = 0;

	for (n = 0; n < NDIM; n++)
    {
    	if(w[n].r < geo.vr_shell_radius[iShell])
    	{
    		w[n].importance = geo.vr_shell_importance[iShell];
    	} 
    	else
    	{
    		iShell++;
    		w[n].importance = geo.vr_shell_importance[iShell];
    	}
	}
	return(0);
}

int 
setup_importance_map (
	WindPtr w)
{
	if(geo.coord_type == SPHERICAL && geo.vr_sphere_shells > 0)
	{
		importance_map_sphere(w);
	}
	else
	{
		importance_map_2d_flat(w);
	}
	geo.vr_np = 0;
	return(0);
}

int dump_imap_info(WindPtr w)
{
	int n, iCells;

	if(geo.coord_type != SPHERICAL)
			iCells = NDIM*MDIM;
	else	iCells = NDIM;

	for(n = 0; n < iCells; n++)
	{
    	if(w[n].inwind == W_ALL_INWIND || w[n].inwind == W_PART_INWIND)
    	{
    	 	printf("WIND: %d %g %g %g %d %g %g %g\n",geo.wcycle-1, w[n].xcen[0], w[n].xcen[2], w[n].importance, 
    				plasmamain[w[n].nplasma].ntot, plasmamain[w[n].nplasma].j, plasmamain[w[n].nplasma].ip, w[n].imp_inc);
		}
	}
	return(0);
}