#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************
             University of Southampton

Synopsis: 
  get_disk_params sets up the disk parameters according to user inputs, 
  e.g. the temperature profile, accretion rate etc.
   
Arguments:		

Returns:
  disk_illum - this is used by python.c and so needs to be returned
 
Description:	

Notes:

History:
	1502  JM 	Moved here from main()
	1510	ksl	Modified to restore illumination
			options, which were brokedn
    1712    ksl Collected input parameters for disk into
                a single routine and modified some of the
                input fil names to be more uniform

**************************************************************/


double
get_disk_params ()
{

  rdpar_comment ("Parameters for the Disk (if there is one)");
  rdint ("Disk.type(0=no.disk,1=standard.flat.disk,2=vertically.extended.disk)", &geo.disk_type);
  if (geo.disk_type != DISK_NONE)
  {
    rdint ("Disk.radiation(y=1)", &geo.disk_radiation);
  }
  else
  {
    geo.disk_radiation = 0;
  }
  get_spectype (geo.disk_radiation, "Disk.rad_type(0=bb,1=models)_to_make_wind", &geo.disk_ion_spectype);


  if (geo.disk_type == DISK_NONE)
  {
    return (0);
  }
  rdint ("Disk.temperature.profile(0=Shakura-Sunyaev;1=readin,2=yso)", &geo.disk_tprofile);
  if (geo.disk_tprofile == DISK_TPROFILE_STANDARD)
  {
    geo.disk_mdot /= (MSOL / YR);       // Convert to msol/yr to simplify input
    rddoub ("Disk.mdot(msol/yr)", &geo.disk_mdot);
    geo.disk_mdot *= (MSOL / YR);
  }
  else if (geo.disk_tprofile == DISK_TPROFILE_READIN)
  {
    rdstr ("Disk.T_profile_file", files.tprofile);
    geo.disk_mdot = 0;
  }
  else
  {
    geo.disk_mdot = 0;
  }

  /* Set a default for diskrad for an AGN */
  if (geo.system_type == SYSTEM_TYPE_AGN)
  {
    geo.diskrad = 100. * geo.r_agn;
  }

  rddoub ("Disk.radmax(cm)", &geo.diskrad);
  Log ("geo.diskrad  %e\n", geo.diskrad);

  geo.diskrad_sq = geo.diskrad * geo.diskrad;

/* If diskrad <= geo.rstar set geo.disk_type = DISK_NONE to make any disk transparent anyway. */

  if (geo.diskrad < geo.rstar)
  {
    Log ("Disk radius is less than star radius, so assuming no disk)\n");
    geo.disk_type = DISK_NONE;
  }

  if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
  {                             /* Get the additional variables need to describe a vertically extended disk */
    rddoub ("Disk.z0(fractional.height.at.diskrad)", &geo.disk_z0);
    rddoub ("Disk.z1(powerlaw.index)", &geo.disk_z1);
  }
  return (0);
}
