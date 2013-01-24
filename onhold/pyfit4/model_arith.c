


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:

  Description:	These are a generic set of routines to carry out 
  model arithmetic

  Arguments:		

  Returns:

  Notes:


  History:
0809	ksl	Began coding	
	

 ************************************************************************/
#include	<math.h>
#include	<stdio.h>
#include        <stdlib.h>
#include	<strings.h>
#include	<string.h>
#include	"log.h"
#include         "models.h"
#define    	BIG 1e32

#define 	LOCATE_ERR 	-99

CompPtr
get_component(model_name)
	char *model_name;
{
	int n;

	n=0;
	while (n<ncomps){
		if(strcmp(model_name,comp[n].var)==0){
			return(&comp[n]);
		}
		n++;
	}
	Error("get_component: There appears to be no component named %s\n",model_name);
	exit(0);


}
