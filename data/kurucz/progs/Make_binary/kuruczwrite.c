
/*-************************************************************************
                    Space Telescope Science Institute


  Synopsis:
		make_binary  file  [yes}

		where file is the complete filename of the binary file, eg ...dat
		and the optional yes invites the program to make an additional
		binary file which only contains the wavelengths.

  Description:	This program reads ascii files from which have been
	put on disk from the kurucz tape received in Jan. 91 and writes 
	out binary files

  Arguments:		

  Returns:

  Notes:
		The program allows for the fact that the various ascii files
		are of different lengths.  The output files put zero in for
		the missing values.

  History:
	
	2/1/91 	Coded ksl

 ************************************************************************/

#include <stdio.h>
#include <sys/file.h>

#define	ARRAYSIZ	1500
#define NUMBER_OF_VALS 	1500

#define	ASCIIDIR        "/home/long/kurucz/ascii/"  
#define BINARYDIR	"/home/long/kurucz/data.91/"
#define ARRAYSIZE 	8 * 1500		/* No. of characters in each output array */
#define LINELEN		100
#define NOK		-1
#define OK		0



main(argc,argv)
int argc;
char *argv[];
{
int i,j,k;
double wave[ARRAYSIZ],flux[ARRAYSIZ];
double w[ARRAYSIZ],f[ARRAYSIZ];
double ww[ARRAYSIZ],ff[ARRAYSIZ];
int no[ARRAYSIZ],n[ARRAYSIZ],nn[ARRAYSIZ];
int npts,nspec;
char	wavefile[LINELEN];


if ( argc < 2 ) {	
	fprintf(stderr,"usage: kuruczlist model#\n [yesno]");
	exit(1);
	}


/*Read a file with all the wavelengths */
npts=raw_kurucz("fp00t50000g50k2.dat",no,wave,flux);
printf("No of points: %d\n",npts);


printf("%d %s\n",argc,argv[2]);

if(argc==3&&strncmp(argv[2],"yes",1)==0){ /*Then make a new wave file */
		writekurucz("kurucz_wavelengths.dat",wave,wave);
		printf("Created a new wavelength file: kurucz_wavelengths");
		}

/*Read the file you want to convert to binary */
nspec=raw_kurucz(argv[1],n,w,f);

/*Initiialize the output arrays */
for(i=0;i<ARRAYSIZ;i++){
	if(i<npts) {
		ww[i]=wave[i];
		nn[i]=no[i];
		}
	ff[i]=0.0;
	}

/* Now put the fluxes in the appropriate place in the
output arrays */

for(i=0;i<nspec;i++){
	j=0;
	while(nn[j]!=n[i])j++;
	ff[j]=f[i];
	}


writekurucz(argv[1],ww,ff);

}


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:
		raw_kurucz(filename,n,wave,flux)

  Description:	This program reads the ascii version of the new
	kurucz files received by ksl (jan 91) and stores the results
	in three arrays n, wave and flux.

  Arguments:		

  Returns:

  Notes:

  This is meant to read the .dat files which came on the tape labelled
  Kurucz P00K2.  Not all the models on this tape have the same
  wavelengths; in particular the red stars do not have the most
  ultraviolet wavelengths.  The wavelength numbers, i.e. the second
  column, in the files seem to always correspond to the same wavelengths,
  however.


  History:

  Coded 1/31/91 by ksl
	

 ************************************************************************/

#include <stdio.h>
#include <strings.h>

#define NSKIP	3		/*No of comment lines at beginning*/
#define	NOK		-1
#define WAVEMAX	10000. /*Maximum wavelength to return */

raw_kurucz(filename,n,wave,flux)
char *filename;
int	n[];
double wave[],flux[];

{

FILE *fopen(),*ptr;
char *fgets();

char name[LINELEN];
char line[100];
int i;

strcpy(name,ASCIIDIR);
strcat(name,filename);
if((ptr=fopen(name,"r"))==0){
	perror("raw_kurucz");
	exit(NOK);
	};

for(i=0;i<NSKIP;i++)fgets(line,100,ptr);
i=0;
while( fgets(line,100,ptr) != NULL ){
	sscanf(line,"%*s %4d %lf %*f %lf",&n[i],&wave[i],&flux[i]);
	wave[i]*=10;
	i++;
	if(wave[i-1]>WAVEMAX) return(i);
	} 

return(i);
}


writekurucz(model,wavelength,flux)
char *model;
double *wavelength;
double *flux;
{
  	int fd;
	int i;
	char fileout[LINELEN],string[LINELEN];
	
	/*Clip the .dat off the filename*/
	strcpy(fileout,BINARYDIR);
	i=strlen(model);
	strncat(fileout,model,i-4);
	printf("fileout %s\n",fileout);
	if ( (fd = open(fileout,O_CREAT|O_WRONLY)) == -1 )
	{	perror("kurucz on output");
		return(NOK);
	}
	if ( write(fd,(char *)flux,ARRAYSIZE) != ARRAYSIZE )
	{	fprintf(stderr,"kurucz: to few flux values\n");
		return(NOK);
	}
	close(fd);

	return(OK);
}
