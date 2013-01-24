

/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:

  Description:	This program 

  Arguments:		

  Returns:

  Notes:

  History:
	

 ************************************************************************/
#include	<math.h>
#include	<stdio.h>
#include	<strings.h>
#define 	LINELENGTH 	132
#define         RAD         57.29578

main (argc, argv)
     int argc;
     char *argv[];
{
double x[3],y[3], z[3];
double phi,theta,psi; // Goldstein p 107
struct basis {
   double a[3][3];
   } base;

phi=30.;
theta=45;

phi/=RAD;
theta/=RAD;

x[0]=cos(phi)*cos(theta);
x[1]=sin(phi)*cos(theta);
x[2]=sin(theta);

y[0]=0;
y[1]=0;
y[2]=1;

create_basis(x,y,&base);

project_from(&base,x,z);

printf("z %f %f %f\n",z[0],z[1],z[2]);

project_to(&base,x,z);

printf("z %f %f %f\n",z[0],z[1],z[2]);

project_to(&base,y,z);

printf("z %f %f %f\n",z[0],z[1],z[2]);

}
