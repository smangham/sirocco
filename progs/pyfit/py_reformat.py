#!/usr/bin/env python 

'''
                    Space Telescope Science Institute


Synopsis:  
	This is a simple program to reformat the models.ls file
	created to describe a set of gridded models into one
	which has all of the variables associated with each
	spectrum on one line so that this can be used for
	interpolating

Description:  

Arguments:  

	None

Returns:
	A new output file 

Notes:
									   
History:
07jan	ksl	Coded


'''

import sys

filein='sscyg_kgrid/models.ls'
fileou='sscyg_kgrid/models_reformat.ls'

try:
	f=open(filein,'r')
except IOError :
	print "This file %s does not exist" % fileine
	sys.exit(0)

g=open(fileou,'w')


# Read all of the header lines 
x=f.readline()

g.write('# Reformated with one spectrum per line. \n')
g.write('# Col. 2 gives the column with wavelengths\n')
g.write('# Col. 3 gives the column with flux\n')
g.write('# Last column has angle of spectrum   \n')

while (x != ''):
	z=x.split()
	if z[0][0]=='#':
		g.write(x)      
	else: 

		try:
			s=open(z[0],'r')
		except IOError :
			print "This file does not exist"
			sys.exit(0)

		line=s.readline()
		zz=line.split('Ang:')
		while len(zz)<2:
			line=s.readline()
			zz=line.split('Ang:')
		s.close()

		j=1
		while j<len(zz):
			print j,float(zz[j])
			g.write("%s 2 %3d "% (z[0],8+j))

			k=1
			while k<len(z):
				g.write(" %s " % z[k])
				k=k+1
				
			g.write("%8.0f\n"% float(zz[j]))
			j=j+1



	x=f.readline()
f.close()
g.close()


