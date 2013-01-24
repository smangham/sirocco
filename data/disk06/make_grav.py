#!/usr/bin/env python 

'''
                    Space Telescope Science Institute


Synopsis:  
	This is a simple program to split disk06.ls into files
	with individual gravities that can be used by kfit

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

def read_file(filename,grav):
	'''
	Read the file, and return all the models with a given gravity.
	
	Note that this is realied tailor to this specific file, and that
	grav is a string
	'''
	try:
		f=open(filename,'r')
	except IOError :
		print "This file does not exist"
		sys.exit(0)
	
	# skip the first two lines
	f.readline()
	f.readline()
	models=[]
	
	xline=f.readline()
	while xline != '':
		print xline
		z=xline.split()
		if z[2]==grav:
			models=models+[[z[0],z[1]]]
		xline=f.readline()
	f.close()
	return models

def write_file(fileroot,grav,models):
	gg=float(grav)*10+0.01  # For roundoff
	gg=int(gg)
	name='%s.g%d.ls' % (fileroot,gg)
	g=open(name,'w')
	i=0
	while i<len(models):
		print '%30s %10s' % (models[i][0],models[i][1])
		g.write('%30s %10s\n' % (models[i][0],models[i][1]))
		i=i+1
	g.close()
	return

def do_one(root,grav):

	filein=root+".ls"
	models=read_file(filein,grav)
	if len(models)>0:
		write_file(root,grav,models)
	else:
		print 'Error: no models of grav %s found' % grav
	return


def do_many():
	'''
	This is just the master run file
	'''

	do_one('disk06','4.0')
	do_one('disk06','5.0')
	do_one('disk06','6.0')
	do_one('disk06','7.0')
	do_one('disk06','8.0')

	return
