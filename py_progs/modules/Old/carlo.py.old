#!/usr/bin/env python 

'''
                    Space Telescope Science Institute


Synopsis:   carlo is a package of python routines intended to
	aid in the interpretation of fitting our python generated
	models 

History:

090116 ksl Coding begun

'''

import sys
import glob




def get_filenames(dirname='.',descriptor='*.pf'):
	'''
	get_filename(dirname,descriptor) locates and returns all 
	the filenames of the type specified by descriptor
	'''
	searchlist=dirname+'/'+descriptor
	names=glob.glob(searchlist)
	return names

def read_pf(filename='test'):
	'''
	read and parse a parameter file.  Note that
	this routine can handle the extension or not,
	but it will always try to read the .pf file
	'''

	# Check whether .pf has been added
	try:
		n=filename.rindex('.')
		name=filename[0:n]+'.pf'
	except ValueError:
		name=filename+'.pf'

	x=[]


	f=open(name,'r')
	lines=f.readlines()
	f.close()
	i=0
	while i<len(lines):
		q=lines[i].strip()
		if q[0]!='#' and len(q)>1:
			q=q.split()
			words=q[0].split('(')
			try:
				value=eval(q[1])
			except NameError:
				value=q[1]
		x=x+[[words[0],value]]
		i=i+1
	return x

def make_table(filelist):
	'''
	Make the man table
	'''
	table=[]
	i=0
	while i<len(filelist):
		print filelist[i]
		table=table+[read_pf(filelist[i])]
		i=i+1
	return table

def get_choices(table):
	'''
	Determine what variables are changed in the grid and what values
	there are of these
	'''

	unique=[]
	
	i=0
	while i<len(table[0]):
		j=0
		values=[]
		while j<len(table):
			record=table[j]
			values=values+[record[i][1]]
			j=j+1
		values.sort()
		v=[record[i][0],values[0]]
		j=1
		while j<len(values):
			if v[len(v)-1]!=values[j]:
				v=v+[values[j]]
			j=j+1
		if len(v)>2:
			unique=unique+[v]

		i=i+1
	
	return unique



		
		




