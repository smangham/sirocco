#!/usr/bin/env python 

'''
                    Space Telescope Science Institute


Synopsis:  
	This is a program that is designed to expand a .pf file into a grid of .pf files.

Description:  

	This program was written for creating grids of models with my radiative transfer
	program python, although it should work with any .pf file in the case where one
	wants to run a grid of models.  As an example this should be a reasonable way
	to do much of the setup work for runing a bunch of ISM calculations.

	To use the program take a .pf file which is correctly formatted to run with python.
	This is going to be the basis of all of the .pf files.  The grid is going to be
	created for the variables in the .pf file which values of $.  

	For example, if your .pf file reads
		variable1  42
		variable2  $
		variable3  55
	
	Run the program. It will ask for the .pf file, herafter grid.pf, to use, 
	and a basename for all of the output files.  It will read through grid.pf and 
	will ask the user to give the values of each of the variables to grid.

	In the example above you will be asked to provide values for variable2 only

	The format of the inputs one must provide are  unusual! (They are read with the input command in
	python, the scripting language'.  The two possibities that are supported by the 
	program are:
	
		a list format which means one must enclose everything in [], e.g.
			[14,27,32]
		which would say that for this variable you want values of 14, 27, and 32.
	
		a special format to create a logarithminc grid of variables, e.g.
			log_ints(xmin,xmax,n)
		where xmin and xmax are the minimum and maximum values of the variable in 
		question and n is the number of variables.  Note that you enter the
		actual values of xmin and xmax and not the log of this.  For example:
			log_ints(1.e16.,1.e18.,3)
		will end up produciing a list with variables [1.e16,1.e17,1.e18]

		another special format of a similar form to create linear intervals
			lin_ints(1.0,2.0,10)

	For the exampe above, if you gave a list format [14,27,32], three pf files
	will be produced
Arguments:  

	None

Returns:
	The program produces a lot of files:

	a file which lists what was varied and the values for them.  The name of this file
	is the same as basename, with a .ls extensition attached, e.g vwhyi.ls
	This is intended for use with a fitting program

		# lum_bl(ergs/s)
		# t_x
		PYGRID/vwhyi_0000.spec         1         3 
		PYGRID/vwhyi_0001.spec         1         4 
	
	a file which can run python for the grid, one model after another.  This has the prefix
	Run_ and the remainder is the basename, e.g something like, Run_vwhyi.  It begins 
	something like
		py vwhyi_0000
		py vwhyi_0001
	A number of .pf files, that are like the original pf files, but now the $ have been
	replaced with the variables associated with the grid.


Notes:
	If the program finds .pf files that seems to have been created earlier, it will ask 
	you if you want to delete them

Bugs:
	The program requires one to enter a lot of inputs by hand that are complicateed.  
	This is not ideal, but the program has to be run so few times that changing this
	seems too much work.  One could fairly easily have the program look in another
	file for the variables.
									   
History:
0709	ksl	Coded and debugged


'''

import sys
import math
import os
import fnmatch
from numpy import *


# Next routine is fro kslio, just to make the installation slightly easier 
def get_input(question,answer):
	""" normal usage answer=get_input(question,answer) """
	string=question+' ('+answer+'):  '
	ss=raw_input(string)
	if ss=='':
		ss=answer
	return ss
#End of code from kslio

def log_ints(xmin,xmax,nints):
	'''
	Create logarithmic intervals between xmin and xmax
	'''
	if xmin<0.0 or xmax <= 0.0:
		print 'Error - Can not take the log of a negative number  %f %f ' % (xmin,xmax)
		return array([xmin,xmax])
	lmin=math.log10(xmin)
	lmax=math.log10(xmax)
	dx=(lmax-lmin)/(nints-1)
	x=[]
	i=0
	while i < nints:
		q=lmin+dx*i
		q=pow(10.,q)
		x=x+[q]
		i=i+1
	return array(x)


def lin_ints(xmin,xmax,nints):
	'''
	Create linear intervals between xmin and xmax
	'''
	lmin=float(xmin)
	lmax=float(xmax)
	dx=(lmax-lmin)/(nints-1)
	x=[]
	i=0
	while i < nints:
		q=lmin+dx*i
		x=x+[q]
		i=i+1
	return array(x)

def expand_array(arrayin,newvar,names_old):
	'''
	Startign with a new variable and the possible 
	values of that variable, calculate all the
	possibilities of the old variables contained
	in array in, and also generate the names
	for each model.

	The routine returns the output array of
	all of the choices and the names
	'''
	names=[]
	if len(arrayin)==0:
		# Then we are just beginning
		outarray=newvar.reshape(len(newvar),1)
		nrows=len(newvar)
		n=0
		while n<nrows:
			names=names+["%02d"%(n)]
			n=n+1
	else:
		inshape=arrayin.shape
		rows_in=inshape[0]
		nmax=newvar.shape[0]
		rows_out=rows_in*nmax
		cols_in=inshape[1] 
		cols_out=cols_in+1 
		outarray=zeros([rows_out,cols_out])
		row=0
		nold=0
		print 'Creating %d output rows' % rows_out,'for ',newvar
		while row<rows_out:
			nold=row/nmax
			fname=names_old[nold]
			fname="%s%02d"% (fname,row%nmax)
			# names=names+[fname]
			names.append(fname)
			col=0
			while col < cols_out:
				if col < cols_in:
					outarray[row,col]=arrayin[nold, col]
				else:
					outarray[row,col]=newvar[row % nmax]
				col=col+1
			if row % 1000 == 0:
				print 'Finished row %d ' % row
			row=row+1
	return outarray, names


# Here names and var_array are parallel arrays.
def export_results(lines, basename,names,var_array):
	'''
	This routine writes out the .ls file

	lines		lines are the variable names (used to create the header).
	basename	The base names
	names		The extension for all of the names
	var_arrya	A numpy.array containing all the values of the parameters that vary

	Note	At this point the variables are somewhat disconnected from the names
		of the variables, so it is too late to do any kind of culling
		of the models

	090531	ksl	This routine or the part the writes the heeader needs
			to be synchronized with define_model (in grail.make_mod)

	'''

	# Write the header
	define_model(basename,lines,opt='new')

	outfile=basename+".ls"
	gout=open(outfile,"a")

	# Write header

	# n=0;
	# while n<len(lines):
	# 	z=lines[n].split()
	# 	gout.write("# VAR %s\n" % z[0])
	# 	n=n+1

	# Write out all of the rows of individual models

	print 'var_array',var_array
	nrows=var_array.shape[0]
	ncols=var_array.shape[1]
	n=0
	while n<nrows:
		filename="PYGRID/%s_%s.spec" % (basename,names[n])
		gout.write("%s " % filename)
		m=0
		while m<ncols:


# Conditions to eliminate specific models ought to go here

# Finish the conditions
			gout.write("%9.3g "% var_array[n,m])
			m=m+1
		gout.write("\n")
		n=n+1
	gout.close()
	print 'The summary file is %s' % outfile


def define_model(modelname,modelvariables,opt='old'):
	'''
	Create  a ls file to store the model results in along with the model variables

	Parameters:
		modelname	the rootname of the model file, e.g test
		modelvariables  A list, with each element of the list defining
				how something should vary, e.g 


	[['alpha',0,4,10,'log']] where alphs is the name of the variable, 0 and 4 are
	lower an upper limits for the variable and 10 number of intervals for delta 
	in the grid.
		
		opt		A word indicating whether one shuld overwrite
				the file, 'old' the default, will try to use 
				the file that exists.

	Usage:define_model('power',[['alpha',0,4,1]],opt='new')


	090527	ksl	Eliminated unused func variable and modified to allow for log and
			linear variables
	'''

	# test if the model already exists and if so check that we have both the grid file and all of the 
	# files that are supposed to exist.  If one is missing attempt to recreate it.

	if os.path.exists(modelname+'.ls') and opt != 'new':
		check_model_list(modelname,func)
		print 'The generic model %s exists already.  Header not updated' % modelname
		return


	# if the model does not exist at all then create the header to the ls file and return

	g=open(modelname+'.ls','w')
	i=0
	while i<len(modelvariables):
		var=modelvariables[i]
		try:
			string='# VAR  %20s   %8.2g    %8.2g   %4d' % (var[0],var[1],var[2],var[3])
		except TypeError:
			print 'Error making VAR string: ',var
		type='linear'
		if len(var)>4:
			type=var[4]

		string=string+' %-8s ' % type
		print string
		g.write('%s\n' % string)
		i=i+1

	g.close()
	print 'define_model: created a new model file: '+modelname+'.ls' 
	return

def create_parameter_files(basename,names,lines,z):
	'''
	names are the number portions of the names
	lines contains the original parmeter fiels
	z contains the values
	'''
	n=0

	# print 'foo ',z

	while n<len(names):
		pffilename=basename+'_'+names[n]+'.pf'
		pffile=open(pffilename,'w')
		m=0
		k=0
		while m<len(lines):
			words=lines[m].split()
			if len(words)>1 and words[1]=='$':
				pffile.write('%s		%8.3g\n' % (words[0],z[n,k]))
				k=k+1
			else :
				pffile.write(lines[m])
			m=m+1
		pffile.close()
		n=n+1


def create_runfile(basename,names):
	'''
	Create a file called Run_basename which can be used to
	run the grid
	'''
	n=0
	runfile='Run_'+basename
	rfile=open(runfile,'w')
	rfile.write("#!/usr/bin/env bash\n")
	while n<len(names):
		pffilename=basename+'_'+names[n]
		rfile.write('py %s\n' % pffilename)
		n=n+1
	rfile.close()
	os.system('chmod +x '+runfile)
	print 'Create a file %s that can be used to run the models one at a time'


def cleanse(basename):
	files=os.listdir('.')
	pattern=basename+'_*.pf'
	pffiles=fnmatch.filter(files,pattern)
	if len(pffiles)> 0:
		print "There are .pf files with this basename."
		print pffiles
		i=int(get_input("Delete them -1, Quit 0, Continue 1","1"))
		if i==1:
			return
		if i==0:
			sys.exit()
		# Delete them
		os.system("rm "+pattern)





def parse_variable_pf_line(line):
	'''
	This is a routine to read a pf line that contains limits on what
	the the variable can be

	The record should have the following format 

	name(extra_stuff)     $ or value  #  min=whatever max=whatever  or n=whatever log or linear

	The portion behind the # is treated as a kind of grammar

	The record that is returned should include the following

	name, current_value, min, max, log or linear


	090531	ksl	Routine was taken from pf.py as part of effort ot homogenize the way 
			this routine and the grail routines operate.  Depending on future
			it's likely these routines should be integrated together
	00531	ksl	This version has been modified to make it directly compatable with 
			expand_var.  The return line is 

			name,min,max,ndetal,log or linear
	'''

	# set up some defaults for the record
	name='unknown'
	value=0
	min=-1e32
	max=1e32
	type='linear'
	ndelta=10

	unknown=[name,value,min,max,ndelta,type]

	q=line.strip()
	if len(q)==0 or q[0]=='#':
		return unknown
	q=q.split()
	name=q[0].split('(')
	name=name[0]  # Really just want the string here

	i=0
	process='no'
	while i<len(q):
		if process=='yes':
			if q[i]=='log':
				type='log'
			elif q[i]=='linear':
				type='linear'
			else:
				pair=q[i].split('=')
				if pair[0]=='min':
					min=eval(pair[1])
				elif pair[0]=='max':
					max=eval(pair[1])
				elif pair[0]=='n':
					ndelta=int(pair[1])
				else:
					print 'Error: parse_variable_pf_line: Unknown choice: %s  Assuming linear' % (pair)
		if q[i]=='#':
			process='yes'  # Start processing after comment
		i=i+1

	if process=='no':
		print 'Error: process_variable_pf:  Could not process: ',line
		raise IOError
		return []

	# Now deal with the initial value of the variable
	# try:
	# 	value=eval(q[1])
	# except SyntaxError:
	# 	if type=='log' and min > 0 and max >0:
	# 			print 'logarithmic'
	# 			xmax=math.log10(max)
	# 			xmin=math.log10(min)
	# 			value=random.uniform(xmin,xmax)
	# 			value=math.pow(10.,value)
	# 	else:
	# 		value=random.uniform(min,max)
	
	# return [name,value,min,max,ndelta,type]
	return [name,min,max,ndelta,type]



def expand_one_var(one_var):
	'''
	Return the allowed values for one variable from that variable_list
	where one_var contains
		[whatever,0,1,10,'linear']

	090531	ksl	Stolen from make_mod in grail.  As with the routine
			above it is possible that gen_grid should be integrated
			into the grail package
	'''
	print 'expand_one_var',one_var
	if one_var[4]=='log':
		vmin=math.log10(one_var[1])
		vmax=math.log10(one_var[2])
	else:
		vmin=one_var[1]
		vmax=one_var[2]

	ndelta=one_var[3]
	delt=(vmax-vmin)/ndelta

	values=[]
	i=0
	while i<=ndelta:
		x=vmin+i*delt
		if one_var[4]=='log':
			x=math.pow(10,x)
		values=values+[x]
		i=i+1
	return values



# Now start the main body of the program

def gen_grid(pfname='grid.pf',basename='vwhyi'):
	'''
	Main routine which generates the grid
	'''
	# input_file=get_input("Starting.parameter.file","grid.pf")
	# basename=get_input("basename for grid","vwhyi")

	cleanse(basename)

	# Read the input parameter file    
	try:
		file_in=open(input_file,"r") 
	except IOError:
		print 'Error: Could not open %s' % input_file
		sys.exit(0)


	lines=file_in.readlines()
	file_in.close()

	print "Enter values for each variable, in bracketed format, e.g. [1, 4, 7,8]"

	# The sum file lists all of the values used for each variable
	sumfile=basename+'.sum'
	sumf=open(sumfile,'w')

	i=0

	names=[]
	olines=[]
	z=[]
	while i<len(lines):
		line=lines[i]
		print 'line: ',line
		words=line.split()
		print words
		if len(words)>1 and words[1]=='$':

			# now check try to parse the line looking at the words after the $
			try:
				xvar=parse_variable_pf_line(line)
				olines=olines+[xvar]
				zz=expand_one_var(xvar)
				zz=array(zz)

			except IOError:
				# try (the old interactive approach)

				print 'Warning: This is an old style input approach and my not be what you want!'

				# Next line is a kluge to allow oldstyle line.  Altering log_ints and lin_ints
				# could ameiliorate this 

				olines=olines+[[words[0],0.0,0.0,0,'unknown']]
				string='%s [values]? ' % words[0]

				# This is the line that generates the values for this line
				# It executes the command string and puts the results in
				# zz
				zz=input(string)
				zz=array(zz)

			# Write the values for this parameter to the sumf file
			sumf.write('%s' %words[0])
			kk=0;
			while kk<len(zz):
				sumf.write('  %8.3g'%zz[kk])
				kk=kk+1
			sumf.write('\n')

			# Next line performs the "multiplication" to keep track of all the possibilties
			# Note z is the old set of choices, and zz is the set one wants to multiply in
			z,names=expand_array(z,zz,names)
		# print 'z xxx',z
		i=i+1

	sumf.close()

	# Note that names here refers to the names of the output pf files

	# This creates a summary of all of the variables that are changed.  It is to be used in
	# fitting programs.  (This is the .ls file)

	# olines	olines are the variable names (used to create the header.
	# basename      The base names
	# names         The extension for all of the names
	# z		A numpy array containing all of the values for the variablesvariables
	print 'IMPORTANT: Creating %d pf files' % len(names)

	export_results(olines,basename,names,z)

	# This creates the individual .pf files
	create_parameter_files(basename,names,lines,z)

	# This creates the file that will run the grid of models
	create_runfile(basename,names)


if __name__ == "__main__":
        import sys
	if len(sys.argv)==1:
		input_file=get_input("Starting.parameter.file","grid.pf")
		basename=get_input("basename for grid","vwhyi")
	elif len(sys.argv)==3:
		input_file=sys.argv[1]
		basename=sys.argv[2]
	else:
                print 'usage: gen_grid.py  filename'
		sys.exit()

	gen_grid(input_file,basename)
