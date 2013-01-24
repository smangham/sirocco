#!/usr/bin/env python 

'''
                    Space Telescope Science Institute


Synopsis:  
	make_royal.py [-t maxtime] [-a] root 
	
Description:  

	This is a program looks for the pf files in a directory
	and generates the commands to run a series of python models 
	on Royal.  

	-t sets the maximum time you want the program to run.  The input can
 		be in seconds or hours.  The routine decides which you have
		entered
	-a says to get all of the .pf files in the directory.  If not, then
		the routine searches for all of the .pf files that match
		root*.pf
	
	To start Royal processing your jobs just exuecute 

	root.msub



Arguments:  

	None

Returns:
	A new output file 

Notes:
									   
History:
08nov	ksl	Coded


'''

import os
import sys
import glob

def pfcheck(pflist):
	i=0
	pfgood=[]
	while i<len(pflist):
		try:
			f=open(pflist[i],'r')
		except IOError :
			print "This file does not exist:", pflist[i]
		        sys.exit(0)
		lines=f.readlines()
		j=0
		good=1
		while j<len(lines):
			z=lines[j]
			# print z
			zz=z.split()
			if len(zz)>0 and zz[0][0]!='#':
				if len(zz)<2:
					good=0
				if zz[1]=='$':
					good=0
			j=j+1
		if good==1:
			pfgood=pfgood+[pflist[i]]
		else:
			print 'Dropping: ',pflist[i]
		f.close()
		i=i+1
	return pfgood

def make_pbs(pffile,maxtime):

	root=get_root(pffile)
	pbsfile=root+'.pbs'

	# Put limits on the values for Royal
	imaxtime=1+maxtime/3600
	if imaxtime>8:
		imaxtime=8
	if imaxtime<1:
		imaxtime=1
	
	if maxtime> 3600* 7:
		maxtime=3600*7


	f=open(pbsfile,'w')
	f.write('#PBS -N %s\n' % root)
	f.write('#PBS -l walltime=%d:00:00\n'% (imaxtime))
	f.write('#PBS -l mem=2048mb\n')
	f.write('#PBS -M long@stsci.edu\n')
	f.write('#PBS -m a\n')    # Onlly send a message if the job aborts
#	f.write('#PBS -m bea\n')
	f.write('echo -n \"The job\'s unique ID: \"\n')
	f.write('echo $PBS_JOBID1 \n')
	f.write('cd  $PBS_O_WORKDIR\n')
	f.write('pwd\n')

	f.write('\npy -r -t %d %s\n' % (maxtime,root))
	f.write('restart.py -n 10 -C "msub %s.pbs"  %s\n\n' % (root,root))
	f.write('exit 0\n')


	f.close()
	return pbsfile



def make_msub(basename,pffiles):
	'''
	Create the initial msub file
	'''
	msubfile=basename+'.msub'
	msub=open(msubfile,'w')

	i=0
	while i<len(pffiles):
		pffile=pffiles[i]
		print 'Start', pffile
		root=get_root(pffile)
		print 'Root ',root

		working_dir='Dir_'+root

		msub.write('mkdir %s\n'%working_dir)
		msub.write('cd %s\n'%working_dir)
		msub.write('Setup_Py_Dir\n')
		msub.write('cp ../%s .\n' % pffile)
		msub.write('cp ../%s.pbs .\n' % root)
		msub.write('msub %s.pbs\n' % root)
		msub.write('cd ..\n\n')

		i=i+1


	msub.close()
	os.system('chmod +x '+msubfile)
	


def get_root(name):
	iend=name.rindex('.')
	return name[0:iend]

# Functions are located above here

argc=len(sys.argv)

maxtime=8*3600

if(argc==1):
	print "Usage: make_royal.py [-t max] [-a] root"
	print "-a all .pf files, -t tmax seconds or hours"
	sys.exit()

i=0
iall=0
while i< argc :
	if(sys.argv[i]=='-t'):
		maxtime=int(sys.argv[i+1])
		if maxtime<=8:
			maxtime = maxtime * 3600
	if(sys.argv[i]=='-a'):
		iall=1
	i=i+1

basename=sys.argv[argc-1]

pfnames=basename+'*pf'
if iall==1:
	pfnames='*.pf'

pflist=glob.glob(pfnames)

print pflist

#eliminate any parameter files that do not pass certain simple tests

pfgood=pfcheck(pflist)

print pfgood

make_msub(basename,pfgood)

i=0
while i<len(pfgood):
	pbsfile=make_pbs(pfgood[i],maxtime)
	i=i+1
