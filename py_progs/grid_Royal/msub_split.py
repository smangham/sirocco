#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Description:  

Primary routines:

Notes:
									   
History:

090606 ksl Coding begun

'''

import sys
import os

def doit(filename='sscyg_sv.msub'):

	try:
		f=open(filename,'r')
		lines=f.readlines()
		f.close()
	except IOError :
		print "This file %s does not exist" % filename
		return []   
	
	i=0
	j=0

	file_no=1
	name=filename+'.%03d' % (file_no)
	print 'Starting file  ', name
	g=open(name,'w')
	os.system('chmod +x %s ' % name)

	while i<len(lines):
		line=lines[i].split()
		if len(line)>0 and line[0]=='mkdir':
			j=j+1
			if j%1000 == 0:
				file_no=file_no+1
				name=filename+'.%03d' % (file_no)
				print 'Starting file  ', name
				g.close()
				g=open(name,'w')
				os.system('chmod +x %s ' % name)
		g.write(lines[i])
		i=i+1
	g.close()

if __name__ == "__main__":
	import sys
	if len(sys.argv)>1:
		doit((sys.argv[1]))
	else:
		print 'usage: msub_split.py  filename'

