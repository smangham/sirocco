#!/usr/bin/python
import math

output=open('Do_modify_py','w')


mdisk_min=math.log10(1.e-08)
mdisk_max=math.log10(2.51e-8)
dmdisk=(mdisk_max-mdisk_min)/2.

print mdisk_min, mdisk_max, dmdisk

mwind_min=math.log10(1.58e-10)
mwind_max=math.log10(1.58e-9)
dmwind=(mwind_max-mwind_min)/2.
print mwind_min, mwind_max, dmwind

d=[8,32]
mexp=[0.0,1.0]
accel_exp=[6.]
vesc=[2.,3.]
accel_len=[5.e09,1.e10]

for i in range (0,3):
	mdisk=mdisk_min+i*dmdisk
	print mdisk
	mdisk=math.pow(10.,mdisk,)
	print mdisk

        for j in range (0,3):
		mwind=mwind_min+j*dmwind
		mwind=math.pow(10.,mwind)

		for k in range (0,len(d)):
			for l in range (0,len(mexp)):
				for m in range (0,len(accel_exp)):
					for n in range (0,len(vesc)) :
						for p in range (0,len(accel_len)):
							if (mwind < (0.11*mdisk)) :
								output.write("Modify_kwd sscyg_base sscyg_%02d_%02d_%02d_%02d_%02d_%02d_%02d  %.2e %.2e %2d %5.2f %5.2f %5.1f  %.2e \n" % (i,j,k,l,m,n,p,mdisk,mwind,d[k],mexp[l],accel_exp[m],vesc[n],accel_len[p]))

