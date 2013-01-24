Notes on python_54b -- ksl

As an initial step to reduce the numter of files in the directory, I
consolidated all of the Readme...txt files from version 46a to version
53 into a single chronological file Readme_46_53.txt.  Old files are
in the directory Readme.hold for now.  I kept the original Readme file,
now Readme.txt that described some of my (ksl's) efforts prior to SS's
entry on the project.

Made minor changes to python.c -- just for better documentation.


Made minor changes to kbf_need and scatter and doppler in resonate.c -- 
	to eliminate some unnecessary steps and add documentation

Made minor changes to trans_phot -- essentially documentation

matom.c -- Deleted checks in photo_gen_kpkt and photo_gen_matom that
	photon is generated inside the disk, and replanced with message
	that will kill the program if this occurs.  it should not be
	happening at all as this is checked previously, and we should not
	be doing it twice.  Checks in matom should eventually be deleted
	as superfluous.


Note -- The changes described here were made in the process of trying to 
understand what is actually happening in the macro atom approach.  Both
Stuart and I were seeing differences in the way the program was working
with macro atoms and with simple atoms that were unanticipated.:


Note: There is a significant problem with trying to run old wind models, in that
a lot of variables that should not be recalculated are being recalculated.
It's tricky because there are some things one might like to change, but many
that should not.  This greatly complicates analyzing problems in detailed
spectra.



