                                                                                                                                 
Python version 54f  -- ksl
-------------

Note that Python_54d addressed some problems which had to do with getting
macro atoms and simple atoms to provide the same results.  54d changes were
made by ssim
                                                                                                                                 
Changed structure geo to add master atomic datafile and fixed concentrations
file.  In the case of fixed concentrations, the .pf file is changed so that
the name of the concentrations file is read in.  Storiing the atomic data file
in the geo structure will allow a change to the programs that read the
windsave file.  Note that the order of some of the input variables in the .pf
file is not optimal for the case where one wants to use an old windsave file.
The reading of the atomic datafile should be moved to after the decision to
read the old windsave file.  

Also modified Makefile so that Version number is determed from the Makefile
itself, and added the version number to various output files, including the
various ascii spectral files and the windsave file.

