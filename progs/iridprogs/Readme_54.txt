Changes in python_54

Python version 54
-----------------

This integrates changes made to version 53 to deal with coordinate systems
(KSL - see version 53b readme) and several minor changes made by SS 
relating to use of the analytic treatment of disk re-radiation (which 
was not treated properly before) and also calculations involving a 
macro atom treatment of hydrogen.

The joint H+He macro atoms have been tested using a CV model and appear to
work (the He recombination lines are reasonably strong). 

BUG: The anisotropic scattering model 1 (i.e. NOT the "thermal trapping model")
doesn't seem to work properly now (something in "reweightwind" is going wrong).
The "thermal trapping" model seems to work, however, and that's the one I'm 
using now.
