Python version 55b - SS Feb 05
------------------------------

At time of writing, it is my plan that this will be the version of Python
used for the initial YSO paper on H I recombination lines. Various minor
bugs have been fixed since Python 54 (many relating to correct application
of the "extract" process in macro atom calculations). One important change
is to the He macro atom lines data (which had incorrect energy levels 
before).

This version has been tested in various ways:

(1) "knigge7.pf" models (i.e. something like Figure 7 in the original paper)
have been done to compare results for "simple" calculations with the new
code to those from python_46. Excellent agreement was found in this 
comparison when small changes in grad_v and the range of frequencies used
by radiation.c for comparison with the older version.

(2) H+He macro atom models for CV winds have been used to test the extract
process (vs live.or.die) for all the scattering modes. With the exception of
the He II Lyman edge, excellent agreement has been found between extract and
live.or.die. Comparing simple and macro atom calculations suggest that the
biggest differences are in the depths of the bf edges. This probably deserved
closer investigation at some point soon.

(3) Two-level atoms for the formation of H I and He II Lyman alpha have been 
made and used to test the macro atom and simple formulations (which should
agree for this set of atomic data). Agreement is not perfect but I suspect that
this is, at least in part, since the level populations are computed 
differently in the two cases. Nothing seems very wrong but again more 
careful investigation would be a good idea sometime. H Lyman alpha is
closer to perfect agreement than He Lyman alpha.

(4) Complete atomic data (H+He macro, everything else simple) calculations 
have been run for a CV model to see what they look like. The Lyman edge does
some interesting things when macro atoms are used.


SS - Mar 05
-----------

One more modification made for the YSO calcualtions - if geo.disk_illum is
set to 1 then the disk becomes transparent - this is intended as a crude
model for a highly reflecting disk.
