SS Sep 2005

Notes on python56f (changes made since python56e):

1) For doing macro atom calculations with elements other than hydorgen
something needs to be done about metastable level populations (e.g. 
in the gound state of something like C II, we expect that the
population ratio of the j = 3/2 to j= 1/2 states should be 2:1 but
there's currently nothing in the code to drive this). Physically this
ratio is presumably forced by rapid collisional transitions between
the levels. So this approach is now adopted in the code. A level is
flagged as metastable by giving it a very large radiative lifetime in
the inout data (> 1.e15) - then when the macro pops routine sees this
level it gives it a false very fast collisional transition to the
ground state to force an LTE population ratio. Of course it would be
better to have the real atomic physics/data for the collisional
process but I think this will be a better approximation that ignoring
this sort of thing altogether (as was the case before).

To do this the following changes have been made and I've tested that
it all seems to work as it should:

- get_atomicdata.c: the radiative lifetime (rl) is now stored rather
                    than set to 0 for macro atom levels

- matom.c:          the macro_pops routine has been modified to
                    the fake fast collisional transitions (see the
                    matom.c file)
                

2) I'm a bit worried out the statistics of the macro atom ionization 
calculation and am experimenting with using approximate estimators
rather than the real ones. At the moment, in estimators.c, the MC
estimators are replaced by approximate ones for bf edges outside the
frequency range 7.5e12 to 1.2e16. The reason for this is that there
are essentially no MC packets at higher energies and so we keep
getting zeroes for the rates up there, which is a problem. Perhaps a 
more sophisticated approach cound be used where the number of packets
contributing is used as a cut to decide whether or not to keep the MC
estimator. But in the end the solution must be more packets when the 
ionizing continuua are important.

3) Bugs fixed

- spherical.c:  corrected lines relating to finding random point in 
		spherical coordinares

- matom.c:      corrected if statement clause on line 2590

