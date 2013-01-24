This file contains the the notes that were writtind on changes that were made to Python 
in 2004 from version python_46a up to the creation of python_54

Changes in python46.a relative to python46

04Mar -- ksl

wind_updates2d.c -- Added code that calls a routine that prints out info
	about how t_r is evolving in the diagnostic file.  The code really
	was adapted from a similar routine used by py_wind.  

define_wind.c -- Attempted to understand a problem in which python
	was generating div v's which were less than zero.  The problem 
	is limited to regions outside the windcone.  I made minor modifications
        to change the way div*v was calculated, and then put in error
	checking that will complain more loudly if div v is less than
	zero inside the wind.

templates.h, wind2d.c, pyray.c, python.c, py_wind.c, windrtheta.c, 
	Changed calls to various routines and eliminated code already commented
	out having to do with a dignostic file called windsumfile.  

bb.c -- Modified the generation of frequencies for photons in the circumstance
	that one wants photons much higher than h nu/ kT and when one wants
	photons only in the range h nu/kT less than 0.2.  

templates.h, recipes.c (qromb, others), pdf.c (pdf_gen_from_func),
anisowind.c, random.c, roche.c --

	I was getting warnings from -Wall on some of the function
	calls.  I believe this is because, e.g. trapzd, the function pointer
	definitions are not really complete.
	     double (*func) ();         /* ANSI: double (*func)(double); */

	See also
		http://www.eskimo.com/~scs/cclass/int/sx10c.html
	I was a bit surprised that I had not seen these warnings before, but
	possibly the compiler has improved in this regard.  

	I did not completely fix all of the errors.  The one associated with
	rtsave looked too hard for me right now


pdf.c -- I found a problem in pdf_gen_from_array in the situation where one
	wanted the cumulative distribution function to have a break at a 
	specific point.  The problem resulted in a non-monotonic cdf, or
	more properly a cdf in which the array elements were out of order.
	Recall that y is the integral of the probability density from the
	minium value to x, and what I actually try to do is calculate the
	values of x_i so that y_i runs at equally spaced intervals from 0 to
	1.  But I also allow for calculating the y_i as specific points (so
	that extremes of the distribution are not completely unsampled). 
	In certain instances the interaction of calculating points at
	specific locations in y and equal sampling in x resulted in a pdf
	that was out of order.  I added some steps to attempt to check whether
	it was out of order, and reorder the cdf if that was the case.  The
	problem seems to be solved.
	The problem is associated with the fact that I attempt to produrce
	a cdf in which one actually attempts to calculate the points

Discussion of problems, not yet resolved.

photon_gen.c, python.h 
	There are several small problems in the way the radial positions of
	photons radiated by the disk are bein generated.  One problem was
	that I was using stratified sampling of the disk annuli, but this can
	result in an odd photon distribution a small number of photons are
	repeatedly generated.  This particular problem is fixed.  

	

Prior to banding one set up a single set of annuli and used them trhoughout
the ionization calculation part of the process.

At present the way photons are generated in the disk is to pick a wavelength
range and regrid the disk such that an equal number of photons is expected
from each disk annulus.  We regrid whenever we chance the frequency range
which is often, as a result of what I sometimes call banding that assures
a certain number of photon bundles in a given frequency range are generated.
In principle this is a good thing, since it adaptively adjusts the boundaries
where most of the photons are coming out.

The result of banding is that variables that parts of the disk structure 
that used to record how many photons were absorbed by the disk are no
longer useful.  As a result, one cannot handle heating of the disk 
by photons in any sensible way.

It is not entirely obvious that the regriding of the disk is not just a bit
too clever.  It has the property that it can cause pretty large gridding
at the outer edge of the disk, and thus the few photons that are out there
are not necessarily produced with the correct temperature distribution.  
It is tempting to change to a scheme where the boundaries are fixed at the beginning,
and simply increase NRINGS to the desired resolution.  This means rewriting 
the code in photon_gen so that it does not assume there is equal flux being 
emitted from each annulus.

Of course, nothing comes for free.  If we grid too finely then heating of the
disk by photons will not be stable.

Perhaps as an intial starting step we should simply create two veesions of the
structure.  The new version could remain fixed and be improved to see how
much the heating is before going further. 

	I have an additional parameter heat to the disk structure to account
	for the fact that all photons do not carry the same energy.  I have
	created a new disk structer called qdisk to recored the heating of
	the disk when banding occurs. Right now this is simply diagnostic. 
	I initialize the qdisk structure in python.c and print the results
	after each ionization cycle to disk.

	The result for a fairly typical model was that abot 10% of the flux
	is hitting the disk.  It seems likely that in fact that in the outer
	edge of the disk, heating from the inner disk may dominate the disk
	emission itself.  In principle, one could modify the temperature
	of the disk to account for this.  Or one might speculate, that 
	this might produce emission lines in the upper atmosphere of the disk
	in low inclination systems, e.g. SS Cyg.

Unresolved errors --- 

Notes for python_47 (version 1) 2 Mar 2004 (S Sim)
-------------
-------------

Files modified since python_46 (version provided by Knox Feb 2004):
--------------

get_atomicdata.c
python.c
radiation.c
resonate.c
trans_phot.c



New files created:
-----------------

matom.c



Modifications made
------------------

get_atomicdata.c:  modifications to read in atomic data for macro atoms
		   hopefully self-explanatory
		   when reading in data the idea is that it is hierarchical
		   -i.e. if macro data is present it must be read first 
		   and other data on the same ions will be ignored. I've 
		   not checked that this is fool proof yet.

		   one issue: at present, in order to deal with recombination
		   the energy structure within a given element must be
		   self-consistent (e.g. in helium if the ground state of
		   He I is set to be energy=0 then the ground state of
		   He II MUST NOT also be =0 but must be 24eV)

python.c:	   only changed freqmin so that Pa alpha is in range.

radiation.c:	   preliminary steps have been taken for what are likely to 
		   be significant revisions of this routine in the future.
		   In the macro atom approach continua will be treated
		   as causing excitation of a macro atom at specific places
		   rather than as gradually reducing the weights of photon
		   packets. For now the reduction in weights in this routine
		   has just been switched off (by brute force) but this needs
		   my attention.

resonate.c:	   this is where the biggest changes have been made.
		   In "calculate_ds" steps have been added to treat bf
		   absorption in the macro atom method: this appears to
		   be working fine. The label "nres" is used to identify what
		   has happened: if it's greater than "nlines" then it
		   means a bf rather than bb event. No provision has yet been
		   made for ff - I intend to put this in very soon.
		   In "scatter" nres is used to identify the type of absn
		   and then "matom" is called to do the macro atom
		   calculation (note that at present the nres < nlines
		   or > nlines choice makes no difference here - making
		   the coding look silly) but in the very near future I will
		   be adding a treatment of k-packets which will require the
		   distinction to be made.
		   Note that in scatter "nres" is now a pointer.
		   In "doppler" the frequency of re-emission in a continuum
		   is dealt with in the same way as for a line.

trans_phot.c	   Very minor changes: "nres" is passed to "scatter" as a 
		   pointer.
		   "line_heat" has been removed (by brute force) so that 
		   photon bundle weights are not altered. As with 
		   changes in radiation this needs some development.


New files - details:
--------------------

matom.c:	this contains the new routines for performing macro atom
		calculations. "matom" is the routine which does the macro
		atom jumps and emission. It is called by scatter and returns
		the label "nres" for the process by which the macro atom
		de-excites. Currently it includes only radiative processes -
		collisions must be added (soon). The basic framework of this
		routine is complete (apart from collisions) but there are
		several improvements to the details that should be considered
		(see "To do" list below).
		"b21" just returns the B21 coeff - based on a12
		"alpha_sp" is used to compute the spontaneous recombination
		coefficient. It is based on similar routines (like 
		sigma_phot_topbase) but is a but different - I decided that it
		would be better to write one new code to do this rather than
		hack around with existing routines. Once the macro atom
		treatment is added fully the older stuff may not be needed.
		"alpha_sp_integrand" is just the function that is integrated
		to get "alpha_sp".



Current state of development ("done" list):
-------------------------------------------

The above modifications mean that the code does now produce emission lines 
by recombination and by non-resonant scattering. I've tested this using a very
simple (n=1 to 4) H model. I have done runs with no iteration (i.e. done a 
single ionisation cycle and looked at the spectrum in the spec_tot file) using
python_46 and the current version: as expected Pa alpha, H alpha and H beta do
appear using the new version but not using the old. Note that these are not
really recombination lines but they test the macro atom mechanism just as
effectively. Soon I will add more levels to my test hydrogen atom so that we
can start to see real recombination lines.



To do list:
----------

This list is rather long. For now I'm confining it to considerations of the
micro-physics, but obviously though is also need for the model structure in the
YSO case. No particular order:

1)
Include k-packets: this should be quite easy - it will allow ff processes and 
the thermal-pool-interaction aspect of bf to be dealt with correctly. 
Ultimately k-packets are the means by which to deal with departures from 
radiative equilibrium should we ever wish to consider that.


2)
Include ff absn. This should be very very easy to do.

3)
Include collision processes. This shouldn't be too much work. Initially
I'll try to put in bb stuff (but for completeness bf can be put in too).

4)
Escape probabilities. At present the matom routine is using approximate
escape probabilities (a routine already exists to generate these). I (we)
need to think about how to get accurate escape probabilities (assuming that
the estimated ones are not good enough). Although an issue, I see this as a
rather low priority for now - not qualitatively important.

5)
Radiation field estimators. Jbar (mean intensity in lines) and gamma
(photoionisation rate) need to be computed for the matom treatment. As a
first step I'll put in estimates for them but in the end they must be iterated.
The coding to do this might take a little time but I don't expect any 
difficulties. For different sorts of gamma are needed by the method (see 
Leon's papers - gamma, gamma_e, alpha_st, alpha_st_e). At present all these
are = 0. This is quite a high priority.

6)
Level populations. I've not thought seriously about how to do this yet, but
in the end we're going to need good populations (probably by iteration) for
all the low levels. Note that the beauty of Leon's method is that we are not
likely to need good populations for the highly excited states - even though we
are interested in modelling lines from those states - the whole macro atom
method is "driven from below".

7)
Heating and cooling rates. With the new macro atom stuff new ways of 
calculating heating and cooling rates will be needed. This will affect the
temperature calculations. I've not looked into this yet (since I've only 
been doing single iteration runs and looking at .spec_tot so far I've not
run into much about this yet).

8)
Crashes - at the time of writing the code is crashing on my machine after
producing the spec_tot file from one iteration. Don't know why this is -
memory issue? I'll look into this...


-------------------------------------------------------------------------
End 2 Mar 2004 Entry
-------------------------------------------------------------------------
-------------------------------------------------------------------------

16 Mar 2004

Minor changes made to 

1) get_atomic_data.c 
2) python.c
3) resonate.c

following suggestions of ksl based on above work.


Changes
-------
To allow code to run using old input fully introduce use of the switch
geo.rt_mode. This is set to 1 for old methods and to 2 for the new 
Macro Atom approach:

modification made:

1) python.h - comment on use of rt_mode and line_mode=6 added.
2) python.c - the way to set geo.rt_mode to 2 is via the input
	      geo.line_mode=6. This is modelled on the previous
	      geo.line_mode=5 case. If it's set =6 at input
	      then python.c will set geo.rt_mode=2 and then
              reset geo.line_mode=3. More options may be needed here later.
	      If geo.rt_mode = 2 then wind radiation is swithed off too
	      (geo.wind_radiation=0) since the macro atom method will not
	      need wind generated photons.
3) radiation.c 
	    - geo.rt_mode switch used to switch on or off the change
              in packet weights. This is not in the least bit elegant 	
	      but improvments are underway...
4) resonate.c
 	    - in calculate_ds use the geo.rt_mode switch to determine whether
	      to treat bf in the old or new ways. A couple of checks have
	      been put in to make sure nothing goes wrong (these can be
	      deleted eventually)
            - in scatter use geo.rt_mode to check if calls to matom
              are required or not.
            - in doppler a check has been added for the non-macro atom 
	      case - it is not needed but will flag if something goes
	      wrong (this should be deleted eventually).

5) trans_phot.c
	    - the removal of the call to line heat is now governed by the
	      geo.rt_mode switch.
-------------------------------------------------------------------------
End 16 Mar 2004 Entry
-------------------------------------------------------------------------
-------------------------------------------------------------------------

31 Mar 2004

Major changes have now been made to the code and all the fundamental steps for 
the implementation of the macro atom stuff have been made:

1) k-packets have now been introduced
2) all the radiation field estimators are now been iterated during the
   ionisation cycles (in much the same way as T_r)
3) the estimators are being used to compute the heating and cooling rates
   for the iteration on T_e

The most important things that are still missing and which I intend to tackle 
next are:

1) At present only bb and bf radiative processes are included. I'll work on
getting ff and collisions to work with the macro atom method.

2) Coding is needed so that a simplified but consistent treatment of minor
elements can be included.


Files modified since last set of changes (above)
------------------------------------------------

emission.c
matom.c
resonate.c
radiation.c
wind_updates2d.c



New files
---------

estimators.c


Details of changes made to files:
---------------------------------

emission.c	in "total_emission" the switch "geo.rt_mode" is used to
		divert the calculation if macro atoms are being
		used. If macro atoms are being used then the
		cooling rates are calculated by calling
		"total_fb_matoms" which computed the fb cooling rates
		using the new monte carlo estimators for stimulated
		recombination. Since the macro atom code only copes
		with fb just now that's the only cooling here. Once 
		I've added collisions and ff they'll get entries here too.


matom.c		I've corrected an error in the value of ALPHA_SP_CONSTANT
		- my fault from before. The estimators used by "matom"
		are now the 'old' values (under the new scheme the
		monte carlo estimators are iterated - the 'old' values
		are the correctly normalised values from a previous
		iteration). GAMMA IS NOT CORRECTED FOR STIMULATED 
		RECOMBINATION - this is probably not very important but
		I'll deal with it at some point.

		An important new routine "kpkt" has been added to this
		file. This routine deals with the destruction of k-packets
		- basically it is called whenever a k-packet is made and it
		decides which cooling process eliminated that k-packet and
		converts it back to a photon bundle. Since only bf processes
		are in so far this is quite simple. In due course ff and 
		collisional processes will be added as alternative means
		of killing k-packets (this shouldn't be hard).


resonate.c	I've corrected a couple of little bugs is handing the
		value of nres for continuua.

		In "calculate_ds" I've added a call to one of the new
		routines "bb_estimators_increment" to add an appropriate
		contribution to the bound-bound monte carlo estimators.
		This comes straight after calculating tau_sobolev.
		
		In "scatter" the possibility of bf processes creating 
		k-packets is allowed. The implementation closely follows
		that in Leon paper (and my own code). If a k-packet is
		made "kpkt" is called but if now a macro atom is activated
		by calling "matom" as before.


radiation.c	In "radiation" I've added a call to the new 
		"bf_estimators_increment" which adds the appropriate 
		contribution to all the bound-free estimators. This step
		is the equivalent of the old step of decreasing the
		packet weights (which was done by "radiation").

wind_updates2d.c       
		In "wind_update" a call to "mc_estimator_normalise" is
		made - this normalises the monte carlo radiation field
		estimators and shifts them into the 'old' slots for
		use in the next iteration.
		In "wind_rat_init" a block has been added to initialise
		the estimators to 0.0 before beginning the calculation.


New file - details
------------------

estimators.c	This file contains four routines:
		The first two "bf_estimators_increment" and 
		"bb_estimators_increment" add contributions to the bf
		and bb estimators as the monte carlo calculation takes
		place. As the estimators are calculated they are stored
		in the variables 'jbar', 'gamma', 'gamma_e', 'alpha_st'
		and 'alpha_st_e'.
		The third routine "mc_estimator_normalise" normalises
		the estimators at the end of the monte carlo calculation
		and shifts them into the 'old' versions of the 
		estimator variables from where they can be read for the
		next next iteration as needed.
		The fourth routine "total_fb_matoms" is used by the
		temperature update step to compute the cooling rate
		from fb processes using the monte carlo estimators for	
		stimulated recombination as well as the usual spontaneous
		emission.
		These routines should all only be called when a macro atom
		approach is being adopted.


Current state of development
----------------------------

The modification discussed above mean that (I think) the code should now
be treating all aspects of a the macro atom method correctly for bb
radiative transitions and bf radiative continuua. Free-free and collisional
processes need to be included before it's really finished, but I think that 
all the conceptual work is now done (bf are the hardest to deal with because
they involve both macro atoms and k-packets). The code runs (at least with 
my default settings) without crashing and produces a quasi-sensible output
(by this I mean that it produces a spectrum that doesn't look completely 
stupid but I've not investigated how sensible it really it yet). 


To do list
----------

1) Include ff and collisional processes - I think this is the number one
priority now.

2) Investigate convergence - I want to look at how well the temperature
(electron) is converging.

3) Generalise the implementation to cope with a simpler treatment of minority
ions that is consistent with the full macro atom treatment of one or two
"important" ions.

4) Non-isotropic scattering - I've not actually tried any runs with
non-isotropic scattering yet - this needs to be addressed eventually. This 
probably related to improvements in the calculation of the escape 
probabilities used by "matom".

5) Level populations. I've thought about this a good deal and am leaning 
towards regarding it as a relatively low priority for the moment - but once
ff and collisional processes are in then all the rates will be known and it 
will be possible (in principle) to use a matrix inversion/iteration to 
improve estimates of the level populations, at least for the low lying levels.


-------------------------------------------------------------------------
End 31 Mar 2004 Entry
-------------------------------------------------------------------------
-------------------------------------------------------------------------
Python_48 is an attempt to integrated python_46a with python_47

April 3 -- 
Copied all of Stuart's code to directory.  Here all all of the files that have
been changed.

                                       python_46a             python_47         
 ** anisowind.c               10958      10895 Mar  14 20:14   Mar  18 15:10   !Use python_46a  
 ** atomic.h                  11474      13005 Jul  14  2002   Mar  29 11:33   !Use python_47
 ** bb.c                      12244       8372 Mar  14 12:27   Mar  18 15:10   !Use python_46a
 ** emission.c                13301      13746 Mar  14 12:27   Mar  30 09:54   !Use python_47
 ** get_atomicdata.c          43988      57979 Mar  14 12:27   Mar  18 15:10   !Use python_47
 ** knigge.c                  10781      10446 Mar  14 12:27   Mar  18 15:10   !Use python_46a
 ** Makefile                   5864       5840 Mar  12 17:41   Apr   2 18:11   !Use python_47
 ** pdf.c                     25442      25403 Mar  14 20:01   Mar  18 15:10   !Use python_46a
 ** photon_gen.c              27686      27408 Mar  14 12:27   Mar  18 15:10   !Use python_46a
 ** py_ray.c                  10767      10859 Mar  14 12:27   Mar  18 15:10   !Use python_46a
 ** python.c                  39453      39567 Mar  14 12:43   Mar  18 15:10   !Changes by both of equal complexity. Used Stuart's version and readded my changes
(Note: There was an Error call when one used macroatoms.  The routine that
should have been used was Log.  I fixed this.)
 ** python.h                  23490      24913 Mar  14 12:25   Mar  29 09:37   !Changes made by both; but most are in python_47
 ** py_wind.c                 31417      31509 Mar  14 12:27   Mar  18 15:10   !Use python_46a
 ** radiation.c               13497      14097 Mar  14 12:27   Mar  26 13:52   !Use python_47.  
(Note: In incoproating this, I fixed a small error in the way bf_estimators
was declared in this routine.  It should be int and not double )
 ** random.c                   4306       4227 Mar  14 20:24   Mar  18 15:10   !Use python_46a
 ** recipes.c                  8166       8087 Mar  14 20:20   Mar  18 15:10   !Use python_46a
 ** resonate.c                18284      29966 Mar  14 12:27   Mar  31 10:08   !Use python_47
 ** roche.c                   19700      19694 Mar  14 20:35   Mar  18 15:10   !Use python_46a
 ** templates.h               17987      18214 Mar  14 20:39   Mar  31 09:08   !Remake from scratch and eliminate main's
 ** trans_phot.c              11541      11851 Mar  14 12:27   Mar  18 15:10   !Changes made by both ksl and ss. 
 ** wind2d.c                  31009      30241 Mar  14 12:27   Mar  18 15:10   !Use python_46a
 ** wind_updates2d.c          12150      13296 Mar  14 12:27   Mar  29 11:00   !Changes made by both ksl and ss, but most are in python_47


Note that python_46a contains an additional file wind_sum.c that I included in the Makefile

The program now compiles.  The warnings that remain are the same as previously in python46a
Changes in python_48a relative to python_48

04 Apr -- ksl

At present the code is set up so that it can run either the old style data inputs or
marcro-ions but not both together in a simple fashion.  I think this is because 
we have not been sophisticated enough in the way we handle the atomic data.  

What follows is some discussion of how the code currently works:


There are a limited number of interaction points between the old code and 
and the new (macro-atom) version

python.c --

Stuart identifies the fact that we are dealing with macro atom by reading in the Line_transfer
variable, as a simple way.  He thens sets

	geo.scatter_mode = 0;   //  isotropic scattering
	geo.line_mode = 3;      //  escape probablilites
	geo.rt_mode = 2;         // Identify macro atom treatment (SS)

I wonder what is the right mixture here, and what the implications of the various choices
would be.


radiation.c  -- This is the part of the code that accounts for continuum 
opacity.  In the old version of the code the photon weight is diminished 
by the continuum opacity.  There is a fair amount of bookeeeping that goes
on as well.  In the macro-atom apporach, continuum heating is treated more
like a scattering process.  

Macro -- calls bf_estimators, based on geo.rt_mode at present.  The bf_estimators code aborts on mixed input.  This is by design, but it is not obvious why this needs to be.

It is not clear to me whether any portion of this code is needed for the
macro atom approach.  radiation is called from photon2d.c 

Would one have mixed photoionization x-sections and macro atoms, or would one
want to define all photoionizations in terms of very simple macro atoms.  In
the old version of the code, the bb transitions and fb conditions are almost
independent of one another!

It's possible our long term plan could be to treat all photoionization in the 
macro atom approach.  In that case, we might want to move the call to 
bf_estimators out of radiation.


trans_phot.c  -- This is the routine that controls the movement of photon bundles
through the grid, and responds to scattering events.

Macro- -- line_heating routine is avoided since this is intended to handle my old atoms.  Note that line_heating reduces the photon weight and puts the extra into heat_lines and heat_tot.

If heat_lines and/or heat_tot are not used in the macro-atom apporach then it may be OK just
to let this run for lines that are not macro atom based.  On the other hand, if memory
serves line heating can be substantial and is important in principle. 

emission.c -- This rotine calculates the overall emission of the wind and includes the 
generation of spontaneously emitted photons in the wind, in my old approach.   At least
for now the macro atom approach turns off the generation of photons in the wind. The portion
of the code that is called is in calculating the luminosity of the wind.  


Macro atoms.  Total emission is still called by the code, but here the only thing that is
called is the free_bound emission of macro atoms.  As noted earlier, bf and fb emission were
treated in very elementary fashion in the code and so perhaps this could be excised in the 
case there are any macro atoms.

balance_sub --- Stuart made a change here.  Note that balance_sub is not actually used by 
python.  It is part of the balnce routine for testing ionizatio.


resonate.c -- This is really where most of the action is.  

There are some changes in calculate_ds that most have to do with the difference in the
way one treats bound free and free bound transitions.


In scatter, one must do a lot more work because of the possibility of frequency shifting.  
matom is called from here.


Other issues:

At present the code does not distinguish betwen macro lines and normal lines.  This would 
clearly need to be fixed if we want to mix and match. Probably could simply add an integer, 
e.g macro_info as Stuart is using for describing whether an ion is macro_ion to each 
transition.  Alternatively could assume macro lines are read first and keep track of that
number separately.  I prefer the former option, I believe.

In matom, Stuart identifies photionizations transitions that exceed nlines.  This could 
be a problem in mixing old-style lines and new style ones, unless Stuart was careful in
get_atomicdata.  More importantly, in matom need to check that it is a transition for a
macro atom.  Note that matom is called from scatter (resonate.c)  and so this check could in principle
be made prior to calling matiom.   


04April -- ksl

get_atomicdata.c, atomic.h -- Revised code to more cleanly separate information for
	macro-atoms from old inputs.  The old inputs are generally referred
	to as "simple" although this is not always the case. 

	The new code is not compatabble with Stuart's version of get_atomicdata
	in the sense that a macro-ion is defined by the keyword IonM, which must
	be read prior to reading any of the configuration information for that
	ion.  

	I have modified the stuctures for lines, configurations, and topbase to include
	a macro_info indicator.  If macro_info is 0, then it is a simple ion or
	line.  If macro_info is 1, it is a line or level for a macro-ion.  A 
	macro_ion should only have macro_configuratons.  

	As written the same is true for lines, so there should be a firm dividing line 
	betwen macro-ions and simple-ions. But I do not ultimately think this is
	what we want for this, as we think we want to migrate completely to 
	macro-ions for calculating photoionziation and recombination, whereas we
	certainly will want to use simple lines for some ions.    

	I decided that nlines_macro should be used to indicate the number of macro lines.
	nlines represents the sum of macro and simple lines

	Revised to assure that all macro-lines are read in before any simple lines,
	or more preciesely that no additional macro-lines can be read once a simple line is read in.


estimators.c -- Made several small changes, and inspected the code to see how
	it would respond to some nlte lines being "simple" lines. The largest
	change was to case bb_estimators_increment to immediately return with
	-1 if the routine was called for a simple ion.  

	NOTE to Stuart, I suspect there is something a little bit wrong with
	bb_estimators_increment.  This is discussed in the code.  I've suggested
	a change but not implemented it.

resonate.c -- Modified to avoid matoms in case where the transition was not a macro-atom
	transition.  This did not solve the jump 1000 times problem that I had in matoms. 
Changes incorporated in python_48b

I was contemplating a significant number of changes to parts of the code Stuart has written, and so
began a new minor revision python_48b

04apr -- ksl  --

get_atomicdata.c -- Revised to assure that all macro-lines are read in before any simple lines,
	or more preciesely that no additional macro-lines can be read once a simple line is read in.
matom.c -- made numerous changes which are documented in matom.  The substantive changes
	were to correct what I believe to be an error calculating one of the upperlevels for 
	a jump, but there were lots of other small changes.
resonate.c -- Modified so that photoionization transitions are represented by *nres = NLINES + ntop_phot,
	rather than something nlines that has to be calculated for each input file.
	Moved the calculation of bf opacities and the selection of the type of continuum scatter to a 
	separate routine.  This reflects a desire to keep it straightforward to identify the old
	structure of the code, and also a desire to make it simple to be able to accommodate the
	old formats of atomicdata.  
radiation.c -- Returned to the old version fo radiation.c.  This routine is almost but not quite a NOP
	in the macro approach.  I put the few bits of code that needed to be accomplished in this
	routine in the calling routine for radiation, namely photon2d.c.


Here are some notes/musings on where I think we need to be headed with python.

For the time, being we need to be able to maintain the ability of python to calculate models as was
done in Long & Knigge.  

This means we need to be able  (a) to force all line radiation to a two-level apprximation, (b) to 
accommodate the various inputs of atomic data, (c) to be able to treat photionization as an absorption
process, and (d) to allow for heating and cooling of the electrons via all the mechanisms as previously.

For the macro atom approach, we need to be able (a) to calculate macro lines in the macro-approach, (b)
to be able to add non-macro lines, (c) to accomodate a simplified macro-approach for photoionization of
high z elements,  
Changes in Python_49

28 Apr 2004, SS

Python version 49
-----------------

This version should now contain all the coding necessary to deal with the
following:

1) collisional processes in bb lines

2) ff processes (H+ and He++)

3) dealing with major and minor elements at the same time ("major" meaning 
that a full macro atom treatment is used, "minor" meaning that only a two
level approximation is used)

4) computation of heating and cooling rates (in order that the temperature
can be determined) but (a) bf processes in major and minor ions; (b) bb
collisional processes in major and minor ions; (c) ff processes in H+ and
He++

5) computation of level populations of macro atoms by inverting a matrix
of the rate equations (populations of minor ions are still taken from a
nebular approximation)

The outstanding issues are:

1) possibly need to improve treatment of anisotropic scattering

2) find a way to compute spectra without the need to compute the
entire frequency range

3) cleaning up: my object has been to get everything in place but there
are several messy things that should be tidied up for clarity.


Details - files changed and explanations of modifications
---------------------------------------------------------

emission.c:  "total_emission": there are now cooling rates for fb
	     (both major and minor elements), bb (total for both
	     major and minor) and ff included. I think that this routine
	     should have its name changed to "total_cooling" since I think
	     this makes more sense in the current context - what do you think?

	     "total_free" and "ff": I've added an "if" statement to allow
	     for the possibility of having runs with only H and therefore
	     only one contribution to ff processes.


estimators.c:
	     "bf_estimators_increment" this routine now includes a fork
	     for dealing with major and minor elements. When minor elements
	     are encountered they just cause the heating rate to be recorded. 
	     When major elements are encountered the full bf estimators
	     are recorded (the major element heating rate is obtained from
	     the estimators at a later stage). The ff heating rate is also
	     recorded here in a similar way.

	     "mc_estimator_normalise" this routine now also makes a call
	     to get the macro atom heating due to bb processes (obtained 
	     from the routines "macro_bb_heating"). I'm not convinced that
	     this is best place to do this but it keeps the Macro Atom
	     bf and bb heating together. This routine also sets the
	     quantity "geo.macro_ioniz_mode = 1": this tells the code that
	     estimators are now available for the computation of macro
	     atom level populations.

	     "total_bb_cooling" this is a new routine to compute the total
	     cooling due to bb collisions. It is used in the heating/
	     cooling calculation of Te.

	     "macro_bb_heating" this is a new routine that computes the
	     heating due to bb collisions in macro atoms. It is called
	     from "mc_estimator_normalise" at the moment.

	     "bb_minor_heat" this is a new routine that computes the
	     heating due to bb collisions in minor atoms. It is called
	     every time that a photon bundle comes into resonance with
	     a minor ion line.


levels.c     "levels" has been slightly modified so that IF macro atom
	     level populations are being computed from MC estimators this
	     routine will not change the level populations of macro atom 
	     levels in any way. Just an extra "if" statement added with 
	     a check on both "macro_info" and the new quantity
	     "macro_ioniz_mode": "levels" only changes populations when
	     both these things are equal to 0.

matom.c	     "matom" now allows for bb collisions in both jumping and
	     emission probabilities. If collisional deactivation occurs
	     then a k-packet is created and followed with "kpkt".

	     "kpkt" cooling processes due to bb collisions (in both
	     major and minor ions) bf recombinations (in both major and
	     minor ions) and ff are now included. Collisions lead to
	     an excited macro atom (and a call to "matom" or else
	     "fake_matom_bb"). bf and ff lead to an r-packet. For ff the 
	     frequency of the r-packet is determined using the routine 
	     "one_ff".

	     "fake_matom_bb" is a new routine that deals with interactions
	     of photon bundles with lines of minor elements. It is called
	     following absn by a minor line. It determined the relative
	     deactivation probabilities in r- and k-packets and then either
	     returns a scattered r-packet or converts the r-packet to a
	     k-packet (and calls "kpkt").

	     "fake_matom_bf" conceptually equivalent to "fake_matom_bb".
	     However, since there are no bf collisions at present the only
	     deactivation process is recombination - so all this routine
	     does is compute the frequency of the recombination photon.

	     "macro_pops" is a new routine supposed to compute the populations 
	     of macro atom levels from the monte carlo estimators for the 
	     rates between levels. As written it includes bb (radiation + 
	     collisions) and bf radiation. It makes use of a standard LU 
	     decomposition matrix inversion routine (from the gsl library): 
	     it's written in such a way that it'll be easy to replace this 
	     with a better matrix inverter when/if necessary. This routine is 
	     called during the wind_update part of the code, specifically 
	     whenever the ionization fractions are recomputed. It gets both 
	     the ion densities and fractional level populations for macro 
	     atoms and stores then in "ion[nion].density" and "levden" 
	     respectively. At the time of writing, I don't really trust this 
	     routine - it needs to be checked.

python.c     I've introduced a new quantity called "geo.macro_ioniz_mode"
	     this controls whether macro atom ionization (and excitation)
	     should be taken from the nebular formula used for minor ions
	     or should be computed explicitly from the mc estimators.
	     In python.c, this quantity is set to 0 (at the start of the 
	     run the estimators are not known so we can't use them to compute
	     populations with "macro_pops"). Once estimators are known the
	     code sets "geo.macro_ioniz_mode=1" and uses "macro_pops" 
	     thereafter.

radiation.c  "kappa_ff": just an "if" loop added to allow for runs in which 
	     there's only H (i.e. no He++ contribution to ff).

recomb.c     "fb": just an "if" loop added so that the contribution of
	     macro atom fb processes is not included here. This routine
	     is used in the computation of the cooling due to minor	
	     ion fb processes (macro atom fb cooling is computed elsewhere).

	     "xinteg_fb": same change as in "fb"

resonate.c   "calculate_ds" has been modified to include ff as a "scattering"
	     process in the same way as bf. Also, when a line comes into
	     resonance in "calculate_ds" the code forks: for major ions the
	     bb estimator is incremented. For minor ions the heating 
	     contribution of the line encountered is recorded by calling
	     "bb_minor_heat".

	     "select_continuum_scattering_process" has been modified to
	     allow for ff. ff scattering is recorded by nres = -2 (note
	     nres = -1 was already used for electron scattering).

	     "scatter" has been changed quite a bit to allow for both major and
	     minor elements and also for ff processes. For bb scattering it
	     now forks to call either "matom" or "fake_matom_bb" as
	     appropriate. For bf it forks for (first) major and then (second)
	     minor elements. In both cases it makes a decision about whether
	     to excite a macro atom or create a k-packet and then calls
	     either "kpkt", "matom" or "fake_matom_bf". After checking 
	     for bb and bf it now checks for ff and, if appropriate, 
	     creates a k-packet with "kpkt". For ff re-emission it is 
	     assumed that the scattering is isotropic. This routine is
	     looking a bit messy and should be tidied up for clarity (I'll
	     try to do this soon).

saha.c	     "concentrations" now follows the call to "saha" with one to
	     "macro_pops" (if geo.macro_ioniz_mode == 1) so that the
	     macro atom ionisation/excitation is obtained from the mc
	     estimators rather than the nebular approximation.

	     "lucy": this routine is called after "concentrations". When
	     geo.macro_ioniz_mode == 1 we don't want it to alter the macro 
	     atom populations so an "if" statement has been added to 
	     prevent this.


--------------------------------------------------------
Extra things: (minor changes make Apr 30 2004)
-----------

Check for population inversions added to the macro_pops routine (in matom.c)

Correction factor for stimulated emission included in the bb estimators
(put in when the estimator is normalised).

Extra case added to python.c: if geo.line_mode is input as 7 it does macro 
atoms with anisotropic scattering (if set to 6 it does macro atoms with
isotropic scattering).
    
Reorganisation of files: I think it might be a good idea to reorganise the 
routines to do with heating and cooling to make a easier to follow what is 
going on. What I have in mind is a single file which contains all the routines
used for computing the heating and cooling rates (these are currently spread 
around, some in emission.c, some in recomb.c, some in estimators.c etc.). 
Other (unused) routines in these files could be archived somewhere in case we
want them at a later date. We could also consider renaming some of the
routines/variables (e.g. I would call "total_emission" "total_cooling")
to make the context clearer. What do you think, Knox?


Additional minor changes:
-------------------------

5/05/2004

Debugging of heating/cooling rates:

collisional cooling/heating rates: these were computing the number of 
collisional cooling/heating events rather than the energy of such events. This
is now corrected.

ff heating: removed multiplication by volume (not needed).
Changes in Python_49a

Notes (SS, May 24, 2004)
------------------------

I have spent most of the last two weeks trying to sort out various things
with the ionization cycles. I've still not looked in detail and the spectrum
calculation cycles: that comes next. In this file are assorted notes and 
comments about the ionization cycles to 1) explain what the macro atom
stuff is doing for heating/cooling and 2) to flag things that may not be
sorted out fully yet.


Heating and Cooling rates
-------------------------

In the current version of python (49) there are five distinct heating/cooling
processes:

1) ff heating and cooling

This is the simplest to deal with: I've not really made any changes to how this
was handled in pre-macro atom versions of the code. In the macro atom version
the ff heating rate is computed at the same time as the bf estimators.

2) bf heating/cooling by simple ions

This is treated in more or less the same way as all bf processes in the 
pre-macro atom version of the code. The heating contribution is recorded
(beside the ff heating rate) along with the bf estimators.

3) bf heating/cooling by macro atoms

Since the macro atom jumping probabilities require detailed information about
the photoionisation / stimulated recombination rates, these rates are computed
using Monte Carlo estimators (see the routines in "estimators.c"). These
Monte Carlo estimators are also used to get the total heating rate for macro
atoms: conceptually they are not being treated any differently from the simple
ions.

4) bb collisional heating/cooling by macro atoms

In macro atoms, scattering probabilities are not used (the concept of 
division into "scattering" and "absorption" is only clear in the
context of a two-level atom anyway). Instead the true total rates of creation
and destruction of thermal energy are considered. Thus the cooling rate
in transition 1 -> 2 is (n1 C12 h nu) and the heating rate it (n2 C21 h nu).
This is as done by Leon in his 2003 paper on macro atoms.

5) bb collisional heating/cooling by simple atoms

In principle, simple bb heating/cooling can also be treated using the toal 
rates as above. However there are two problems with this: first we prefer not
to just use (n2 C21 h nu) as the heating rate since for simple ions we will
never know n2 very accurately. Secondly, if all the simple lines are treated 
this way then we can waste a lot of time destroying k-packets via bb 
collisions in a fake two-level-atom only to get a k-packet back by collisional
de-excitation of the two-level-atom. So to avoid both these issues we use 
instead a "scattering probability" approach (as discussed in KSL's notes on 
Python) for the cooling/heating of simple atom bb transitions. This approach
means that "effective" heating/cooling rates are used which avoid the need
to ever deal explicitly with the k-pkt -> excited macro atom -> k-pkt process
in simple ions, which is a big time saver. The heating rate is computed using
a Monte Carlo estimator in much the same way as in earlier versions of the 
code.

Testing: 

I've done many runs to test the heating/cooling rates computed as discussed 
above. My primary goal has been to get good agreement between the older 
versions of the code (which do not use macro atoms) and the new version when
it treats everything as simple ions. These calculations should give identical 
results since they both use two-level-atoms for everything. At the time of
writing, it seems that this is now working (I've checked the heating rates in 
detail i.e. process by process, cell by cell using a twenty level H atom
and the ixvel CV model). The simple-only and "old-style" computations seem
to agree very well. 


Bug(s):

When running 100 grid cells (10 x 10) 35 or the 38 cells converge but the 
remaining 3 do not. These 3 are all along the outer edge of the wind and I 
think the problem is actually with the density in these cells, not the 
computation of the heating/cooling rates or photon propagation. This
needs further investigation. When adiabatic cooling is included, it masks
this effect by providing an extra cooling term, however, I think there is a 
real underlying problem. For now I've left adiabatic cooling switched off, 
not because I think it is wrong but because it complicates the heating/cooling
tests and seems to inhibit convergence.



Atomic data
-----------

In the data/atomic directory there are now files for running four different 
sets of "test" atomic data. 

h20 : this is the atomic data set for a 20 level + continuum H atom. All
line for n<15 are included and the Ly, Ba and Pa series go up to n=20.
Photoionisation from all levels is included.

h20_he20  : combined data for H and He. Same 20 level + continuum H model
as above. The He model has 10 levels of He I, 10 levels of He II and the
He III continuum state. Lines and photoionisation processes are included as 
above. I've not actually had a chance to test that this set of data works yet.

h20_standard : this is the H 20 level model + atomic data for all other 
elements as taken from the usual standard39 set. Only topbase photoionisation
processes are included since the macro atom method currently only deals with
these

h20_he20_standard : this is the 20 level H model + the 20 level He model +
"standard39" data for all the other elements. My expectation is that this
data set should ultimately by the one that gets used for "real" calculations.
I've not tested that the macro atom method deals with this data properly yet.

[I also have an atomic data set "h20_test" which I've been using to 
investigate the heating/cooling due to individual atomic processes.]
Notes on Python version 49b: S Sim 14 June 2004
----------------------------------------------

The main difference between Python 49b and Python 49a (the last Readme file)
is that the new version allows the macro atom approach to be used in the
spectral synthesis of limited frequency ranges. 

The operation of this is quite straight forward: during the ionization cycles
the code records the emissivity in the frequency range of interest due to 
each excited macro atom level and destruction of k-packets in each cell. 
In the spectral cycles a switch is flipped (geo.matom_radiation) so that the
emission of r-packets due macro atoms / k-packets is accounted for using these
emissivities rather than by scattering of r-packets. In these cycles 
r-packets are absorbed (not scattered). I've tested that this works and it 
looks like it is pretty good.

There are various minor changes to the code (described in detail below), 
ordered by the file that has been changed.


Current issues:
---------------

1) Adiabatic cooling is still switched off. It can probably be switched on
again now with no problems, but I've not tried.

2) The only strange thing I've noticed in the computed spectra is that
the Lyman continuum is rather higher when the full macro atom treatment is used
than when not. This may be a bug - needs investigation.

3) I've still not tested anisotropic scattering or the use of multiple macro
atoms at the same time (e.g. H and He together).


Modifications to the files:
---------------------------

python.h:     New quantity "geo.matom_radiation" to control creation of
	      photons by deactivation macro atoms and k-packets. It should
	      be set to 0 during ionization cycles and to 1 during spectrum
	      cycles IF the macro atom method is in use.
	      "geo.f_matom" and "geo.k_kpkt" contain the specific energy
	      emitted by deactivation of macro atoms and loss of k-packets
	      respectively.
	      New structure "emiss_range" to contain the frequency range
	      for which the macro_atom / k-packet emissivities need to be
	      recorded.
	      Also reduced LPDF to 3 so that it can cope with calculations
	      with small numbers of lines (for test purposes).

atomic.h:     NLEVELS_MACRO defined as maximum number of macro atom levels
	      (needed in the storing of level emissivities). "nlevels_macro"
	      is now stored here and must be less than "NLEVELS_MACRO" or else
	      get_atomicdata will complain.


emission.c:   In adiabatic cooling - have forced return of 0 for the moment
	      (i.e. no adiabatic cooling). It may be possible to switch this
	      back on now, but at the time of writing I've not check it.

estimators.c: Some rearranging of code in mc_estimator_normalise
	      to get the call to line heating for macro atoms out of the main 
	      loop (which was wrong before). Use of macro_pops to get level
	      populations is now switched on. 
	      
extract.c:    Modification of if statement on line 260 to allow for 
	      nres > NLINES

get_atomicdata.c:
	      nlevels_macro is no longer a local variable but is in atomic.h 
	      so that the new routines that use it (for spectral synthesis)
	      can get it easily. A couple of extra (probably unnecessary) 
	      check are also added to make sure that the level data has
	      all the macro atom levels first and that the number of macro 
	      atom levels is not too big.

ionization.c: 
	      adiabatic cooling is now switched off in emission.c rather than
	      here. Comment about this added to code.

matom.c:      There have been quite a few changes here. Among the old 
	      routines, the main difference is that matom and kpkt no longer
	      interact with each other directly: they are both governed by a
	      new higher level routine called "macro_gov" which is called
	      by "scatter". Macro_gov deals with the conversion of k-packets
	      and excited macro atoms into r-packets and prevents the need
	      for matom to call kpkt (which can then call matom again etc.)
	      which leads to very complex nesting of routine calls.
	      Macro_gov also records the macro atom/k-packet emissivities
	      during the ionization cycles in a specific frequency range which
	      are used in the detailed calculation of the spectrum. During the
	      spectral synthesis step, Macro_gov throws away photons as they
	      interact with macro atoms (or become k-packets) as appropriate.

	      New routines for spectral synthesis:
	      
	      get_kpkt_f: gets the specific emissivity in the frequency
	      range for which the spectrum is required in the entire volume.

	      get_matom_f: does the same as get_kpkt_f but for macro atom 
	      levels rather than k-packets

	      photo_gen_kpkt: generates packets from k-packet destruction
	      for the spectral cycles. Closely modelled on the routine for
	      generating wind photons in earlier versions of the code. Calls
	      kpkt.

	      photo_gen_matom: generates packets from macro atom de-activations
	      for the spectral cycles. Closely modelled on the routine for
	      generating wind photons in earlier versions of the code.

	      emit_matom: called by photo_gen_matom to select the frequency 
	      for a de-activating macro atom that makes an r-packet. Closely
	      based on matom.
	      
pdf.c	      Tiny modification so that when it says
	      "pdf_gen_from_array: all y's were zero or xmin xmax out of range 
	      of array x-- returning uniform distribution " it says whether it
	      was "allzero" or not as well.

photon_gen.c:
	      Various modifications to deal with the creating of macro atom
	      deactivation / k-packet destruction photons in the spectral
	      synthesis cycles.

	      During the spectral synthesis steps it created photons based on
	      the (already recorded) emissivities due to k-packet destruction
	      and macro atom de-activation. (The variable geo.matom_radiation
	      tells it whether or not to do this). The logic closely follows
	      the earlier treatment of wind generated photons. The various
	      new routines in matom.c are called from here.

python.c:     The ranges for the spectral synthesis (which must be known during
	      the ionization cycles so that the necessary k-packet and
	      macro atom emissivities can be recorded) are stored in the
	      new variables em_rnge.fmin and em_rnge.fmax.

	      geo.matom_radiation = 0; is set for the ionization calculations
	      and then (assuming the macro atom treatment is required)	
	      geo.matom_radiation = 1 for the spectrum cycles.


radiation.c   Rather than cutting off at CR have changed this to CR/100
	       - this makes a small difference to the heating/cooling rates

resonate.c     In "calculate_ds" I've added an error message to flag when it
	       finds vol = 0. I think this should never happen but it does 
	       - I think there may be something wrong with the photon movement
	       here.
	       Also, due to the use of a push-through distance, it checks that
	       we are in the correct cell before adding contributions to 
	       bb estimators and heating rates.

	       Some simplification of "scatter" by use of the macro_gov
	       routine: matom and kpkt are no longer called directly.

trans_phot.c   Added loop to kill off photons that macro_gov kills by setting
	       weights to 0.0 during the spectrum cycles. 
	       (Pretty minor change.)

wind_updates2d.c:
	       The condition for moving densities between cells is now 
	       based on vol > 0 rather than "inwind" (as suggested by KSL)
	       Proposed copying of ne added (although probably not important
	       as discussed with KSL).
	     
Changes in Python_50

Overview of Changes

In python, the specific subroutines to calculate the velocities and densities of various
models are contained in files such as sv.c, stellar.c.  In most recent model types, which
I had created the collection of input variables were in subroutines in the same files.  But
this was not the case with a few types, e.g sv.  So I have now made this consistent.  

The variable wind.mdot was read in python.c.  This has been moved to the appropriate subroutines
that describe the individual models, since (a) not all wind models have a wind mdot, and (b)
the yso model is going to have two.  NOTE: THIS CHANGE RESULTS IN A MODIFICATION OF THE
ORDER IN WHICH INPUT VARIABLES ARE READ.  UNLESS KPAR IS MODIFIED TO ALLOW READING OF THE
.PF FILE MULTIPLE TIMES, THEN OLD .PF FILES WILL NOT WORK, UNLESS ONE MODIFIES THEM TO PLACE
WIND.MDOT AFTER WIND.RMAX.

Old order:
    rddoub ("wind.mdot(msol/yr)", &geo.wind_mdot);
    rddoub ("wind.radmax(cm)", &geo.rmax);
    rddoub ("wind.t.init", &geo.twind);

New order
    rddoub ("wind.radmax(cm)", &geo.rmax);
    rddoub ("wind.t.init", &geo.twind);
    rddoub ("wind.mdot(msol/yr)", &geo.wind_mdot);

For stellar models, I have changed the variable geo.wind_mdot to geo.stellar_wind_mdot.  This allows
for two mass loss rates.  NOTE: THIS CHANGE IMPLIES THAT SIMPLE STELLAR WIND .PF FILES WILL HAVE
TO BE MODIFIED FOR PYTHON_50.  FORTUNATELY, THERE ARE NOT MANY OF THESE. HOWEVER, NEED TO CHECK
THAT THIS MODEL TYPE STILL WORKS AT ALL.

The yso model is implemented in a very simplistic fashion.  It is completely
contained in the the initialization of the grid which establishes the
kinematic model for the wind.  This wind region as far as the program is
concerned extends from the polar axis (z) to the outer wind cone of the Knigge
wind.  The space outside the outer windcone is treated as empty (although in
some cases there will initialized there (as is true of the Knigge wind).  I
have checked more or less that the velocities appear to be correct.  For 
example the stellar wind is not rotating, while the knigge wind is.  And I
have checked that changes in mdot of the two components cause densities to
change roughly as one would expect.  



Modifications
python.h: Added a parameter geo.stellar_wind_mdot to allow for two types of winds in the same
	model.
python.c:
	Modified to allow for call for yso style wind. Moved lines that get wind mdot to
	associated routines that read in the specific subroutines in yso.c, stellar.c, and sv.c
	that read the input data for those models.  Reworked some of inputs, involving the disk, 
        to allow for combination of yso-style model. Most of these involved when you ask for
	a disk, as had previsouly assumed there would be no disk when using spherically symmetric
	models.  The way to avoid a disk is to set the disk radius to anything smaller than
	the radius of the star.   

yso.c: Routines that are specific to a yso style wind are here.  Mainly this
	calls routines from knigge.c and stellar_wind.c.  

stellar.c: Moved subroutine get_stellar_wind_params from python.c to this file.  Reworked
	some of inputs, mainly those having to do with a disk to allow for combination
	stellar + disk wind models.

corona.c: Moved subroutine get_corona_params from python.c to this file

progaa.c: Moved subroutine get_proga_params from python.c to this file

sv.c: Moved subroutine get_sv_wind_params from python.c to this file
Changes in Python50a

S Sim June 24 2004
------------------

Python 50a is an update of Python 50 that includes various correction I've made
in the last couple of weeks. 

Differences from Python 50:

bb.c - taken from my python_49b
estimators.c - taken from my python_49b
ionization.c - taken from my python_49b
matom.c - taken from my python_49b
photon2d.c -  taken from my python_49b
python.c - frequency limit changed to same as my python_49b.
resonate.c - removed the ds_calculate outside wind error message. Not resolved.
windsum.c - taken from my python_49b
wind_updates2d.c - taken from my python_49b
Changes in python_51

SS Jul 20 2004 - notes on Python 51
-----------------------------------

Following on from visit from KSL and conference poster on YSOs I'm now
beginning what will hopefully be the last set of improvements on Python
before writing something on YSOs. The modifications from Python 50a
(versions used for conference poster) are discussed below:

speed up of code
----------------

The most obvious thing to do here is replace repeated computation of 
alpha_sp with single calculation and then storage. There's now a
quantity w[n].recomb_sp which stores the spontaneous recombination rate for
each macro atom transition (and also a related w[n].recomb_sp_e). 
Using these rather than calling alpha_sp during the cycles of the code
speed it up quite a bit.

Also, during the spectrum cycles the estimators are not recomputed. Therefore
the calls to increment the estimators are now not made at this point 
(controlled by geo.do_estimators). This also speeds it up a good bit.

Together these two changes seem to speed up the code by nearly a factor of
2 on my laptop for pure macro atom calculations. (The computation of the bf 
opacity is the limiting factor now in such cases).

There is still a bit of a problem with simple atoms - for these the 
recombination coefficients are not stored but are computed which slows us down 
again. Perhaps these should be stored too - but that will need a bit more 
memory. I've not had time to try doing this yet but will try it soon.

YSO models run much faster than CV models which I suspect is because of 
different densities. There may be ways to speed up k-packet processing which
might help the CV models to go faster. I'll have a bit more of a thing/look 
into this when I get time.



reducing memory requirements
----------------------------

I've reduced the dimension of the estimator arrays to NLEVELS_MACRO rather 
than NLTE_LEVELS. To allow for H and He macro atoms at the moment, NLTE_LEVELS
is 50. This permits a 100 x 100 grid to run on my laptop (which has 750 Mb of 
RAM) - but only just! It should be fine on a desktop with ~2Gb. Memory may 
become an issue again if we want to include macro atoms for C, N and O in the
future, but for now it should not be a worry. Perhaps for CVs it may be 
possible to use a simpler H model (or no H macro atom) and then include C
instead.



adiabatic cooling
-----------------

This is now switched back on. It doesn't seem to cause any problems to the
running of the code now.


anisotropic scattering
----------------------

Several lines of code have been moved from trans_phot.c to "scatter" in 
resonate.c in order that the "thermal trapping" model for anisotropic 
scattering can be easily used for macro atom calculations. In the near future
this should become the default scattering option for macro atoms but I've not
tested it yet (or compared to the alternative anisotropic scattering 
approach).

Changes in Python_51a

Received version 51 from SS (30July).  The intent for 51a was
just to make minor changes prior to restructuring.  

O4 Jul30 -- ksl

Initially the code did not run at all on cvs, becuase of problems
with negative bf in some cases.  Fixed this by changing alpha_sp.
See notes there.

However, the code is very slow for CVs, basically impossible to use in macro
mode, with the full macro atom files that Stuart assumed.  
Even with all macro atoms turned off using atomic/standard39
it is 50% slower than python_45 on the same input files.  

Most of the time in this case is spent in radiation, and for some reason
the time spent there is significantly more than previously. There
is not much difference between the code and so it is not obvious
what the difference was. Furthermore, the number of calls to 
radiation is not significantly greater, it just takex radiation
longer to run.  It is as if the photon distributions that radiation
sees are the difference. 

The problem predates this version of the code.  It looks like it occurred
sometime between 49 and 49b.  There were considerable changes between the
two
atomic.h
emission.c
esimators.c
extract.c
get_atomicdata.c
ionization.c
levels.c
lines.c
matom.c
pdf.c
photon2d.c
ython.c
python.h
photon_gen.c
radiation.c
recomb.c
resonate.c
saha.c
templates.h
trans_phot.c
wind2d.c
wind_sum.c
wind_updates2d.c


Turning on macro atoms causes with a file that
just has the topbase photoionization x-sections pushes the time up by
another factor of 3.  The program spends much of the time doing 
integrations of the bf x-sections (as SS has noted).
Changes in python_52

Tasks:

1. Generalize the types of source models which can be used to descibe
the disk, e.g. kurucz, bb, hubeny, by creating a generalized ascii
reading routines.  Replace kucucz.c (which reads binary files) and hub.c
(which read a single file.


Notes:

The meaning of spectype in the geo structure was changed. So that the
user is not confused the input spectypes look similar, e.g mainly 1 for
bb, but these get changed as soon as they are read into the program,
so that -1 corresponds to bb, -2 to uniform, and 0, 1, 2 ... are the
spectypes read in as models.  (Note that the routine for reading in the
input models is sufficienctly smart that it recognizes when you try to
read in the same list of models multiple times.)

Routines affected:

kurucz.c and hub.c were eliminated and replaced by two other files
	get_models.c and continuum.c  (get_models is generic and contains
	routines that read a list of files and can interpolate between
	them.)

python.c:  Various changes to the inputs

photon_gen.c : Changes were made to eliminate the routines within
	kurucz.c and hub.c and replace them with routines that call the
	two routines one_continuum and emitance_continuum


Whats left to do on this.

In my specialized routines, I had ways to deal with the possibility that
one requested a temperature out of the grid. Basically I scaled the fluxes
up or down by the ratio of fluxes one would expect from blackbodies,
rather than simply using the highest or lowest temperature in the grid.

Also, there may be difficulties when one requests a photon for wavelengths
above or below the model wavelengths.  One probably should scale these
with a BB as well.

There are some repetitive bits of code in photon_gen.  This is not good
practice.  It is possible that some of the functionality should be moved
to continuum.c (i.e the choice of BB or not could be there).

The headers, and particularly the interaction, between models.h and
templates.h has become tricky.	This is because one of the routines
passses structure.  It's probably that one should do more to assure that
the theaders are not getting messed up by assuring that theire are no
duplicate variables in the headers.

Finally, it seems to me one should be able to avoid some of the the
input variables by letting spectype = -3 when have radiation from that
particular source.


Task2 -- Add the possibility of a disk with vertical extent.

Notes: The routine that checks whether the photon has hit the disk (among
other things, is walls.  It is located in phtoton2d.c.	Walls is called
from trans_phot.c, extract.c and resonate.c  The call from resonate.c
seems to be a macro-atom addition; there was no such call in python_45.
Why is this neccessary

I decided to allow the disk height to vary something to scale as r**alpha.  T
his will allow for a simple opening angle as well as slightly more 
complicated geometries.  The input variables are the fractional 
height at the outer edge of the disk, and the power law. 

The wind  calculation is fairly isolated, so it would not be too hard
to change.  For example, one could imagine disks where the scale height 
is in the middle. Basically what one needs is an analytical expression
for the zheight of the disk.  

These calculations of when a ray hits the wind are surely fairly compute 
intensive. To address this I changed the definitions associated with w->inwind 
so that cells that are completely in the wind can be identified easily.
in the wind this is not actually necessary:

There are no coordinate system dependent portions of this part of the code, 
fortunately.

Bug? -- Several of the numerical recipes routines, most notable rtsafe and
golden are givin incompatiable pointer type warnings.  This is a problem
of long standing, e.g from at least python_45, but it is still a worry.

Files modified:

python.h: Changed the variable geo.disk_absorbs to geo.disk_type.
	To keep the program running, I changed the variable names in the
	following files as well photon2d.c, plot_roche.c (Checked that
	plot_roche does still compile),  thierry.c.  At this point I just
	changed the variables and did nothing else.  Added parameters
	geo.disk_z0 and geo.disk_z1 to allow for description of vertical
	extent of disk.  I also rearranged some of geo to collect the
	various disk parameters together a bit more.

photon2d.c:  Removed the calculation of how far the photon can travel in a
	cylindrical cell to a separate routine in cylinder.c.  This makes
	photon2d.c coordinate system independent, except for the call
	itself.   Rearranged the logic of the way smax is calculated to
	avoid calculating the distance to the edge of the wind when you
	know that the cell does not contain any of these boundaries.

python.c: Modified the input variables to allow for vertically extended
	disks.	Modified what's written out in to diak.daiag.

disk.c: Created routines to calculate the distance to the disk along the
	line of sight to a photon, and to calculate the zheight of a disk.

photon_gen.c: In photo_gen_disk, made chanages so that photons
	emerge above the disk surface when the disk is vertically
	extended.

wind2d.c -- in where_in_wind, added a new possiblity for a position
	to be located inside the disk (returns -5) in this case. Began
	modification to change the way inwind worked BUT am concerned
	about edges of the wind in some of this.  Modified walls so
	that the vertical height of the disk is taken into account.
	Reworked the meaning of inwind.  The purpose of this was to know
	when there was no danger of encountering the edge of the wind,
	or any other boundary of the disk.  This allowed one to avoid
	some calculations regarding smax in trans_phot

roche.c -- Moved numerical recipes routines contained here to recipses.c

recipes.c -- Accepted numerical recipes routines from roche.c


knigge.c, sv.c -- Modified the routines so the wind starts at the edge
	of the disk.  Basically, I calculated where the footpoints of
	the streamlines are and worked from there.  I also modified the
	routines so that they would calculate the xyz velocity in all
	cases, rather than the velocity in a cylindrical corrdinate 
	system.  THESE REALLY NEED TO BE CHECKED.

corona.c -- Made the moodification so that velocities are in xyz coords,
	and added an error that stops the program for a vertically
	extended disk since one needs to think about this before using
	it again.

Whats left to do:  

Some of the logic in trans_phot assocated with moving photons is not clear to me,
in particular it's not obvious that we need both to calculate how far
the photon can travel and then check it at the end of the calculation.

For calculating disk models, that use stellar spectra one probably wants 
to calculate the gravity given the vertical extent.  This would require a 
change.

Although photons now start a the disk surface, the direction of photon
emission is still based on the assumption that the disk is flat.  This
probably does not matter too much since our prescriptions of disks to have
very large gradients with radius, but it is a shortcoming.


Task3 -- Add illumination of the disk
Notes: Basically what is done is to accumulate information about the
illumination in each ionization cycle. For option 2, the heating is added as
an extra source of energy to be reradiated in the next cycle.  For the
detailed spectrum generation phase one must (should check this) use the
heating from the last ionization cycle.  The heating is recorded in the
wavelength independent disk structure so it has to be interpolated when one is
mapping onto the wavelengthe dependent ones that are calculated.  (This might
be done better I expect.)

For an albedo of 0, one would presumably keep the photon wavelength the same
and allow it to reradiate with the angular distribution that you expect from a
photon radiated by the disk.

Routines affected:

python.h -- Added new variable geo.disk_illum to handle the various
	possibilities with repect to disk illumination.  The following
	possibilites should be implemented

	0 -- no heating effect 
	1 -- a disk with an albedo of 1 
	2 -- a disk with an albedo of 0 -- heating changes the T of the disk
	3 -- analytical approximation to the heating.  This is additive
		to the energy generated within the disk

Options 0, 2 and 3 been implemented.



python.c -- modified input variables to generalize the disk illumination
possibilites.  A yso model does not automatically used Proga's analytic
approximation anymore

trans_phot.c -- Fixed the way heating of the disk was being accumulated

disk.c -- Modified to allow disk heating to affect the temperature dist


What's left to do: Test it, and add remaining options


Task4 -- Change the coordinate system (or more properly incorporate an
additional coordinate system

Notes:
There are not too many areas of the program that are coordinate system
dependent.  They are mainly associated with initialization (wind2d)
and issues associated with where a photon is and how to interpolate
in the grid. The problem of interpolation is made slighly complicated
because some parameters, e.g, densiites, are interpolated from the
mid-points of cells, while others are interpolated from the edges
of cells.

As a first step, I am trying to isolate all of the coordinate
system dependent routines.  As a result, I have moved the
routines that calculate distances in cylindrical cells to the file
cylindrical.c

An additional issue associated with coordinate systems is going to be how to
display the results.  Py_wind takes advantage of the fact that the grid was
essentially rectangular in the yz plane.

Routines affected:


python.h: Changed variable geo_type --> coord_type

python.c: Added user input to set coordinate system type.  Hopefully
	this will not have to change further.

translate_in_wind: Added checks for coordinate systems program supports

wind2d.c:  Nearly all of this routine is coordinate system dependent.
	Started moving the bits that are independent of coordinate system
	to a new file wind.c

wind.c:  New file for receiving coordinate system independent portions
	of wind2d.c and perhaps other routines.  Added several new 
	subroutines, model_velocity which will retrun the xyz velocity
	for whatever is the current wind model.  Also, created a routine
	for calculating the velocity gradients generically.

cylindrical.c:  This is a new file, that is intended to receive routines
	that are very dependent on this specific geometry.

knigge.c, stellar.c, corona.c, sv.c, proga.c: Modified the velocity routines associated
	with each of these routines so that the xyz velocity is returned as
	opposed to the velocity at a position in the xz plane.  At present,
	these velocities are only used in setting up coordinate grids, but
	this means that one could in principle calculate the velocity at any
	point analytically. (For proga.c this is not true really since the
	data are provided in a grid.)  Removed the velocity gradient routines
	for each model.

What's left to do: Lots


General issues:

1. Having changed the definition of inwind.  I am getting some velocity
divergence errors, which may or may not be significant.

2. The size of the wind_save file is becoming unmanageable (300 Mbytes).
Changes in Python_52a

ksl

Note: This is a continuation of the notes in Readme_52.txt since some tasks
were not completed.

Tasks:

1. Generalize the types of source models which can be used to descibe
the disk, e.g. kurucz, bb, hubeny, by creating a generalized ascii
reading routines.  Replace kucucz.c (which reads binary files) and hub.c
(which read a single file.


Notes:

The meaning of spectype in the geo structure was changed. So that the
user is not confused the input spectypes look similar, e.g mainly 1 for
bb, but these get changed as soon as they are read into the program,
so that -1 corresponds to bb, -2 to uniform, and 0, 1, 2 ... are the
spectypes read in as models.  (Note that the routine for reading in the
input models is sufficienctly smart that it recognizes when you try to
read in the same list of models multiple times.)

Routines affected:

kurucz.c and hub.c were eliminated and replaced by two other files
	get_models.c and continuum.c  (get_models is generic and contains
	routines that read a list of files and can interpolate between
	them.)

python.c:  Various changes to the inputs

photon_gen.c : Changes were made to eliminate the routines within
	kurucz.c and hub.c and replace them with routines that call the
	two routines one_continuum and emitance_continuum


Whats left to do on this.

In my specialized routines, I had ways to deal with the possibility that
one requested a temperature out of the grid. Basically I scaled the fluxes
up or down by the ratio of fluxes one would expect from blackbodies,
rather than simply using the highest or lowest temperature in the grid.

Also, there may be difficulties when one requests a photon for wavelengths
above or below the model wavelengths.  One probably should scale these
with a BB as well.

There are some repetitive bits of code in photon_gen.  This is not good
practice.  It is possible that some of the functionality should be moved
to continuum.c (i.e the choice of BB or not could be there).

The headers, and particularly the interaction, between models.h and
templates.h has become tricky.	This is because one of the routines
passses structure.  It's probably that one should do more to assure that
the theaders are not getting messed up by assuring that theire are no
duplicate variables in the headers.

Finally, it seems to me one should be able to avoid some of the the
input variables by letting spectype = -3 when have radiation from that
particular source.


Task2 -- Add the possibility of a disk with vertical extent.

Notes: The routine that checks whether the photon has hit the disk (among
other things, is walls.  It is located in phtoton2d.c.	Walls is called
from trans_phot.c, extract.c and resonate.c  The call from resonate.c
seems to be a macro-atom addition; there was no such call in python_45.
Why is this neccessary

I decided to allow the disk height to vary something to scale as r**alpha.  T
his will allow for a simple opening angle as well as slightly more 
complicated geometries.  The input variables are the fractional 
height at the outer edge of the disk, and the power law. 

The wind  calculation is fairly isolated, so it would not be too hard
to change.  For example, one could imagine disks where the scale height 
is in the middle. Basically what one needs is an analytical expression
for the zheight of the disk.  

These calculations of when a ray hits the wind are surely fairly compute 
intensive. To address this I changed the definitions associated with w->inwind 
so that cells that are completely in the wind can be identified easily.
in the wind this is not actually necessary:

There are no coordinate system dependent portions of this part of the code, 
fortunately.

Bug? -- Several of the numerical recipes routines, most notable rtsafe and
golden are givin incompatiable pointer type warnings.  This is a problem
of long standing, e.g from at least python_45, but it is still a worry.

Files modified:

python.h: Changed the variable geo.disk_absorbs to geo.disk_type.
	To keep the program running, I changed the variable names in the
	following files as well photon2d.c, plot_roche.c (Checked that
	plot_roche does still compile),  thierry.c.  At this point I just
	changed the variables and did nothing else.  Added parameters
	geo.disk_z0 and geo.disk_z1 to allow for description of vertical
	extent of disk.  I also rearranged some of geo to collect the
	various disk parameters together a bit more.

photon2d.c:  Removed the calculation of how far the photon can travel in a
	cylindrical cell to a separate routine in cylinder.c.  This makes
	photon2d.c coordinate system independent, except for the call
	itself.   Rearranged the logic of the way smax is calculated to
	avoid calculating the distance to the edge of the wind when you
	know that the cell does not contain any of these boundaries.

python.c: Modified the input variables to allow for vertically extended
	disks.	Modified what's written out in to diak.daiag.

disk.c: Created routines to calculate the distance to the disk along the
	line of sight to a photon, and to calculate the zheight of a disk.

photon_gen.c: In photo_gen_disk, made chanages so that photons
	emerge above the disk surface when the disk is vertically
	extended.

wind2d.c -- in where_in_wind, added a new possiblity for a position
	to be located inside the disk (returns -5) in this case. Began
	modification to change the way inwind worked BUT am concerned
	about edges of the wind in some of this.  Modified walls so
	that the vertical height of the disk is taken into account.
	Reworked the meaning of inwind.  The purpose of this was to know
	when there was no danger of encountering the edge of the wind,
	or any other boundary of the disk.  This allowed one to avoid
	some calculations regarding smax in trans_phot

roche.c -- Moved numerical recipes routines contained here to recipses.c

recipes.c -- Accepted numerical recipes routines from roche.c


knigge.c, sv.c -- Modified the routines so the wind starts at the edge
	of the disk.  Basically, I calculated where the footpoints of
	the streamlines are and worked from there.  I also modified the
	routines so that they would calculate the xyz velocity in all
	cases, rather than the velocity in a cylindrical corrdinate 
	system.  THESE REALLY NEED TO BE CHECKED.

corona.c -- Made the moodification so that velocities are in xyz coords,
	and added an error that stops the program for a vertically
	extended disk since one needs to think about this before using
	it again.

Whats left to do:  

Some of the logic in trans_phot assocated with moving photons is not clear to me,
in particular it's not obvious that we need both to calculate how far
the photon can travel and then check it at the end of the calculation.

For calculating disk models, that use stellar spectra one probably wants 
to calculate the gravity given the vertical extent.  This would require a 
change.

Although photons now start a the disk surface, the direction of photon
emission is still based on the assumption that the disk is flat.  This
probably does not matter too much since our prescriptions of disks to have
very large gradients with radius, but it is a shortcoming.


Task3 -- Add illumination of the disk
Notes: Basically what is done is to accumulate information about the
illumination in each ionization cycle. For option 2, the heating is added as
an extra source of energy to be reradiated in the next cycle.  For the
detailed spectrum generation phase one must (should check this) use the
heating from the last ionization cycle.  The heating is recorded in the
wavelength independent disk structure so it has to be interpolated when one is
mapping onto the wavelengthe dependent ones that are calculated.  (This might
be done better I expect.)

For an albedo of 0, one would presumably keep the photon wavelength the same
and allow it to reradiate with the angular distribution that you expect from a
photon radiated by the disk.

Routines affected:

python.h -- Added new variable geo.disk_illum to handle the various
	possibilities with repect to disk illumination.  The following
	possibilites should be implemented

	0 -- no heating effect 
	1 -- a disk with an albedo of 1 
	2 -- a disk with an albedo of 0 -- heating changes the T of the disk
	3 -- analytical approximation to the heating.  This is additive
		to the energy generated within the disk

Options 0, 2 and 3 been implemented.



python.c -- modified input variables to generalize the disk illumination
possibilites.  A yso model does not automatically used Proga's analytic
approximation anymore

trans_phot.c -- Fixed the way heating of the disk was being accumulated

disk.c -- Modified to allow disk heating to affect the temperature dist


What's left to do: Test it, and add remaining options


Task4 -- Change the coordinate system (or more properly incorporate an
additional coordinate system

Notes:
There are not too many areas of the program that are coordinate system
dependent.  They are mainly associated with initialization (wind2d)
and issues associated with where a photon is and how to interpolate
in the grid. The problem of interpolation is made slighly complicated
because some parameters, e.g, densiites, are interpolated from the
mid-points of cells, while others are interpolated from the edges
of cells.

As a first step, I am trying to isolate all of the coordinate
system dependent routines.  As a result, I have moved the
routines that calculate distances in cylindrical cells to the file
cylindrical.c

An additional issue associated with coordinate systems is going to be how to
display the results.  Py_wind takes advantage of the fact that the grid was
essentially rectangular in the yz plane.

Note that while the coordinate systems have been changed, there are
some aspects of the models that do not, mainly that we assume that 
flow is bipolar and that it is azimuthally symmetric.

Routines affected:


python.h: Changed variable geo_type --> coord_type. Added a variable w.xcen to
represent the center of each cell; this makes some of the initializatio easier
when considering more than one coordinate system.   Created the infrastructure
for the wind_cones associated with rtheta geometry

python.c: Added user input to set coordinate system type.  Hopefully
	this will not have to change further.

translate_in_wind: Added checks for coordinate systems program supports

wind2d.c:  Nearly all of this routine is coordinate system dependent.
	Started moving the bits that are independent of coordinate system
	to a new file wind.c.  04aug21--Moved coordinate system setup up to
cylindrical.c

wind.c:  New file for receiving coordinate system independent portions
	of wind2d.c and perhaps other routines.  Added several new 
	subroutines, model_velocity which will retrun the xyz velocity
	for whatever is the current wind model.  Also, created a routine
	for calculating the velocity gradients generically.  Have gone 
        through this now and removed all coordinate system dependences.
	Drive routine get_random_location to create a photon in a cell
	is now here.
	

cylindrical.c:  This is a new file, that is intended to receive routines
	that are very dependent on this specific geometry.

rtheta.c:  This is a new file that contains routines that are very specific
	to an rtheta geometry.  

knigge.c, stellar.c, corona.c, sv.c, proga.c: Modified the velocity routines associated
	with each of these routines so that the xyz velocity is returned as
	opposed to the velocity at a position in the xz plane.  At present,
	these velocities are only used in setting up coordinate grids, but
	this means that one could in principle calculate the velocity at any
	point analytically. (For proga.c this is not true really since the
	data are provided in a grid.)  Removed the velocity gradient routines
	for each model.

util.c: Created a routine that can coord_frac to calculate the fractional
position of a vector in the grid for use in inerpolation

gradv.c: Modified routines to use coord_frac.  This file should eventually be
absorbed into wind.c or wind2d or something.  Also fixed an error in dvds_ave,
which may be significant!  dvds_ave was a factor of 10 too large previously.

emission.c & matom.c: Replaced portions of routines that attempt to generate
	photons randomly within a cell, with calls to get_random_location
	which does this generically for any coordinate system.

What's left to do: Lots

The parameters geo.log etc are set up in the knigge.c, sv.c etc., or at least 
are modified there. This needs to be looked into.  The parameters need to be
set to something plausible in the main program and then modified by the
individual calling routines, or checked before the gridding is done.

matom has some coordinate sytem dependencies currently, in that photon packets
are generated.  This needs to be pulled out, and attached to a more general
routine for creating a photon packet in a cell.

gradv.c has a routine get_fractional_position that needs to be replaced.

An issue which needs to be understood is how close to the edge of the grid
one can go.  It seems to defined as MDIM-3 which seems much to much.  The
array itself is defined to MDIM-1, and one can imagine that the last complete
grid cell is MDIM-2, becase you need to be interpolate to it's edge somehow.

General issues:

1. Having changed the definition of inwind.  I am getting some velocity
divergence errors, which may or may not be significant.

2. The size of the wind_save file is becoming unmanageable (300 Mbytes).
Changes in python v52b
-----------

SS August 27 04 - minor modifications from v52a to get yso finite disk 
calculations underway. At the time of writing I've not tested if 
calculations with finite disks are really working, but they are at least
running.

Modification made since 52a:
----------------------------

disk.c:  "* geo.diskrad" added to "zdisk"
	 z-components of disk planes changed to +1 in "ds_to_disk"
	 modified loop that assigns "x2" in "ds_to_disk"
	 modified choice of which distance is "smax" and which is "smin" in "ds_to_disk"
	 modification to disk_deriv just in case z1 = 0

knigge.c: To force it to work I've removed the new loop in "kn_velocity" by changing the
	  condition on line 199 to ==3 rather than ==2. This avoids the loop which is
	  causing problems (see below).

photon_gen.c: test inserted to "photon_gen_star" to ensure that stellar photons are 
	      launched with z > height of disk at stellar edge.


python.c: block removed (L654) so that finite disk is not switched off here
	  block added (L668) just in case diskrad = 0 is ever used to mean "no disk"

wind2d.c: Modified the determination of whether a cell is in the wind or not - with a finite disk there
	  were problems with cells that have some corners in the disk and some outside the wind (generally	
	  along the xy plane) - such cells do contain wind material but were being discounted because of
	  the positions of the corners.

wind.c:	  Removed "else" from where_in_wind (L94) - don't think we wanted it.

yso.c:    In "get_yso_wind_params" I changed "geo.wind_rho_min" to = 0. Is this correct?


Current problems with finite disk
---------------------------------

To get it to work I've forced out a loop in "kn_velocity" (as described under "knigge.c" above).
This loop is correcting the velocities for the finite disk thickness but it causes problems
and seems to result in velocities greater than C. The problem arises when ds_to_disk (L210 of
knigge.c) returns an infinity. There should be an easy way round this but I've not had time
to sort it out yet.
Changes in python_52e

python v52e
-----------

The significant changes made since the last Readme file (version 52b) are:

(1) finite disk - the boundary of the wind with a finite disk is now the
outer corner of the disk rather than the point where the disk crosses the 
xy-plane. This removes the problems I was having with finite disks before.

(2) macro atom level emissivities - up to version 52d these were computed by
counting emission events in the ionization cycles. However the MC noise makes
this unsuitable. Now, absorption rates are computed during the ionization 
cycles and converted to emission rates by "get_matom_f" in "matom.c". The
accuracy of this conversion is determined by the parameter "n_tries" which 
is currently set to 1000. I suspect that 3000, or (ideally) 10000 would be 
better but very time consuming. I'm investigating this at the moment.

Various specific comments on the files are given below.

SS Oct 04



emission.c: noted a potential minor bug at line 315 

estimators.c: modifications to record the energy absorption rate by both
              macro atom levels and k-packets

knigge.c: made modifications so that with finite disk the wind is cut off
          at corner of disk rather than at intersection of disk with xy
	  plane

matom.c: modifications to incorporate new calculation of macro atom 
         level emissivities. Major changes to "get_matom_f". Also
	 macro_pops now lags the populations.

wind_updates2d.c: zeroes the macro atom/k-packet absorption rates

yso.c:  modifications similar to knigge.c for finite disk cases

Changes in python_53a -- ksl

This is a consolidation that Stuart and I carried out 04nov19-21.  His
efforts culminated in python_53 whereas I was working in python_52f.

Here I have

(a) Eliminated all of the old files from previous builds
(b) Copied my newer versions of py_wind.c and py_write.c to python_53a
(c) Added my changes for py_wind.h to Stuart's version in python_53
(d) Added stdlib.h to get_models.c so that exit(0) would be defined.
(e) Fixed the adiabatic_cooling routine in emission.c.  Note that SS
    had simply commented this out.  I put it back, but created a new
    variable geo.adiabatic that controls whether adiabatic_cooling
    is ever called in wind2d or wind2d_updates.  Even if it is called
    it doesn't appear to do anything.  But now, I can calculate the
    adiabatic cooling in py_wind.
(f) Eliminated unused variables in emission.c, disk.c
(g) Cleaned up a routines that had no defined returns in wind2d.c, windsave.c
(h) Changed version number to 53.0 in python.c
(i) Fixed claim that an operation on x might be undefined in pdf.c.  It
    should have been OK previously, and it was simply a logging operation. 
(j) Fixed prototype for rtsafe so that it corresponds to the ascii standard.

Changes in python_53b --ksl

This is mainly an effort to understand differences that have appeared as a result of 
the addition of new coordinate systems:

py_wind -- Added option X to print out various parameters at a specific place in the grid
rtheta.c -- Fixed problem with defining the center of grid cells.  This got the densities of the cylindrical and
	rtheta grid to agree more closely, but did not completely fix the differences in CIV profiles for a
	fixed abundance model with no wind emission.  

py_ray -- This is a program that calculates basic parameters along a ray.  I made minor modifications to it
	to get it running again, and added a few notes about what I remember it doing.

py_ion -- This routine had not been made coordinate independent.  When I ran py_ray I found that for a rtheta
	grid, the densities were larger than a cylindrical grid.  So I fixed it, and it seemed to produce the
	correct result, as indicated by py_ray.  It also produced a very similar result when python rtheta
	was run with fixed abundances.  04nov23  ksl
