These are all of the notes made to python (in Readme files) up to
and including python_46 .... namely all of the notes made prior
to the time when Stuart Sim joined the project.

History file:

July 6, 1997--Moved everything to Python 11 on MAC at home

I believe this is the latest version of Python, but that needs 
to be checked against the verstions
at STScI.  This version does include Eddington limb darkening.  
I believe I had rearranged some of
the routines within the programs.

!!Note: I looked at extract however and there it is not
obvious that the limb darkening is properly calculated !!

The immediate task is to incorporate the possibility of 
calculating eclipses which I have been working on.

To do this I made modifications to python to get the appropriate
parameters, and had python call binary_basics.  This seemed to work.

I then made a modifications to extract and to init_spec to see if
I could get it to work.  The first attempte appeared to produce
about what I expected!  Amazing

1997 Jul 12-- Although the results of the eclipse appear basically
	correct, there is some fine scale structure in the features
	whose phyisical origin is hard to understand.  It may be necessary
	to do some detailed work to understand where this originiates.
	Possibly one needs to create a frequency shift array
	

Started Python 12 at home


1997 Aug 25--  Modified to use rdpar and mac_file in the rdpar directory and
	updated so that the input and output files on the Mac are
	in the directory "models".
	
	Added a basic capability to construct spectra from Kurucz models.  I now
	need to update the spectral summary section so I can see the underlying
	continuum.  Note that it is not entirely obvious how to do this easily
	as the unscattered continuum depends on the inclination.  Also I have
	not yet doppler shifted the continuum, which needs to be done for the 
	Kurucz models.  
	
1997 Aug 26--  Worked to make the Kurucz models more robust.  In the process
	I rearranged and added to the pf file inputs to python.  Old versions
	of the pf files will no longer work.  Specifically what I added was the
	capability to easily define what type of spectrum you want to use for the 
	star and for the disk in the ionization determination and spectrum
	determination stages.  (The affected sections of the code are python.c
	and photon_gen.c.)
	
	There is some kind of bug in one_kurucz such that it very occassionally
	was obtaining a frequency way out of range.  Because of the way the 
	frequencies are used in determining t_eff and w in the individual cells,
	the INF frequency was messing up these values in radiation.c.  One bad
	frequency could mess up the program.  I fixed the effects of the bug
	by checking that the frequencies were in range and forcing them to the
	maximum or minimum value before the return from one_kurucz.  A warning
	message is printed when this happens to remind me that the bug still
	really exists.
	
	Once this problem had been fixed I was able to consistently generate
	spectra with any of the options involving Kurucz spectra.
	
	py_smooth:  I updated py_smooth to operate in the "models" directory
	on the Mac.  It is now parallel to the other programs.  I also modified
	it so that it could handle fewer than the maximum number of spectra.

1997 Aug 27 --  saha:  I updated the Saha calculations so that the 
	partition function is partly taken into account in the calculation of
	the LTE distributions.   This is somewhat better than the old situation
	where the aprtion functions were ignored entirely.  I also added fairly
	extensive comments to this routine.
1997 Aug 27 --  permitted:  I updated permitted.c to add lines from more
	ions.  Mg, Al and several others were missing. There are now about 
	3400 lines in all.  Making these changes
	necessitated a small change in "index_lines" which is one of the 
	subsidiary routines for get_atomic_data.  The size of the structure
	describing the lines in atomic.h was also increased to 3500.
1997 Aug 27 -- photon_gen & python :  I added a very crude boundary layer
	model to python.  Basically, I assumed that a boundary layer is defined
	by a temperature and a luminosity and that the boundary layer photons
	are generated uniformly from the WD surface.  This is a cheat. 
	
	Vitello & Schlossman have a different approach.  They model the boundary 
	layer as a band around the equator of the star; to get the angular 
	distribution they assume the photons are radiated from a 45 degree surface.  
	That would be harder and I am not convinced it is any better.  One could
	check the sensitivity of wind models by considering two other alternatives,
	which might be less trouble to program--one in which the radiation
	arises from the star but you limit it to a ring of some height or 
	one in which the radiation comes from an inner zone of the disk.  

1997 Aug 28 --  I was concerned that the profiles for SV inputs did not look
	as much like SV as I would like.  I was not really sure if this is a
	problem, but insofar as possible I would like to understand it.  
	
	The main problem turned was that there were sharp features in the
	spectra of CIV which I could not understand.  The reason for this
	turned out to be that I was operating python in a diagnostic mode
	with only a 10 x 10 grid for the wind lines.  There were big changes
	in the ionization parameters from grid point to grid point and
	when I increased this to 30 x 30 the effect became much less noticeable.
	(Even at 30 x 30 there are sometimes small what look like systematic
	bumps in some spectra which still seem to me hard to understand.)
	
	A problem which I have not resolved is that at low inclinations SV
	(for their standard model) get a profile with two symmetric emission
	features split by a fairly narrow absorption feature.  In python,
	at 10 deg, the emission features always seem to be asymmetric, with
	the blue side feature stronger than the weak side feature.   I
	checked that the live or die extraction produced the same thing; the
	assymetry may have been a little less but it was stil there.
	
	

!!	A useful diagnostic would be to be able to separate out which photons
	were the result of unscatterd, singel scatttered, twice scattered
	etc, photons.  (This has been done; see below.)
	
	
!!	Another interesting diagnostic would be one in which the ionization
	abundances in the wind were fixed.  This should allow a more
	direct comparison with some of Christian's models.  Presumably
	this could be done in define_wind (part of wind_2d.c)


August 29--python

	Python and "extract" were updated to allow one to use the number of
	scatters as one of the selection criteria in generating detailed spectra.
	This allows you to see the spectrum of, for example, singly scattered
	photons, doubly scattered photons, or photons scattered more than a 
	certain amount.  This modification was only carried out for "normal"
	extractions, although it would be straightforward to do it for "live
	or die" extractions as well.

	Python was updated (simply) to write most of the files to default names.
	The root names for the default files are all "python".
	Except for diagnostic purposes I was not looking at the majority of the
	files and they were cluttering up the directories.	
	
	I have run several various of the cases, i.e. various mdot and mdot_wind
	in SV93 in figures 7, 8, and 9.  I am seeing most of the features that 
	were seen by them, although not necessarily with the same value of mdot
	or mdot wind.  It looks to me as if my ionization for CIV may be somewhat 
	higher than they report.  Also I still have the problem at low inclination.
	I have not been able to produce this.

August 31

!!	It should be possible to dynamically allocate the sizes of the wind arrays
	and the spectrum arrays since I am using calloc to allocate the arrays.  
	(Changes to "wind_save" and to ".pf" files would be required and to 
	py_wind as well as python. It may also require changes to py_smooth.	
	If I have done everything correctly, these all should be straightforward,
	but this is one change where I expect one should switch to a new version
	to try, because it might not be easy to back out)

Started Python_13

Sept 13
	Corrected a problem which made it hard to increase the number of 
	extracted spectra beyond 4.  At present, the number of spectra which
	can be extracted is 10.  If one desires to increase or decrease this
	the only thing which needs to be changed is the value of NSPEC in python.h
	
	Affected routines were python.h, spectrum.c, python.c	
	
	Added the capability to fix the ionization state of every cell in the
	wind.  This was done by adding a new input variable to choose among
	the various "ionization options", adding new choices in "define_wind"
	and "wind_update" and adding a new routine "fix_concentrations" to 
	"saha.c".  See the discussion under notes.
	
	Affect files were python.h, python.c, saha.c, wind2d.c, wind_update2d.c
	
	I also made some minor changes to py_wind.c and py_ion.c to reflect
	the fact that I had changed the behavior of rdpar so that EOF in
	interactive mode exits the program immediately and to change the
	way in which the relative concentrations are calculated in py_ion.c.
	Basically, what I did was to replace a line which calculated the
	hydrogen density.  The new line assumes that the first to ions in
	the density array are HI and HII; previously I assumed I could calculate
	nh from ne by assuming it was a cosmic plasma completely ionized.
	
Sept. 21
	Fixed the normalization of spectra in the Live or Die option so that
	it, like the extract option, should be the flux at 100 pc.  This was
	done simply by calculating the fraction of 4 PI into which one was
	accepting photons for each angle.
	
	Affected files were python.h, python.c and spectrum.c
	
	Yesterday and today I carried out some tests to search for unexpected
	features in CIV profiles.  I did not find any, either in the Live
	or DIe option or extract, either with fixed or variable densities.
	There does seem to be a difference in what Live or Die produces
	vs extract for the spectrum i=10.

Sept. 28
	I am still having problems producing the line shapes seen by SV.  I
	found an error which certainly contributes to the problem, and that is
	that I was using the incorrect relative waitings in extract for continuum
	and scattered programs.  I fixed this to correspond to the weightings
	described by Christian (in his thesis equation 2.19).  Note that if
	one changes from the Eddington approximation, one must make changes
	both in the photon generation routines and in extract.  
	
	This error is likely to be the reason Christian was unable to produce
	the correct line rations in some of his routines, especially at high
	inclination angle.  It was not immediately obvious that this would solve
	the problem I have been seeing at low inclination angle.
	
	Affected files were extract.c

Started Python 14
	
Oct. 2
	Christian has been working to improve the ionization equilibrium calculation.
	He has been comparing the results of nebular_concentrations to the same
	thing in Cloudy.  Nebular_concentrations as currently constructed appears
	to produce a plasma which is too highly ionized, and so we (he) have been
	working to improve this.
	
	He finds that assuming Te=0.9 is not always a good
	assumption.  He has fit a grid of Cloudy models and found the following
	emperical relation.
		te=3.409*pow(tr,0.9292)*pow(nh,-0.01992)*pow(w,0.1317);
	I have not incorporated this into wind_updates_2d.c
		
	!!Christian, while this may be reasonble, I would prefer some
	kind of curve which produces results which yield a physical result
	in limits we think we understad, i.e. I find it unsatisfactory
	that it te!=tr in the limit where w==1.
	
	Affected file was wind_updates_2d.c
	
	Previously, the ionization relation we had been using did not properly
	take the recombination term into account. Effectively, I had assumed
	that all recombinations went to the ground state.  
	
	The new version of nebular_concentrations implements equan 11 of
	Mazzali & Lucy.
	
	n_j+1 * n_e/n_j = W (xsi+W (1-xsi)) (te/tr)**0.5 (n_j+1 n_e/n_j)_LTE
	
	where xsi is the fraction of recombinations going to the ground state.
	
	Modifications to get_atomic_data and atomic.h were required to include
	the additional information concerning recombination rates.  (At present
	these are calculated outside of python.)
	
	Modifications to saha were also required to implemented the modified
	algorithm.
	
	In summary the affect routines were python.c,atomic.h, get_atomic_data.c
	and saha.c
	
	I also updated py_wind so that one could separately print out t_e and
	t_r.  CK's changes do make a difference in the temperatures of the
	wind far from the star.
	
Oct. 5:
	I have been having some trouble on the Mac at home with kurucz models.
	What I see is a non-reapeatable problem where in the ionization phase
	saha occassional does not converge.  The problem only occurs with Kurucz
	models.  I find the w[n].ave_freq has numbers which are not recognizatble.
	My suspicion is that this has something to do with the way one_kurucz, or
	more likely some file in pdf.c is working.  I suspect it is possible to 
	try and read a section of an array which does not exist 
	(probably just by one element).  There are some
	funny checks in wind_update.c currently as a result.

	I moved python back to STScI today.  I ran into a small problem when
	I tried to run it here on the Mac.  It made a small change in "mac_files.c".
	to correct this.  The problem was that it was not doing the right thing
	if a file did not exist.

	I checked that python did work with the kurucz grid on my Sun.  It does
	apparently, though I would still be suspicious that there is an error
	which is possibly benign on the Sun.

Oct. 8:
	I believe I found the error on the Mac which was due to a problem in rand_pdf
	(in pdf.c) where I was running off the edge of an array of structures (by 1.)

	Affected routines: pdf.c

	I added doppler shifts to the calculation of the emitted spectrum for Kurucz
	models.  This required a change to the disk structure to include the velocity
	of material in the disk at radius r.

	Affected routines: python.h, photon_gen.c

Started Python 15

Nov 19-- Made small modification to account for the fact that mass is lost from
	both sides of the disk.  This will result in a factor of 2 decrease in the
	density of the wind compared to models calculated with earlier versions.  
	Christian believes that Shlossman and Vitello did not contain this error in
	their code, even though that is what the paper says.
	
Nov. 23-- Christian has concluded and I tend to agree that there really is a
	problem at low and perhaps high inclination angles.  He has created a version
	of Python which exactly duplicates his Morgana code and still does not
	see much.  So I have been doing some tests:
	
	1. I reran the standard sv inputs and compared the results in live or die
		vs extract.  As far as I can tell the results are identical.
	2. I reran the standard sv (not changing the wind calculation) but with
		a disk radius of half the old value.  There was essentially no change.
	3. I reran the standard model but set the diskrad to zero so that photons
		should be able to pass through the disk.  This produced a peak at
		the center of about the right amplitude, though it looked too narrow.
		Exactly what happened in this case is not clear though because the disk
		rad here also affect where the photons are generated.
	4. I reran the disk but made the radius of the disk equal to the radius
	   of the wind.  This produced essentially no emission peaks.
		
	Partly the question here is what is the problem.  One possibility is that
	the scattered to unscattered fraction is too low.  If this were the
	answer then one would expect that the shape of the scattered continuum would
	be correct, i.e with about equal blue and red emission wings, but there 
	just would not be enough scattered photons.  The fact that the live or
	die option gives almost exactly the same answer as the "extract" option
	suggests the problem is not due to the reweighting of photons which goes
	on in the live or die case.  It could have to do with the original
	distributions I suppose
	
	A second possibility is that somehow the photons below the disk are not
	getting out.  Maybe too many photons scattered below the disk are not
	getting through the disk plane.  That was what my second test above was
	designed to consider.  
	
	A third possibility I guess is that photons are being created in
	the wrong place perhaps too close to the center.  In that case the photons
	are being absorbed too far on the inside and so one does not get
	photons which are scattered in the outside regions of the disk.

	...
	
	I really haven't made much progress here so I begin changing the way in
	which errors and diagnositics are logged.  Now one can log conveniently
	to a diagnositic file or both the screen and the file.  This will ultimately
	affect a lot of print.f statements.  The new routine is log.c
	
	Affected routines: python.c, spectrum.c ....
	
Nov 26, 1997
	I am still trying to understand the spectrum at low inclinations.
	
	I tried setting the electron scattering cross section to very low levels;
	there was little difference.  I do notice that there are a number of
	photons which are absorbed in extract (i.e tau>20)
	and that the vertical velocity gradients are fairly small.  Thus photons 
	which resonate are really dimenished in strength and this could be one reason 
	that the scattered flux is so small.   I ran a model in which I just looked at
	the scattered fluxes.  At 10deg, the scattered fluxes are very asymmetric blue strong
	and red weak.
	
	I then checked to assure myself that the code produces the exact same
	line shapes on opposite sides of the disk by running a model with i=170 in
	extract mode.  It does!
	
	Next, I made a modification so that instead of actually calculating the
	opacity and thereby diminishing the effect of photons which go through the
	wind, I just counted photons.  This is a simply one-line change to extract.c.
	This increased the strength of the scattered contribution to the flux
	significantly, but did not change the basic fact that the profiles were
	very asymmetric with a very weak red wing.
	
	This pushed me back toward questioning whether the accounting of photons
	hitting the disk was being done correctly.  I made a simple change to
	the program which made the geometrical size of the disk small (geo.diskrad_sq)
	without affecting anything else.  A run of the program with a vanishing
	small disk gave exactly the same profiles at i=10 deg (both for scattered
	and unsscattered light).  This was not plausible.
	
	I reinspected the z_plane_intercept routine and ran the program stepping
	through it.  Despite the fact that the algorithm was correct, it was
	giving incorrect answers in some cases.  The algorithm works basically
	by defining a line between the new photon position and the old photon
	position, and then calculating the position along the line where the
	intercept is.  The problem was that the new photon position could be a
	very large distance from the origin and I was using the new photon position
	instead of the old photon position as the starting point for my line.  The
	new photon position could be so far away from the system (10**50 cm)
	that one could not properly extend the line back to the disk.  
	
	A simple fix was simply to start from the old photon position which is 
	always quite close to the disk plane when z_plane_intercept is called.
	This finally produced a profile at 10 deg inclination which was double
	peaked (even including the effects of absorption).  It is just possible
	we have solved the problem.  It is worth noting that the profile is 
	created in a rather complicated fashion.  The single scattered photons
	have a very pronounced asymmetry, with the blue wing dominating, but
	the multiply scattered photons have a rather complicated shape and if
	one looks far from line center then the red wing is stronger.
	
	Since z_plane_intercept is only used in this calculation (called by
	walls) and since there is another generic routine which calculates
	intercepts to planes, it would probably be better (and faster) to
	eliminate z_plane_intercept altoghether and put the necessary part
	of the code directly into walls.
	
Nov 27, 1997

	z_plane_intercept was deleted and walls appropriated updated.  I checked
	that this produced exactly the same distribution at i=10 deg as before.

	Routines affected--photon_2d.c, vector.c
	
Nov 28, 1997
	Combined the various kurucz related subroutines into one file kurucz.c
	and changed the names of some of the subroutines so that I could seriously
	incorporate the hubeny routines.  Their appear to be additional kurucz
	related routines in pdf.c which I have not moved yet.
	
December 15, 1997
	There was a problem discovered when Christian and I began to run some
	detailed models at angles which are inside the wind-cone.  One would
	observe fairly regular features in the blue wing of the line in the
	standard model, especially at 40deg.
	
	The problem is due to the fact that in the gridding which I have been
	using, the grid cells get very large compared to the inner disk far
	from the star.  This is a problem because the physical conditions from
	cell to cell vary considerably.  Since the optical depth in CIV is
	of order 1, there can be a jump in the effective tau as the resonant
	surface moves across the a grid boundary because as python_15 was
	written, the densities were not interpolated within the cell.
	
	I have tried two solutions to this problem and both work at some level.
	The first is to change the gridding in the vertical direction.  In the
	program as written the spacing was logarithmic and the minimum vertical
	step size was only 10**7.  Increasing this to 10**10, reduces the observed
	problem because now the steps well away from the disk are smaller.  The
	problem with this approach of course is that now I do not really have very
	good resolution close to the disk where the conditions are changing rapidly
	(but which as it turns out don't affect the CIV profile very much.
	
	The second apprach is to calculate the densities at the vertices of the
	grid.  This causes the physical conditions to change continuously across
	the boundary of the grid.  Right now, the way this is implemented is to 
	calculate the ion densities at the center of the cells, and then to
	average 4 adjacent cells to find the density at the vertices.  At present,
	I have just added another set of densities to the wind structure, which
	doubles the size of this section.   There are several problems with
	what I have done presently (although qualitatively it is a step in the
	right direction.)  Mainly this has to do with the fact that the wind has
	edges.  There also some problems due to the fact that the grid has edges.
	
	Routines affected: wind2d.c wind_update.c resonate.c

December 19, 1997
	I have revised the procedure above, removing the densities calculated
	at the vertices.  Instead, I assume that the densities calculated refer
	to the centers of the grids, and interpolated between the centers to
	get the density at a particular point.  I have dealt with the problem
	of how to interpolate the grid on the skin of the wind by copying the
	calculated densities inside the wind (at the edges) to the outside.  
	
	I have still not made allowance for the fact that there are grid cells
	which are partially in the wind, although I believe it can now be
	done simply by finding the geometric center of the wind in each grid
	cell.  (These would have to be stored in the wind pointer however).
	
	Routines affected: wind2d.c (just to remove the code earlier in the
	    week, wind_update (to cause the wind densities to be updated outside
	    of the wind itself) and resonate to do the interpolation.

February 15, 1998
	There is a new version of rdpar that should be incorporated into python.	    

Started python_17

February 19, 1998
	I have now resynched all the versions of python at home and at work.
	The new version of rdpar has also been incorporated into python.  This
	version is designed to make it easier to use the same rdpar file with
	multiple programs.  The main point however is that if you use the same
	names in multiple programs they must be accessed in the same order.
	In principle this will allow one to use one parameter file for python,
	py_wind and py_smooth.

	Note that "python", the scripting tool, has now been put on the Suns and 
	interferes with my naming convention.  So the name of the binary has now 
	been changed to "py".
	
February 20, 1998

	I have added the capability to keep track of how many photons are scattered
	back into the disk. The results are printed to the file "python.phot". 
	For the standard SV model, there are sections of the disk where the flux of
	photons impinging on the disk are 10% of the photons being emitted by the
	disk.  I wonder how much this would change a stellar atmosphere!
	
	Routines affected: trans_phot, disk_init, phot_gen_sum, python.h
	
	
	I have added the capability to keep track of the momentum absorbed by the
	wind as a result of scattering of photons.  These are stored in dmo_dt which
	is part of the wind_ptr array.  Presumably one could now calculate the 
	radiative acceleration of the wind.  vector was changed to incorporate 
	a new routine which projects from xyz to cylintrical coordinates.
	
	Corresponding changes were made in py_wind so one could display the results.
	
	I ran this a few times.  It was not entirely clear to me that the momentum
	transfer was in the outgoing direction always.  One needs to run a longer
	cycle to make sure that it really is correct.
	
	Routines affected: wind_save,wind_read,resonate,python.h,vector,py_wind
	
February 21, 1998

	I want to attack the problem of collisionally excited emission from the wind.
	As a first step, I added a routine line_summary (in py_ion.c), so that
	py_wind could be used to display the luminosity of CIV 1550 in the wind
	due to collisional excitations and assuming the wind were optically thin.
	
	The surprising result for the SV model was that CIV emissivity in the
	wind was about 10% of the total luminosity of the system and much
	greater than the energy absorbed by the wind.  Scary!!!
	
	Routines affected: py_wind, py_ion.c

February 22, 1998

	I wondered whether an alternative from the collisional calculation I made
	yesterday would be to assume that the rate was the rate calculated assuming
	LTE and spontaneous emission, e.g.
		E = N_2 A_21 h nu
	It looks like to me a21 is about 2e9
	
	The collisional rate
		E = Ne N_2 (R_gk B_k) h nu
		
	and thus to compare the rates, one must compare
	
	A_21 to Ne (R_bk B_k)  (there is an exponential in the latter that
	  should not be included in this)
	
	For an electron density of 10**14 and T of 30,000, I found the collisional
	rate would be about 100 times less than the rate calculated by spontaneous
	emission.

March 1, 1998
	I spent time over the weekend trying to (a) improve the determination
	of the amount of energy being deposited into the plasma by photoionization,
	(b) creating a way of calculating the ionization equilibrium that would
	use the energy absorbed and energy emitted by collisional lines to
	establish the temperature of the plasma.
	
	There are still inconsistences in the results at this stage.  The program
	seems to drive to very low ionization states, but the temperatures are
	high.  This must be an error!

March 6, 1998
	I have been thinking about the problems I was having.  There were two
	problems really.
	
	1. The program seemed to produce on occassion very low ionization fractions.
	This seems to be due to the extreme sensitivity of the collisional 
	energy loss to the electron density.  If one does not calculate the temperature
	exactly, you find that the plasma has completely recombined.  I believe
	I have fixed this problem by replacing my dumb search for t with one based
	on the numerical recipes routine zbrent (which finds zeroes) for you. In
	the latest standard SV model I now get reasonable ionzation fractions and
	electron temperatures which are all near 10,000 K.
	
	2. The second problem is the very large energy loss due to collisional
	effects in the optically thick calculation.  I knew this was an overestimate
	but I did not know why.  The reason is that one has to take into
	account (a) how close the continuum field is to a BB and (b) the probability
	that the emitted line photon will be reabsorbed locally.  The escape
	probability formalism helps one deal with this particular problem.  I hve
	not implemented anything having to do with this yet, but in principle I
	now know how to proceed.
	
	The next question, once on has solved the ionization balance issue is how
	to incorporate self-emission into the program.  A way to do this is as
	follows.  Calculate the emission in a cell.  For the collisinally 
	induced line emission, add a field to the pointer which keeps track of
	this policy.  We know the typical luminosity of a MC photon (the photon
	weight.) If the total emission in a cell corresponds to say 0.1 standard
	weights one could generate a random number between 0 and 1, and if the
	random number is less than 0.1 you emit a photon with W = L*.  If not,
	no photon.  If the total emission is 5.5, you can either emit 5 photons
	with weight 1.1 or 5 with weight 1 and then generate a random number
	and choose whether to emit a sixth.

March 11, 1998
	I have been working on the emission question.  It is still not done
	though I have started to install something into python.  
	
	The main issues right now are that unlike the disk and star it is
	not easy to predict how many photons will be generated by the wind.
	The current (last working version) of emission.c could not handle
	the possibility that more than NPHOT photons need to be created.
	
	Some of the other routines in python, e.g. spectrum_create have
	not been fixed since the number of photons in the photon structure will
	not always add up to nphot.  Need to seaarch on NPHOT to do this.
	
	The current version of the program will certainly not compile.
	
March 13, 1998
	I have rethought what I have been trying to do with emission 
	and it now seems to me that I ought to be able to carry out a
	modification to python which treats the wind_radiation more
	like the other forms of radiation in the disk.  It will require
	one to estimate the luminosity of the wind, as we estimate the
	luminosity of the star in photo_gen.
	
	Specifically there will need to be a routine which estimates the
	luminosity of the wind, or the specific luminosity as the case
	may be.  From the relatively luminosity of the wind and disk and
	star ... etc. one should then be able to determine the fraction
	of the photons to generate from each component.  This will allow
	one to create a full flight of photons.  The two difficulties I
	see with this approach are that (a)	one must carry out the luminosity
	calculation twice (once to get the total and then in each cell, and
	(b) since we always want to create an exact number of photons in the
	wind, then I believe the wind photon generting program will first
	need to generate a random number to decide in which cell to generate
	a photon.  Then the routine will need to generate the photon. If you
	are generating most of the photons in a few cells this may be inefficient.
	
	
	So what did I do:
	
	1. Rolled back to older versions of python.c and trans_phot.c which
	I had modified.
	
	2. Added parameter lum_wind to geometry struct.  This will store the
	specific luminosity of the wind at the appropriate time.
	
Started python_18

	March 14, 1998--Yesterday I coded in changes to emission.c and photo_gen.c
	and emission.c to implemented a scheme where self emission from the wind
	can is included.  Today, I added changes to trans_phot so that it would
	recognize the fact that some photons are created in the wind.  Here is
	a summary of how it is supposed to work.
	
	In the main routine, one is asked whether you want to include wind
	radiation in the generation of spectra. (Note:  At present, wind cooling
	is hardcoded into the calculation of t_e and hence the ionization
	structure.  It may be desirable to make this an optional approcah.)
	Assuming you wnat wind radiation to be included in the radiative
	transfer calculation it begins to affect the ionization calculation
	after the first full cycle.  (This is because I expect the initial
	calculation of the t_e is far from correct and I wanted the initial
	cycle to be dominated by emission from everything but the wind. Parameters
	have been added to the call to define_phot to tell it whether to
	include wind radiation or not, and also whether to force a recalculation
	of the "specific luminosity" of the wind.  ("Specific Luminosity" in
	this context means the luminosity between freqmin and freqmax"). If one
	wants to include the wind radiation, the number of photon bundles to
	be generated from the wind is obtained from the fractional luminosity
	of the wind (compared to the luminosity of all the components star, bb,
	disk, and wind).  
	
	The wind photons are generated in the routine photo_gen_wind.  The basic
	approach is as follows.  As part of the intialization of the wind
	for radiation which is carried out by total_emission (which in turn
	is called by total_luminosity) the specific luminosity of each wind
	element is calculated and stored in the wind element ww->lum_rad.  To
	generate a photon the first thing that is done is to use a random number
	and the luminosities of the individual wind elements to calculate in
	which cell a photon must be created.  From this and the conditions in
	this cell, the luminosities of each line are recalculated and stored
	in xcol->pow for the lines.  (We have not stored this information in
	the wind ptr to save space.  It should be trivial from the programming
	perspective to add additional radiation mechanisms by incorporating
	them in total_emission.) Knowing the power in each line (or component
	of the radiation field, we now use a random number process to find
	out which line to create for this photon bundle.  The photon is created
	with random direction and random position in the cell.  (I do ensure
	that the photon is created in the part of the cell which contains
	the wind.)  (NOTE: At present the photons are not properly shifted
	from the local rest frame which needs to be corrected.)
	
	With this approach, photons generated within the wind can be treated
	almost like other photons.  The only thing that needed to be done
	is to assue the appropriate treatment of the photons when extracting
	them in a specific direction.  (Wind photons do not suffer limb darkening).
	This is done in trans_phot by the artifice of the variable wind_phot_start
	which is defined within phton_gen.
	
	All that remains is to test it!	
	
	Routines affected by all the changes: python.c, photo_gen.c, emission.c,
	trans_phot.c, python.h
	
	March 15, 1998 -- Today I did test it! I believe it does work at this
	point.  At least, as far as I can tell it does produce plausible results
	given the assumptions that are made.  The next steps would seem to
	be the following:
	
	1.  Introduce additional cooling terms.  Right now not only are we not
	using the escape probability formalism but the only form of cooling
	is (a small number of) lines.  Thus CIV is going to be too important.
	
	2. Fix wind photons so they are produced at the appropriate Doppler shift.
	
	3. Add more collisional lines; this can help the cooling and diminish
	the effect of CIV.
	
	4. Although f_abs is close to the radiated flux it is only close within
	a factor of 2.  Improve the search for the best fit.
	
	5. Continue to think about the ionization equilibrium.  It's got to be
	the weakest part of all of this!

March 20, 1998
	I have begun to add PdV cooling to the code.
	
	I found an error in the way wind_midx and wind_midz was calculated.  
	Previously it was not really the midpoint of the cell but the edge.
	I had put wind_midx and wind_midz in to deal with the problem
	of edges in wind profiles created by a small number of cells.  So
	one should look to see if there is a reappearance of this problem.
	Note that to fix this problem I created a new routine wind_complete
	and have it called by wind_init and by wind_read....rather than
	putting the same file in the same place.
	
	Affected routines: wind_read, wind_init	
	
March, 23, 1998
	I spent time changing various pieces of the code at home over the
	weekend, but I did not have the current version of python there so
	the first task is to integrate the edits:
	
	1. Combined bb_gen.c and planck.c to make bb.c since both of these
	files dealt with bb and it was not obvious why the two had been
	separated.
	
	2. Incorporated new pdf routines into python.  These are well-described
	in pdf.c. Basically, in various stages of writing python I had
	created two ways for storing the pdfs, and in the case of bb I 
	had some routines written into planck which allowed one to look
	at a portion of the pdf.  I have eliminated the multiple ways to
	do this, and provided general purpose tools for this.  The new
	tools allow for edges in pdfs.  This was all in preparation for
	adding bremsstrahlung and recombination emission to the code. 
	I have done a limited test to see that this works.
	
	Affected routines:  pdf.c, bb.c, random.c
	
Started python_19

April 6, 1998
	A lot was done over this weekend.  Some of it may even be correct.

	1. I removed the collisional line radiation code based on G&S
	and replaced it with a poorer approximation, but one which  permits 
	one to handle all of the lines in the line list.  This allowed me to 
	include collisional dexcitation in a straightforward manner.  
	Note that I have not yet included gaunt factors which one should
	do in this approach.

	2. I replaced the code for generating ff and fb photons.  This
	was because the existing code was very slow.  The basic change
	was to use pdf_gen_from_array, instead of pdf_gen_from_func.
	This is advantageous because pdf_gen_from_func calculates the
	function a large number of times, a number which is set globally,
	while pdf_gen_from_array allows you to set exactly which points
	you want to use.  As a result you do not need to calculate the
	function (in this case several exponentials as often).

	3. Other changes were made in subroutines in emission.c to
	store data in the appropriate structures and to modify slightly
	how photon_gen_wind and other routines worked.

	Most of this work was in various parts of emission.c

	4. I began to separate the SV dependent parts of wind2d.c out
	from the more generic bits.  This is to allow us to choose
	different wind models more easily.  This work is not complete.
	All I did was move the most obviously SV specific parts off
	to a separate file sv.c.  Some more work is need to split
	up pieces of define_wind into appropriate parts.


April 10, 1998
	
	I have begun work on incorporating escape probabilities into the code.
	
	I uncovered a problem with the use of inwind.  In some cases in emission.c
	I had assumed that p[].inwind indicated whether any part of a cell was in
	the wind.  In fact, this is not the case.  Inwind refers to whether the
	specific position defined by p[].x is in the wind.  The way to know whether
	a cell has any part in the wind is to look a p[].vol.
	
	Affected routines: emission.c, py_ion.c.
	
	I have coded the changes I believe are needed to incorporate escape probabilites.
	
	(a). I changed the inputs to python so that one could turn on de-excitation and
	escape probabilities.
	
	(b)  I added w[].one_over_dvds_ave to the wind structure
	
	(c)  I calculate one_over_dvds_ave in a subroutine of define_wind in wind2d.c.  
	A Monte Carlo integration is carried out in the subroutine one_over_dvds_ave().
	This subroutine is also in wind2d.
	
	(e)  I modified trans_phot so that it would call a routine to reweight photons which
	were resonately scattered.
	
	(f) I created the routine reweight (in emission.c) which calculates the weights.
	
	(g) I added code to radiation_fraction (in emission.c) to calculate the effect of
	the escape probabilities on collisionally excited atomes
	
	What has not been done is to thoroughly debug the code or to add code which accounts
	for the differences in the escape probabiltiy which depend on the direction in which
	the photon finally escapes.

April 18, 1998
	
	Last weekend, I ran several models with various versions of the escape probabilities.
	Christian has also tried to run some models on the unix side.  He found several
	problems with a boundaray layer model and a large wavelength range which I do
	not yet understand.  The first problem was crashing in wind_update.  This problem
	seemed similar to one I thought I had fixed previously, one in which wind update
	failed as it tried to update the the wind parameters due to divide by zeros.  The
	second problem was new in which one saw some lines that looked very strong and
	narrow in the spectrum.  It was not obvious that these lines were real at all
	and suprisingly they did not appear in all spectra of the same run.	  I need to
	keep a look out for these but so far have not seen them.
	
	There was a problem which I knew existed which I believe is now fixed.  This is
	a problem in total_fb which arose because when I calculate absorption in the
	wind due to photoionization, I have only been calculating the energy which
	goes into electrons and not the energy that went directly into ionizing an atom,
	i.e. the absorbed energy was recorded has h (nu-nu_threshold).
	
	Today I modified total_fb so that it calculates only the luminosity which
	is removed from the electrons on recombination.
	
	Routines effected: total_fb (in emission.c)
	
April 25, 1998
	
	I had found some problems when I tried to calculate a model in the EUVE.  All
	the problems appear traceable to situations where the frequencies were very
	much higher than the effective temperature.  Accordingly I added fixes in several
	places to check when freqmax < freqmin.
	
	Routines effected: pieces of bb.c and emission.c
	
April 26, 1998
	Overnight, I ran a model with EUVE.  It did run properly which in this case means
	that the program did not crash.
	
	Christian had complained that it looked like the absorbed fluxes were too high in
	some models.  I took a look at this problem this morning, and found first of all that
	the absorbed flux seemed to increase in every cycle which was suspicious.  There
	was a problem in wind_rad_init.  Not everything was being reinitialized properly; this
	has been fixed.
	
	Routines affected: wind_rad_init (in wind_updates2d.c)
	
	I also worked a bit on the formatting of errors so some messages were more clear
	or better formatted.
	
May 3, 1998
	I fixed a bug in wind_updates2d.c that was causing incorrect densities to be
	written into the wind structure on one boundary.  The problem was actually due
	to the fact that I had not completely fixed an earlier bug in the same program.
	The section of the code that was affected was the part which put densities into
	regions immediately outside the wind.  The densities there are used for intepolation
	purposes.
	
	Routine affected: wind_updates2d.c
	
	Note that there is some kind of problem, or at least a suspicios behavior, affecting
	a small number of photons in that they appear in the detailed line calculation to
	get absorbed by the secondary, even though the phase is 0.5.  The problem appears
	to affect photons at low (pole on) inclinations, more than high inclinations.

Started python_20
	
May 17, 1998
	Christian had wondered how many iterations one needed to reach convergence in the
	wind, so yesterday I made some changes to wind_updates2d.c to see how much the
	t_r and t_e were changing in each windupdate.  The results were discouraging
	because it looks as if some cells have big temperature changes even after many
	iterations (at least in the option where one is recalculating the electron temperature
	based on the absorbed emission).  Looking at the absorbed fluxes it seems clear that
	this in part due to the fact that from cell to cell, the absorbed fluxes are changing
	considerably.
	
	Today, I added more tracking to see how the average temperatures are
	evolving and also to find out what cells it is happening in. 
	
	Routine affected: wind_updates2d.c
	
	Yesterday, I also made changes to kappa_ds and to emission to try to begin to
	track how many photoionizations and recombinations of H should be occuring in the
	wind.  This is not yet working, though I don't entirely understand why this is the
	case. For the recombinations the problem is most likely due to the fact that I
	have tried to do this in the photon_gen portion of the program, but there are
	in fact very few recombination photons created there.  For the absorption, it is
	likely due to the fact that there are so few hydrogen atoms that they almost
	never create any absorption.
	
Started python_21

May 24, 1998
	Last week I made changes to wind_update so that it used inwind to determine 
	whether to update a grid cell.  This is clearly incorrect since w[].inwind
	refers to the position of one vertex of the cell, rather than whether any
	portion of that cell is in the wind.  As a result, I need to go back to
	using w[].vol as the basic discriminator of whether to attempt a wind_update.
	But in that case, I need to make sure that the volumes of the two outer edges
	of the grid are zero, because one makes a linear interpolation in the velocity
	which implies one cannot work in the outer grid cell.
	
	So I corrected wind_update to go back to using vol for deciding whether to
	attempt to update the wind ionization, and I changed the portion of wind2d
	which calculates the volume to set the volume in the edge regions to zero.

May 28, 1998
	Last weekend I began to try to install a new way of trying to calculate
	the ionization equations.  Basically, the problem I wanted to begin to
	address was the apparent overionization but low temperature of much of
	the innner wind.  Ultimately, what I wanted to do was to introduce the
	possibility of calculating the ionization of H and He more exactly.  To
	do this a major reorganization of wind_update and saha was first effected.
	Basicallly what was done was to move some of the logic for calculating
	ionization equilibrium from wind_updeate2d to separate routines (all now
	in ionization.c) and then to modify some of the ionization calculations,
	e.g. the one which does the Saha calculation itself so that it could
	be restricted to a range of elements.  The point of this is to allow one
	to use an "exact" calculation for H, say, but a saha calculation for 
	other elements.  In addition, kappa_ds was modified so that one calculated
	the number of ionizations (per atom) for H and He.  Note that kappa_ds
	was moved to become part of radiation.c.  Finally, a routine num_recomb
	was added to emission.c which calculates the recombination rate for H
	and He given an electron temperature.
	
	Changes were also made to the wind stucture to add parameters for storing
	the ionizing fluxes and number of ionizing photon bundles.
	
	I ran into a big problem when I tried all of this. What happened was that
	it appeared that there were situations where the number of ionizing photons
	in a cell went to 0.  I believe, currently, this is due to the initial
	conditions chosen.  If one sets t_wind == 10000 initially, then H is partially
	neutral, sufficiently so that there are cells which are very optically thick
	in the Lyman continuum.  In this case, the photons in the interior of the
	wind are effectively shielded by regions on the the outer edges of the wind.
	Whether one actually gets convergence if the wind starts at a sufficiently
	high temperature is yet to be determined.  (Actually it looks like the answer
	is no!)

June 7, 1998
	I have made a number of changes over the weekend, many of which are simply
	diagnostic or limit some error cases.  The main change however was in saha
	where in concen_exact and neb_con_exact, I have limited the density of neutral
	H and neutral He so that tau cannot exceed TAU_MAX in transiting a single cell.
	This prevents (from a computational sense) some of the problems that were
	causing unreasonable results.  Of course, it also implies that with these
	options the program only gives reasonable results if it settles to a place
	where one is not testing these boundaries.	
	
July 17, 1998  
	While I was in Chile, Christian made a number of changes to the
        program.  He discovered several errors, but he also created at
	least one new problem.  So  the first thing I need to to while at
	Oxford is to homogenize the code. All of Christian's changes are
        included in python_22; my latest version is in python_21.

	The biggest error was in random where I was incorrectly creating
        an isotropically directed vector.  (I was creating vectors in a
        box rather than a sphere.)  This has now been incorporated.
July 20, 1998
	
	This version of the python now incorporates most of the errors that 
	Christian had found.  It does not incoporate his formulation of
	the escape probabilities because that was too computer intensive.
	Instead I now calculate the <dvds> and use that for the escape
	probabilites (rather than <1/dvds>.  This at least avoids the
	problem of a few very small values of dvds messing up the result.

	Christian's version also includes a number of other heating and cooling
	processes.  I am not yet convinced these are important to the main
	problem at hand which is why the temperatures in a portion of the
	disk are cool.

Started python_24

Oct 4, 1998
	I have not been good about keeping track of changes here for the past several
	months.  This is because I (we) have been embarked on a major diversion
	which tries to balance heating and cooling and to obtain a better understanding
	of how the radiative transfer equation relates to the code.  Some of this is contained
	in notes.tex.  There are at present major changes in the way line heating and cooling
	are calculated and in the way in which stimulated emission is included.  In addition
	there is a new test code called balance which really is designed to understand how
	the code works for a single cell.
	
	As I leave for a vacation in France, I only note that one of the problems in
	the way python tries to do the ionization balance is that we calculate the heating
	for a fixed set of conditions (as photons traverse the grid).  We then calculate
	how one has to change the t_e so that heating and cooling are matched.  In this process
	we modify t_e and recalculating the cooling, but we do not recalculate the heating.  This
	is OK sometimes, but in real life changing the te, which can change the abundances, can
	result in very different heating, especially in regions where H is recombining.  Cooling
	on the other hand is a relatively slow function of te, slower than one might expect,
	especially when the w are near 1.  (This is because in the current scheme, n2 will be
	radiation dominated when W is near one.  In the on-the-spot approximation it is
	also true that it is almost impossible to get the ion abundances to change with T.)  This
	means that I now think that a lot of the problem associated with the electron temperature
	going to very small values is due to the fact that heating rises rapidly as H recombines.
	
	[There is also a second difficulty that in neb_con it may be hard to make H recombine).
	
	Anyway, as I leave there is an option d for dumb which allows one to take steps
	in t to head toward balance (without moving the plasma conditions so far that it does
	not converge.)
	
	[Note that one can still not work at too low a temperature in t_r in balance because
	in that case one has very few ionizing photons.]
	
Started python_25 -- Oct 17

	Python_25 is simply changed so that my home version of python and
	my Sun version are identical.  There were no obvious substantive
	differences but I was afraid they were moving apart.

Oct. 27th -- Modified pdf.c somewhat to rationalize the error returns
	and then softened some of the checking in pdf_check to prevent
	it reporting so many errors.  I have some residual concerns about
	"softening" the error checking because there are pdf's produced
	that, although they are monotonic, can have two entries the same.
	This comes about due primarily to the fact that there are certain
	values which must be in the pdf, e.g are specified by jumps.  Probably
	this is not a significant error and I wanted to worry more about
	why the wind did not seem to be cooler than it was, but it is
	one that should be revisited.

	I also continued the general program of eliminating warnings produced
	by gcc, mainly defining functions as integer.

	Changed routines include pdf.c and emission.c

Oct. 28th -- I have continued to try to get python_25 to work sensibly.  The
    problem has been that I still get errors in the pdf files that cause
    the program to halt.  The current problem is in the fb process and looks
    to be due to the fact that one was not using the same frequency limits
    in one_fb and total_fb.  I have all the routines in emission.c so that
    the frequency limits are no modified here.  This may fix my pdf problem,
    but may have increased problems in balance.
    
    But there is another problem with free bound specifically and that is
    that one must be careful with the frequency limits because a cool, but
    ionized plasma will still have free bound emission.  This is more of
    a problem for balance because in fact in python the temperature limits
    are now mostly pretty large (after the fix above)
    
    Affected routines: one_fb, one_ff, and total_emission in emission.c
    
Oct 29th -- The pdf problem did appear to get fixed, at least in the
	sense that the program now runs, and produces results for t_e which
	were similar to before.
	
	py_wind was modified so that it can display the various forms of heating
	and cooling and to include logging of the files so that one can
	print out the results.
	
	get_te and the associated routine used by zbrent were moved from emission.c
	to ionization.c since it is currently called by this routine.
	
	I am very suspicious that ionization_mod_on_the_spot does not converge
	properly.  This is evident using balance.  To duplicate the problem
	use option (c) for say t_e and t_r = 30,000 K.  Heating and cooling
	are pretty well balance.  Now run m (mod_on_the_spot) and it will drop
	t_e to the minimum and the heating and cooling balance are not quite bad.
	It is clearly not converging in all cases.
	
	I also wonder whether it is appropriate to modify the abundances at the
	same time as one is modifying t_e.  The reason I wonder this is that
	the abundances do not change much and what really would happen of course
	is that the heating is also strongly affected by the abundances.  Maybe,
	it would be simpler only to modify t and keep the abundances the same
	
	Affected routines: py_wind, emission.c, ionization.c

November 27th, 1998

	It is Thanksgiving and I have finally discovered the problem which was 
	afllicting my attempt to get the line cooling to behave in a reasonable
	fashion.  The problem was that it is not true that the line emission
	is q n_2 A21 where q is the absorption fraction in the two level approximation.
	If one uses the incorrect expression then the emissivity of the plasma
	actually appears to increase as the temperature drops in a radiation dominated
	situatuation.  This is now fixed and I have moved on to other problems
	in balance and python.
	
	Added back tracking of h and he ionizations in radiation(w,p,ds). These
	had been left out when I rewrote that routine over the summer.
	
	I fixed the approximation that previously had been a simple conversion
	factor between rho and nh by adding a simple calculation of this factor
	based on the abundances.  The new variable rho2nh is now part of atomic.h.
	The calculation is carried out in get_atomic.c.  
	
	Affected routines also include wind_update2d.c, wind2d.c
	
	balance contained many parallel subroutines to calculate ion abundances
	in various ways, e.g saha, on-the-spot etc.  I added a new routine
	ion_abundances to ionization.c which uses a variable mode to determine
	the various ways to calculate ion abundances.  This allowed me to
	consolidate (eliminate) various parts of balance.c.  It also allowed
	modifications of to python which also focussed any ion abundance 
	recalculation through "ion_abundances(w,mode)".
	
	Affected routines ionization.c,balance.c,wind2d.c,wind_update2d.c
	
	I then reran python several times against the "standard model", using
	mod_on_the_spot recalculation of the energy balance to set the electron
	temperature.  The first time I used a single scattering approximation
	and the second time I used an escape probability approximation.  From
	the ionization perspective these runs were not too different from before
	in the sense that there were some portions of the wind which had an
	electron temperature of 2000 despite fairly high radiation temperatures.
	Also there remain problems in the sense that the electron temperature
	did not converge from iteration to iteration in some portions of the
	wind.  On the other hand, the line emissivity was considerably reduced
	from previously and did not seem to affect the observed spectra too
	much except in the face on case.  I consider this progress.

python_26

December 6
	I am trying to do several things with both python and balance.  For
	balance I want to be able to implement a Sobolev radiative transfer
	using more of the python code.  For python, I would like to be able
	to implement other wind laws, e.g. radial expansion. The two problems
	are related because one needs to incorporate multiple wind array elements
	into balance which might be done through "define_wind"
	
	I have eliminated the aray element w[n].rzero.  It does not appear to
	be used anywhere.  I have added an element of geo, namely log_linear.
	It is used to change the stepping of the intervals, 0=log, 1=linear.
	
	I have changed the calls to where_in_wind and where_in_grid so they
	are simple vectors.  Thus where_in_wind tells you whether the position
	is inside or outside the wind, while where_in_grid gives you an
	array element.  It seems to me these two routines could be combined!
	
December 7
	At this point I have created in balance at least subroutines which
	should allow multiple wind laws, or at least the easy installation
	of such laws.  In the process, I have removed a section of wind_update
	whicm made an estimate of t_e based on t_r, w, etc. The formula
	used there was a parameterization of CK's I have stored the formula in the
	notes on ionization on the spot but it is not currently used anywhere.
	
	Affected routines: wind_updates.c
	
December 9:
	I have implemented the simplest Sobolev radiative transfer in balance.
	It appears to get the correct answer in the sense that heating is close
	to cooling for lines.  This version does not have any radiative transfer
	in it really and is quite independent of the relevant parts of python.
	
	The only interesting thing which I (re)discovered in doing this is that
	one needs to be careful in ones use of limit lines.  If you are going
	to sum from nline_min to nline_max, it needs to be inclusive. The python
	code looks to be correct in this regard.  The next step is to try to
	use more of the python code.
	
December 22:
	A lot of undocumented work has gone on.  The most major change is that
	I discovered the primary problem with the energy calculation which
	was due to the fact that I was not including the escape probability into
	the emission calculation of the lines.  This solved a lot of problems.
	I am now changing the code so that it can also allow a spherical model.
	Along the way I have made many small changes which I need to be documenting
	better.  
	
	One of the generic problems in python right now is that logical choices
	such as whether there is a disk or not are not done in the c convention.
	It is important to remember that in temrs of a logic (e.g in an if statement),
	0 is false and non-zero is true.
	
	Changed the name geo.idisk to geo.disk_absorbs and changed the logic associaated
	with this variable so that 0 implies no disk, 1 implies there is a disk. Fixed
	various other variables eg. bl_radiation, disk_radiation, star_radiation
	so that the logical choices conformed to the c conventions
	
	Routines affected: python.c and photon2d.c (walls), python.h (geo)	photon_gen.c
	
	Replaced check on ave_freq and j in radiation.c with sane_check
	
	Routine affected: radiation.c (radiation)
	
	I made significant modifications to the io parts of python.c and
	correspondingly to geo to accomodate a spherical wind.  In addition,
	I made changes to ds_to_wind to allow for the fact that the wind
	could either be defined as a bi-cone or as a spherical wind with
	an inner and outer radius.  I also made changes to wind_define,
	especially in terms of how the gridding was done.  (Note that although
	at this point a spherical grid has the same spacings in x and z, the
	volumes of the say 5,3 cell will not be the same as the 3,5 cell because
	the coordinate system and the cells are really defined in cylintrical
	coordinates and therefore cells close to the z axis have lower volumes.
	
	Routines affected: python.c, wind2d.c

Python 28
Jan 4, 1998 --
	Spent some time trying to speed up the program.  The main routines affected
	were in lines.c or radiation.c.  The program was about a factor of
	two quicker at the end of the day than at the beginning.  The two
	techniques I used were attempt to not calculate line emission when
	one could predict that the fluxes would be very low and to store the
	results of calculations, eg two_level_atom, so that if one asked
	for the same calculation again with exactly the same variables then
	one short circuited the calculation.  It is likely that some
	additional savings could be made in this manner.  Also, storing 
	the A21 looks like it might increase the speed by 5%.  Something
	that I did not try, would be to save values of exp, but this might
	also help especially on the MAC.  Note that I need to do a
	more careful comparison between the old and new version.

October 8, 1999
        In preparation for Janet Drew's visit, I tried to confirm that python
	still was operational and that furthermore it would operate on linux as
	well as in a sun enviornment.  I have clearly not documented any changes
	to python recently.  I believe that all of the changes that
	have taken place were designed either to make the program work with
	gcc (everywhere) or to introduce and upgraded version of rdpar, namely
	kpar.  This then is mainly a statement of how little has been done.

October 13, 1999
	I have been trying to understand the behavior of ds_to_wind which 
	calculates the distance to the edge of the wind (from outside the
	wind.  There were two problems I believe neither of which is solved
	to my satisfaction.  The first, simpler, problem is that at very
        large distances adding DFUDGE to ds does not "punch through" since 
	DFUDGE is so small relative to DS that it is lost in roundoff errors.
	The second problem is that ds_to_wind includes calculates the distance
	to the outer radius of the wind, but that where_in_wind says that
	the photon is still in the wind at the outer boundary.  where_in_wind
	does not have a separate return value to say that the photon is still
	in the wind_cone but beyond the outer boundary.  What to do needs 
	researching because a photon can be eclipsed even if it is beyond
	the outer boundary.  AT present ds_to_wind is set up so that it 
	it keeps calculating ds_to_wind until it punches thourh this outer
	sphere.  I have simplified the logic for this from the previous program
	but exactly what should be done needs to be determined! ????

October 19, 1999
	Yesterday, I believe that I fixed the problem with ds_to_wind which
	was arising in part because one needs to calculate this distance both
  	from inside and outside the wind.  Also, in retrospect it is clear that
	part of the problem was due to confusion generate by the desire to allow
	both spherical and cone-shaped winds.

	Routines affected: ds_to_wind (in photon2d.c)

	Today, I believe I have straightened out some of the problems with the
	way in which python.c was interactiong with the new version of rdpar.
	In particular, the program should once again create mod.pf and name all
	files mod.pf in those cases in which one inputs conditions from the
	command line.  I also solved what appeared to be a problem with the ispy
	routines when called in those times when one is not calculating a new
	wind structure.  This problem likely crept into the code when I tried to
	incorporate spherical models.  It is not entirely clear that ispy does
	very much that is useful and it could be EXPUNGED from the code altogether.

	Routines affected: python.c, ispy.c

	The main remaining problem of which I am aware in the code (at least as
	it runs on walleye is that the code creates a very narrow component
	to lines, indipendent of inclination.  This can't really be correct. It
	looks as if this is happening to the same photon, regardless of inclination.
	These are photons in the "emitted" spectrum.  They do not appear in the
	so-called scattered spectrum.

	The problem is due to the fact that the while statement in photo_gen_wind
	which creates the photon position in the wind is failing to work properly.
	Based on a very short run, I belive this is now fixed, and the narrow line
	features have disappeared.

	Routines affected: photo_gen_wind (in emission.c)

Dec 29, 1999 -- Fixed a problem in defining the wind radius for stellar winds.  
	Also added additional spectra to be tabulated, so that I could see what the
	spectrum of the wind, the disk, and the star looked like separately.  It
	should be possible to add additional global spectra by altering MSPEC in
	python.h and by adding the appropriate code to spectra.c
	
`	Routines affected: python.h, extract.c, python.c, specta.c

Dec 31, 1999 -- Made it possible to read some of Daniel Proga's wind models into the
	program.  In the process, I created a separate routine wind_check in wind2d.c
	to check the input varariables in the wind, and updated sane_check so it
	really works.

	Routines affected: proga.c, get_proga(in python.c), sane_check in log.c

*** Python 30

Jan 1, 2000 --- It's a new day in a New Millenium and the world will be grand.

	In the old version of Python, I made a very simple approximation for calculating the
	partition function, namely I used the multiplicty of the ground state.  In addition, 
	the old version of the code can only treat a two level atom, and hence one cannot
	for example handle HeII 1640 at all.  I would like to remove these limitations,
	so that I can more closely approximate Leon Lucy's new appraches.

	To this end, I have begun to change the structure of the atomic data.  The first step
	was to create a config structure to hold the configurations of the ions.  I had previously
	used set of parallel arrays.  Note that the configurations are not used anywhere in
	python at present.

	Routines affected: atomic.h, get_atomicdata.c

Jan 2, 2000 --- Removed the code which attempted to calculate the ionization equilibrium
	from H and He differently from the heavier elements.  This was done because I
	intend to modify some of the calls in both ionization.c and saha.c and did not
	want to continue to carry the old code along.

March 7, 2000 -- Fixed problem in concentrations which caused nonsensical results if the
	temperature was too great.  In this situation, the inital estimage of ne was
	incorrect.
	Routines affected: saha.c
	
Sept 27, 2000 -- Began work on incorporating a model which would allow one to place a corona above the disk.
	WHEN I QUIT, SOMETHING WAS GOING ON BUT IT WAS CLEAR THAT THE CORONA WAS NOT GOING UP IN THE Z DIRECTION
	PROPERLY.  A new routine corona.c was added.

	Routines affected: python.c

October 6, 2000 -- Modified details of when python saves a wind file.  This was done for diagnostic purposes in
	trying to understand the corona option.  Fixed various problems with implementing the "coronal" model.  
	One "fix" was to prevent the density from dropping too low anywhere, but I did this in the initialization
	for corona rather than in the wind2d itself.  The problem manifested itself in the ionization calculations
	which failed because nh was extremely low.  Ultimately, one would prefer such checks to be made in wind2d.

	Routines affected: python.c

November 10, 2000 -- Tried to see whether balance still works.  The answer is no.

Nov 17, 2000 -- Fixed up yet more problems with corona.  It was not setting the outside edge of the wind correctly.
	Also changed some of the order in which the scattering routines are called in trans_phot.  This was done
	so that one did not scatter a photon with zero weight (after pure absorption).  Note that in the pure
	absorption mode, electrons scattering is still calculated!

	Routines affected: python.c, trans_phot.c

Nov 29, 00 -- I made modifications to get_atomicdata so that it could handle line lists involving levels in which
	the lower level was not the ground state.  I also modified python so that different sets of line files
	could be utilized.  This means a new input variable, unfortunately.  Finally, I made modifications ot
	the routine which calculates level populations so that it could handle lines which connect two excited
	states.  This still needs work from a theoretical as well as a computational perspective.

	Routines affected: python.c, get_atomicdata.c lines.c

Dec 4, 00 -- Found an error in sobolev which made excited state lines too strong.  Have fixed to properly calculate
	ions in lower and upper states and updated sobolev to use this.  The program seems to run slower
	as a result.  I eliminated calculation of line opacity when the density of the ion
	was very small.  Also ran into a problem which produced an infinite loop in generating wind photons, 
	when it tried to create a photon in cell 0.  I am surprised
	this has not been seen previously.  I updated define_wind to set the volume of this cell to zero.

        Routines affected: calculate_ds, sobolev in resonate.c, define_wind in wind2d.c

Python 33 (shad)

June 21, 01 -- I made some modifications to the corona option so that the scale height was calculated from that
	predicted by simple disk theory.  This may not be a very good idea and may need to be removed.  In the
	process I discovered that generating Kurucz models did not include the gravity term for the disk.  I have
	fixed this (or at least attempted to... though it needs to be checked (once ddd is installed on shad.  
         It probably would make since to recheck the heating of the disk, and to move the calulation of the 
	parameters of the disk to the file disk.c.   

	Routines affected: photon_gen, disk, corona

June 24, 01 -- I checked the gravity calculation and fixed it sort of. But kurucz has no gravities larger than 
       5 so it doesn't help much

Python_34 (shad)

Aug 10, 01 -- I changed the punch through distance DFUDGE from 10**6 --> 10**5.  At 10**6 it was 10% of the 
	distance of the smallest shell. It is possible that a more sophisticated version of punch through 
	should be executged which uses a fraction of the size of a grid cell to keep the number of "failures"
	of this type to a minimum.

Python_35

Aug	I added partition functions back in using values and a routine that Ivan provided from the Opacity Project

Aug 21, 01 -- Fixed a problem with the normalizations which produced differences between disk models as
	calculated by cv_hub and python.

Notes:

1. What lines are included?  The current linelist is from Verner's list of 
resonance absorption lines.  All of the transitions therefore are connected
to the ground state.  

There is a program "permitted.c" which will read the files I got from 
the Kentucky website and reformats them so that get_atomic_data can read them.
The program requires you to specify which elements one wants; currently (9/27/97)
the following elements are included:
				"H",1,2,2,
				"He",2,1,2,
				"C",6,2,6,
				"N",7,3,7,
				"O",8,3,7,
				"Ne",10,1,8,
				"Na",11,2,9,
				"Mg",12,2,8,
				"Al",13,2,8,
				"Si",14,2,7,
				"S",16,2,7,
				"Ar",18,1,9,
				"Ca",20,1,10,
				"Fe",26,2,12,
				"Ni",28,1,11
	
The above structure describes the name,z,first ion and last ion
for which lines are extracted.
				
Note the Verner does not have Ni lines in his compendium so this is
a dummy.


2. How is the Saha equation implemented?  (The Saha calculations are in
	the file saha.c)


The Saha eqn. is  (e.g. Motz, p 102)

   N_r+1 N_e / N_r =  (2*pi*m*kT)**1.5/h**3 2 Z_r+1(T)/Z_r * exp (-X/kT)
   
where N_r+1 is the number density of the r+1 ionization state
and   Z_r+1 is the partition function of that state.  

In the current version of the program (8/27/97) the partition functions
are all set to the multiplicity of the ground state.

3. Is a boundary layer included?  Yes. It is currently (8/27/97) defined
	by a temperature and a luminosity.  The radiation pattern is identical
	to that of the WD.  

4. What are the possibilities for calculating the ionization of the wind?
	
	As of 1997 Sept. 13, there are 4 options for calculating the ionization
	of the wind.
	
	The first possibility is to use the on the spot approximation, i.e.
	to use Lucy's modified SAHA equation.  In this case, one needs to
	set "wind_ionization_mode" to 0 and needs some number of "wind_cycles".
	
	The second possibility is to use LTE.  In this case, one needs to
	set "wind_ionization_mode" to 1.  If wind_cycles are set to 0, 
	then one obtains LTE for the initial temperature of the wind. If
	wind_cycles are > 0 then it will still be LTE but the temperature
	of the photon field will be calculated from the photons which penetrate
	into the shell.
	
	The third is to fix the concentrations of specific ions.  In this
	case setting wind_cycles > 0 is a waste of time.  You will always
	get the concentrations specified in fix_concentrations.  (As of
	now the most famous ions are set to 1).
	
	The fourth is to read an old ionization calculation.  Set
	"wind_ionization_mode" to 3 in this case.
	
5.	How are the spectra which python produces normalized?
		The spectra (both in the extract and live or die options)
		are supposed to be normalized to the flux one would get
		@ 100 pc.
