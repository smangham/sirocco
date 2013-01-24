05apr -- ksl
I would like to fully implement spherical coordinates within python.  This
is not, however, as straightfoward as spherical polar coordinates, since
it involves a dimension change.  I also want to keep in mind that Stuart
found he could not make use of the vertically extended disk option because
the grid cells did not match the boundary conditions well.  

The problem is partially expressed in the way interpolation is done.  In
2d azimuthally symmetric coordinate systems the general interpolation 
approach is to project a vector back to the xz plan and then to interpolate
in two dimension.  Finally one projects back to to rtheta space.  Our grid
then is defined in the xz plane data points

In two dimensions we have boundary grid cells all the way around the outside
The grid is only fully defined to MDIM-3 and NDIM-3 th cells.  (In
there are operations performend on the cells where there are summations over
MDIM and NDIM separately, that might be eliminated if something internal to
the grid cell were used to detemine whether a calcuation is done. I think
volume can be used for this.

Many of the interpolations are made in cells on the x-grid points and the y grid
points.  

For spherical coordinates we only have two boundaries.


There are a lot of places where MDIM is used explicitly so that data can be
printed out in a sensible way.  It's possible that these can be combined into
a single function.  A simple example of this is spectra.c, but it needs to be
dealt with in py_wind as well.  Fortunately, I appear to that collected the
display routines in py_wind in a single routine py_write!


What is currently calculated for the grid?

position -- entire grid

velocity -- entire grid

velocity gradient -- entire grid

one-d positions used for interpolation

wind_x 

wind_midx appear to be from 0 to NDIM-1 or MDIM-1, in other words eveywhere. 

Routines changed:

Initially, I tried to get a clearer picture of when MDIM was really needed in
the 2d calculations, with a bias toward isolate this to the routines rtheta and
cylindrical if at all possilble, since these are where the coordiante system
specific routines are supposed to be.

wind_updats2d -- Fixed basic update within the cells in the wind. But there
remains a section that is very 2d dependent, having to do with how to extend
some variables outside of the wind for interpolation purposes.  

wind2d.c -- Removed Stuart's calculation of whether something was inwind, and
put the same functionality in the volume calculations, since it was
essentially the same calculation.  Fixed a problem with wind_div_v that was
giving the incorect div_v in the rtheta system  

cylindrical.c, rtheta.c -- Incorporated the inwind calculation.  Note that
as far as I can tell inwind only refers to whether something is inwind or
not. python.h indicates that at one time I had intended to say whether the
point was inside the wind or outside. I believe this "feature" is not utilized
but I need to continue to check.  If I am coorect, the comments in python.h
need cleaning up.

spectra.c -- Eliminated sections which summarized positions of last scatter
in the logfile, since something complicated had to be done about it for spherical
coordinates, and it was not clear that this was being used for anything.

ionization.c -- eliminated trivial dependence on MDIM, in favor of assuming
that if the volume is non zero then the that cell is at least partially in
the wind.

util.c -- Added spherical coordinates to the interpolation routine 
coord_fraction.  coord_fraction returns iz and dz of 0 in the case of
spherical coordinates.

photon2d.c -- We have moved with time towards a system where being inwind is 
completedly determined by the volume intergrals.  Previously, we used the
vertices of the grid cells to determine this.  This means it is now possible
for the program not to "misclassify" a cell because it has a very small
volume in the wind.  It's not obvious what to do with this because very few
photons should pass through this volume, and hence any physical parameters
derived by including this in the wind will be incorrect.  On a computational
level, though this creates a problem since a photon can hit one of these
regions, and previously this was trapping the photon in an infinite loop.  So
I modified translate_in_wind to recognize this situation and translate it
through the cell.  An error is printed.

