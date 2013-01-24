Python_55d ended with being able to actually generate models in all three coordinate systems.
There are however a few remaining issues, which I continue to deal with associated
with this, and I have hardly touched the problems of py_wind.  Furthermore, I want
to make py_wind produce files that I can plot more easily with tecplot S/W which
I purchased, since smongo really was not good enough for many of the things
I wanted to do.


wind_updates2d -- Moved the functionality that existed there to extend the ion
densities into cells just outside the wind into the routines cylind.c, spherical.c
and rtheta.c  This extension of the density is still one of the weaker points
in python. 

py_wind.c & py_write.c -- Multiple changes to accommodate new output files for
tecplot and also to deal with multiple geometries. Also removed the capability
to create fits files, since I the main reason I had gone to different plotting
routines was because I was not happy with the plots I could produce with the
fits files.   See the .c files themselves, for mroe information on this.

python.h -- Added a new variable geo.binary_system to make it clearer how to
deal correctly with single star systems. 

python.c -- Changed inputs to specifically ask whether this was a binarry
system.

extract.c -- Made minor mod to avoid checking whether one has hit the
secondary star in a situation where the secondary does not exist. 


The biggest conceptual problem remaining with coordinate systems in python is
how to deal with vertically extended disks.  The problem is that to believe
the results, one would like to have fine resolution above the disk at all
points in the disk, but that currently none of the coordinate systems allows
one to do this.  There are three issues really.  

1. We need to be able to determine what cell a photon is in very easily.

2. We need to be able to determine how far a photon can go in a cell before it
hits the wall of the cell.  

3. We need to be able to interpolate variables within the grid.  

All of the coordinate system grids that are currently implemented make this process
very easy, since they are regular orthogonal grids.  By regular, we mean that
the cell position in the z does not depend on x.

One conceptually simple approach to deal with the vertically extended disk
problem would be to offset the z grid as a function of x, really rho.  The
advantage of this approach is that issue 1 is then easily answered.
Presumably we would also be able to find the intersection of rays fairly
straightfowarwardly (using cones that are different for each cell in the
vertical direction), and presumably some kind of normalized position would
work as well.

An alternative to this approach would be to adopt curved grid lines.  This is
possible as long as one can easily invert from real space to grid space.  For
example, suppose z=a+r**alpha.  One can fairly easily solve for r(z) just as
easily as z(r) and therfore one can create a grid system based on this.  The
issue is that once you have defined the cell walls, how do you solve for the
distance a photon can travel within the cell.  One would actaully need to see 
whether this can be done in closed form as it can for spherical systems.
