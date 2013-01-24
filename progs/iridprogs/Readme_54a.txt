Changes in python_54a -- ksl

Stuart had made a good deal of progress in python54 in getting a more
complete set of macro atoms to work, and had combined my work on diffent
coordinate systems as well.

In python_54a -- Knox has

A. Eliminated all of the warnings that arose from -O3 compilations in python 
py_wind , essentially variables that were not necessarily defined.  
None of the warnings was really serious.  Note that some warnings still
exist in the diagnostic programs.

B. Modified the macro portion of get_atomic_data so that it followed the 
philosophy of simple ions in the sense that one should simply be able to
change the number of levels in the elements_ions sections of the inputs,
and leave the rest of the files alone. 
