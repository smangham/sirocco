06oct -- ksl Python_58 is intended to be a new baseline version of
Python, replacing py57h.  It differs from python_57h only in that it
incorporates changes from python_57ib, which affect the way random
numbers are generated.

As normally run, python58 uses a MC generated that allows for a
non-uniform probability density within the intervals at which the cdf
is stored.

For reasons that are not entirely clear, this version seems to run about
as fast as python_57h where one assumes uniform density in the intervals.
(If one generates, a huge number of random numbers with the exactly the
same parameters, the new version of pdf_get_rand and pdf_get_rand_limit
is 30% slower (as in test_bb), but in python the two versions run at
about the same speed.)  Indeed pyton_57ic which uses the old formulation
is actually slower than python_57ib (and python_58).  I was not able to
understnd this; the old verson of pdf.c is called pdf_

In the process of working on this, I wrote a random number code which
does not interpolate a cdf.  This should be extremely accurate, but it
is several orders of magnitude slower than the appraoch we normally use.
The code is in rand_func.c and is reasonably self explanatory
