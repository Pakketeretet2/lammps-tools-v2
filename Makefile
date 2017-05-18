## Change these options to whatever works for your system/setup:
# For delaunay triangulations:
HAVE_LIB_CGAL = 1
# For more output:
VERBOSE_LIB   = 0
# For stuff that needs a LAMMPS lib:
HAVE_LIB_LAMMPS = 1
# For armadillo:
HAVE_LIB_ARMADILLO = 1
# For ICP:
HAVE_LIB_ICP = 0
# For writing/converting gsd format.
HAVE_LIB_GSD = 1
# For reading in gzipped dump files.
HAVE_BOOST_GZIP = 1

# Some routines are parallelised with OpenMP for speed.
# This enables those:
USE_OMP = 0

# Stuff

# Makes the code that does assertions and error checking
# use exceptions rather than terminate:
USE_EXCEPTIONS = 1


include Makefile.common

EXE = liblammpstools.so

include Makefile.recipes
