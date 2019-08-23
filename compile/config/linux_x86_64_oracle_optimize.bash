# $Id$
# Makefile definitions for Sun Studio's Fortran compiler on GNU/Linux

# Fortran 95 compiler and linker
FC=sunf95
LD=$FC

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside

# == Optimization ==
OPTIMIZE="-C -xO4 -xvector=simd -libmil -ftrap=%none"

# == Warnings ==
# This is the preferred warning level when compiling CLUBB with Sun Studio.
#   Level -w4 produces a warning message for each call to an overloaded 
#   operator, resulting in a huge number of warnings.
WARNINGS="-w3 -ansi"

# == Machine specific flags ==
# Note that when linking to sunperf (for LAPACK) you must use -dalign
# The -g option allows for viewing the source code in the analyzer program
ARCH="-g -xtarget=native -m64 -dalign"

# == Used to promote all real's to double precision ==
DOUBLE_PRECISION="-xtypemap=real:64,double:64,integer:32"

# == NetCDF Location ==
NETCDF="/usr/local/netcdf-sun64"

# == LAPACK libraries ==
#LAPACK="-xlic_lib=sunperf" # Sun performance library
#LAPACK="-L/usr/local/lib -llapack -L/usr/local/atlas/lib -lf77blas -lcblas -latlas"

# == Linking Flags ==
# Use -s to strip (no debugging); 
# -Bstatic; for static linking of libraries
# -xlibmopt; link the optimized version of libm (this seems not to work on GNU/Linux)
LDFLAGS="-L$NETCDF/lib -lnetcdf $LAPACK -xlibmopt"

# == Compiler flags ==
# You will need to `make clean' if you change these
FFLAGS="$ARCH $OPTIMIZE"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
# Define include directories. 
# Need location of include and *.mod files for the netcdf library
CPPDEFS="-DNETCDF -DNO_LAPACK_ISNAN -DCLUBB_REAL_TYPE=8"
CPPFLAGS="-M$NETCDF/include"

# == Static library processing ==
AR=ar
ARFLAGS=cru
RANLIB=ranlib

# == Shared library processing ==
SHARED=$FC
SHAREDFLAGS="-G"

# Location of 'mkmf' utility
mkmf=$dir/mkmf

# gmake command to use and options: '-j 2' enables parallel compilation
gmake="make -j1"
