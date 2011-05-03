# $Id$
# Configuration file for Mac computers in W434.  It may work on other Macs
# as well but is not tested.

# Fortran 95 compiler and linker
FC=gfortran
LD=gfortran

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside


# It is sometimes helpful to turn on floating-point trapping for the 
#  standalone program, but this will not work when using the tuner.
# These are the options for debugging symbols, bounds checking & IEEE-754 
# standard floating-point arithmetic (for consistency).
DEBUG="-g -fbounds-check -mieee-fp"

# == Warnings ==
WARNINGS="-Wall -pedantic" # This enabled most compiler warnings

# == Machine specific flags ==
# Note: If the native architecture is 64 bit (most newer Mac's) then the 
# netCDF and LAPACK libraries used must be 64 bit too
ARCH="-march=native -msse3 -mfpmath=sse -fopenmp"

# == Used to promote all real's to double precision ==
DOUBLE_PRECISION="-fdefault-real-8"

# == Optimization ==
OPTIMIZE="-O2"

# == NetCDF Location ==
# There's currently no netCDF for the Macs in W343.
# Point this to the location of your copy of netCDF if you're using a 
# different computer and use the lines below that have "netCDF v4" after them.
#NETCDF="/Users/dschanen/netcdf-gfortran" # netCDF v4

# == LAPACK libraries ==
LAPACK="-llapack -lblas" # The netlib reference LAPACK/BLAS
#LAPACK="-L/usr/lib64 -llapack -L/usr/local/atlas/lib -lf77blas -lcblas -latlas" # ATLAS BLAS (faster)

# == Linking Flags ==
# Use -s to strip (no debugging); 
# Use -L<library path> -l<lib> to link in an external library
# Use -Wl,-rpath <library path> to set a search path for shared libs
#LDFLAGS="$ARCH $LAPACK -L$NETCDF/lib -Wl,-rpath $NETCDF/lib -lnetcdf -lnetcdff" # netCDF v4
LDFLAGS="$ARCH $LAPACK" # No netCDF

# == Compiler flags ==
# You will need to `make clean' if you change these.
FFLAGS="$ARCH $OPTIMIZE $DEBUG"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
#   -Dradoffline and -Dnooverlap (see BUGSrad documentation)
# You will need to `make clean' if you change these
# Use -I<include path> to set a module or header file directory
#NETCDF_INCLUDE="$NETCDF/include/" # netCDF v4
#CPPFLAGS="-I$NETCDF_INCLUDE" # netCDF v4
#CPPDEFS="-DNETCDF -Dnooverlap -Dradoffline" # netCDF v4
CPPDEFS="-Dradoffline -Dnooverlap" # MacOS X (no netCDF)

# == Static library processing ==
AR=ar
ARFLAGS=cru
RANLIB=ranlib

# == Shared library processing ==
SHARED=$FC
SHAREDFLAGS="-fPIC -shared"

# Location of 'mkmf' utility
mkmf=$dir/mkmf

# gmake command to use and options: '-j 2' enables parallel compilation
gmake="make"
